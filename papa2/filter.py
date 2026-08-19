"""FASTQ quality filtering and trimming.

Python port of dada2's ``filterAndTrim``, ``fastqFilter``, and
``fastqPairedFilter`` functions.  The per-read decisions reproduce R's
implementation exactly — same order of operations, same boundary
conditions, same floating-point accumulation for expected errors — so the
filtered output is byte-identical to R's (after decompression).

Reads are processed in chunks of ``n`` records (default 1e6, R's default)
using vectorized numpy operations.
"""

from __future__ import annotations

import gzip
import multiprocessing
import os
import zlib
from concurrent.futures import ProcessPoolExecutor

try:
    # Spawn-safe worker pools that do not require a __main__ guard in the
    # calling script (forking after an OpenMP region deadlocks children).
    from loky import get_reusable_executor as _loky_executor
except ImportError:
    _loky_executor = None

try:
    # python-isal: SIMD-accelerated gzip.  Decompressed content is
    # identical; compressed files are slightly larger than zlib level 6.
    from isal import igzip as _igzip
    from isal import isal_zlib as _isal_zlib
except ImportError:
    _igzip = None
    _isal_zlib = None
from typing import Optional, Sequence, Tuple, Union

import numpy as np

from .utils import is_phix as _is_phix
from .utils import seq_complexity as _seq_complexity

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# Expected-error lookup by raw quality byte: lut[b] = 10^(-(b-33)/10).
# Index 0 (used for padding) contributes exactly 0.0, so padded positions
# do not perturb the sequential left-to-right sum (matching R's C_matrixEE,
# which stops at the end of each read).
_EE_LUT = np.zeros(256, dtype=np.float64)
_EE_LUT[33:127] = 10.0 ** (-(np.arange(33, 127, dtype=np.float64) - 33.0) / 10.0)

# Non-ACGT indicator by sequence byte (case-insensitive, matching R's
# alphabetFrequency(baseOnly=TRUE)[,"other"]).
_IS_OTHER = np.ones(256, dtype=bool)
for _b in b"ACGTacgt":
    _IS_OTHER[_b] = False

_RC_TABLE = bytes.maketrans(b"ACGTRYSWKMBDHVNacgtryswkmbdhvn",
                            b"TGCAYRSWMKVHDBNtgcayrswmkvhdbn")


def _rc_bytes(seq: bytes) -> bytes:
    return seq.translate(_RC_TABLE)[::-1]


def _open_fq(path: str, mode: str = "rb"):
    """Open a FASTQ file for binary reading, transparently handling gzip."""
    if path.endswith(".gz") or path.endswith(".gzip"):
        if _igzip is not None:
            return _igzip.open(path, mode)
        return gzip.open(path, mode)
    return open(path, mode)


class _GzipRecordWriter:
    """Buffered gzip (or plain) writer for assembled FASTQ bytes."""

    def __init__(self, path: str, compress: bool):
        if compress and not (path.endswith(".gz") or path.endswith(".gzip")):
            path += ".gz"
        self.path = path
        out_dir = os.path.dirname(path)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
        self._fh = open(path, "wb")
        if compress:
            if _isal_zlib is not None:
                self._gz = _isal_zlib.compressobj(
                    3, _isal_zlib.DEFLATED, 16 + zlib.MAX_WBITS)
            else:
                self._gz = zlib.compressobj(6, zlib.DEFLATED,
                                            16 + zlib.MAX_WBITS)
        else:
            self._gz = None

    def write(self, data: bytes):
        if self._gz is not None:
            self._fh.write(self._gz.compress(data))
        else:
            self._fh.write(data)

    def close(self):
        if self._gz is not None:
            self._fh.write(self._gz.flush())
        self._fh.close()

    def remove(self):
        """Delete the (closed) output file — R does this when no reads pass."""
        if os.path.exists(self.path):
            os.remove(self.path)


class _FastqChunkReader:
    """Yield chunks of complete FASTQ records as (buffer, line_starts, line_ends).

    ``line_starts``/``line_ends`` are int64 arrays over all lines in the
    chunk (4 per record), with CR stripped from line ends.
    """

    _READ_SIZE = 1 << 24  # 16 MB of decompressed data per read() call

    def __init__(self, path: str, chunk_records: int):
        self._fh = _open_fq(path)
        self._chunk = chunk_records
        self._leftover = b""
        self._eof = False

    def __iter__(self):
        return self

    def __next__(self):
        target_lines = self._chunk * 4
        parts = [self._leftover]
        n_newlines = self._leftover.count(b"\n")
        while n_newlines < target_lines and not self._eof:
            blob = self._fh.read(self._READ_SIZE)
            if not blob:
                self._eof = True
                break
            parts.append(blob)
            n_newlines += blob.count(b"\n")
        buf = b"".join(parts)
        if not buf:
            self._fh.close()
            raise StopIteration

        arr = np.frombuffer(buf, dtype=np.uint8)
        nl = np.flatnonzero(arr == 10)
        # Trailing bytes without a newline: at EOF treat them as a final
        # line; otherwise carry them to the next chunk.
        n_lines_avail = len(nl)
        trailing = len(buf) - (int(nl[-1]) + 1 if n_lines_avail else 0)
        if self._eof and trailing > 0:
            n_lines_avail += 1

        nrec = min(n_lines_avail // 4, self._chunk)
        if nrec == 0:
            if self._eof:
                if n_lines_avail > 0:
                    raise ValueError("Truncated FASTQ record at end of file")
                self._fh.close()
                raise StopIteration
            raise ValueError("FASTQ chunk contains no complete record")

        nl_used = nl[: nrec * 4]
        if len(nl_used) < nrec * 4:
            # Final line lacks a newline (EOF without trailing \n)
            nl_used = np.append(nl_used, len(buf))
        cut = int(nl_used[-1]) + 1
        self._leftover = buf[cut:] if cut < len(buf) else b""
        if not self._eof and self._leftover == b"":
            self._leftover = b""

        starts = np.empty(nrec * 4, dtype=np.int64)
        starts[0] = 0
        starts[1:] = nl_used[:-1] + 1
        ends = nl_used.astype(np.int64).copy()
        # Strip \r from CRLF line endings
        cr = (ends > starts) & (arr[np.minimum(ends - 1, len(buf) - 1)] == 13) \
            if len(buf) else np.zeros(0, bool)
        ends[cr] -= 1
        return buf, starts, ends


def _decode_windows(buf, starts, ends, max_w):
    """Gather variable-length byte windows into a padded 2D uint8 matrix.

    Returns (matrix, valid_mask): positions beyond each window's end hold 0
    in ``matrix`` and False in ``valid_mask``.
    """
    arr = np.frombuffer(buf, dtype=np.uint8)
    n = len(starts)
    widths = ends - starts
    col = np.arange(max_w, dtype=np.int64)
    idx = starts[:, None] + col[None, :]
    valid = col[None, :] < widths[:, None]
    np.minimum(idx, len(arr) - 1, out=idx)
    mat = arr[idx]
    mat[~valid] = 0
    return mat, valid


def _first_low_qual(qual_start, qual_end, low_pos):
    """For each read window, the absolute position of its first low-quality
    base, or -1.  ``low_pos`` is a sorted array of absolute positions."""
    n = len(qual_start)
    out = np.full(n, -1, dtype=np.int64)
    if len(low_pos) == 0 or n == 0:
        return out
    first_idx = np.searchsorted(low_pos, qual_start, side="left")
    has = first_idx < len(low_pos)
    cand = np.where(has, low_pos[np.minimum(first_idx, len(low_pos) - 1)], -1)
    hit = has & (cand < qual_end)
    out[hit] = cand[hit]
    return out


class _ReadBatch:
    """One chunk's worth of reads from a single FASTQ file, with mutable
    per-read trim windows.  All arrays are per-read (length n)."""

    def __init__(self, buf, line_starts, line_ends):
        self.buf = buf
        nrec = len(line_starts) // 4
        self.n = nrec
        ls = line_starts.reshape(nrec, 4)
        le = line_ends.reshape(nrec, 4)
        self.head_s, self.head_e = ls[:, 0].copy(), le[:, 0].copy()
        self.seq_s, self.seq_e = ls[:, 1].copy(), le[:, 1].copy()
        self.qual_s, self.qual_e = ls[:, 3].copy(), le[:, 3].copy()
        if not np.array_equal(self.seq_e - self.seq_s,
                              self.qual_e - self.qual_s):
            raise ValueError("FASTQ record with mismatched sequence and "
                             "quality lengths")
        # RC'd reads (orient_fwd) are stored out-of-band as full strings
        # plus the window start at creation time; the live trim window
        # (seq_s/seq_e) is interpreted relative to that origin.
        self.override = {}  # read index -> (seq bytes, qual bytes, orig_s)

    def widths(self):
        return self.seq_e - self.seq_s

    def trim_left_(self, k):
        self.seq_s = self.seq_s + k
        self.qual_s = self.qual_s + k

    def trim_right_(self, k):
        self.seq_e = self.seq_e - k
        self.qual_e = self.qual_e - k

    def truncate_(self, length, mask=None):
        """Truncate reads to ``length`` (post-trim) where mask is True."""
        new_e = self.seq_s + length
        if mask is None:
            take = new_e < self.seq_e
        else:
            take = mask & (new_e < self.seq_e)
        self.seq_e = np.where(take, new_e, self.seq_e)
        self.qual_e = np.where(take, self.qual_s + length, self.qual_e)

    def apply_override_rc(self, indices):
        """Reverse-complement the given reads (stored as byte overrides)."""
        for i in indices:
            i = int(i)
            s = self.buf[self.seq_s[i]:self.seq_e[i]]
            q = self.buf[self.qual_s[i]:self.qual_e[i]]
            self.override[i] = (_rc_bytes(s), q[::-1], int(self.seq_s[i]))

    def get_seq_qual(self, i):
        ov = self.override.get(i)
        if ov is not None:
            s, q, orig = ov
            rel_s = int(self.seq_s[i]) - orig
            rel_e = int(self.seq_e[i]) - orig
            return s[rel_s:rel_e], q[rel_s:rel_e]
        return (self.buf[self.seq_s[i]:self.seq_e[i]],
                self.buf[self.qual_s[i]:self.qual_e[i]])

    def get_header(self, i):
        return self.buf[self.head_s[i]:self.head_e[i]]

    def seq_matrix(self, indices, max_w):
        """Padded seq/qual matrices for the given reads.  Overridden reads
        (RC'd by orient_fwd) are patched in afterwards."""
        smat, valid = _decode_windows(self.buf, self.seq_s[indices],
                                      self.seq_e[indices], max_w)
        qmat, _ = _decode_windows(self.buf, self.qual_s[indices],
                                  self.qual_e[indices], max_w)
        if self.override:
            for row, i in enumerate(indices):
                if int(i) in self.override:
                    s, q = self.get_seq_qual(int(i))
                    w = len(s)
                    smat[row, :w] = np.frombuffer(s, dtype=np.uint8)
                    smat[row, w:] = 0
                    qmat[row, :w] = np.frombuffer(q, dtype=np.uint8)
                    qmat[row, w:] = 0
                    valid[row, :w] = True
                    valid[row, w:] = False
        return smat, qmat, valid

    def truncq_cut_(self, tq_byte, active):
        """R's trimTails(k=1): truncate each active read at its first base
        with quality byte <= tq_byte (within the current trim window)."""
        arr = np.frombuffer(self.buf, dtype=np.uint8)
        low = np.flatnonzero(arr <= tq_byte)
        first = _first_low_qual(self.qual_s, self.qual_e, low)
        cut = active & (first >= 0)
        new_qe = np.where(cut, first, self.qual_e)
        self.seq_e = self.seq_s + (new_qe - self.qual_s)
        self.qual_e = new_qe
        # Overridden (RC'd) reads scan their own byte strings within the
        # live window
        for i in list(self.override):
            if not active[i]:
                continue
            _, q = self.get_seq_qual(i)
            for p, b in enumerate(q):
                if b <= tq_byte:
                    self.seq_e[i] = self.seq_s[i] + p
                    self.qual_e[i] = self.qual_s[i] + p
                    break


def _as_pair(val):
    """Normalise a scalar or length-2 sequence to a 2-tuple."""
    if isinstance(val, (list, tuple)):
        if len(val) == 2:
            return tuple(val)
        if len(val) == 1:
            return (val[0], val[0])
        raise ValueError(f"Expected length-2 tuple, got {len(val)}")
    return (val, val)


def _resolve_trunc_ends(trunc_len, trim_left):
    """R's endF computation: truncation length measured from the original
    5' end.  Returns the post-trim target length, or None when disabled."""
    start = max(1, trim_left + 1)
    end = trunc_len
    if end < start:
        return None
    return end - start + 1


def _phix_flags(batch: _ReadBatch, indices) -> np.ndarray:
    """isPhiX flags for the given (trimmed) reads of a batch."""
    if batch.override and any(int(i) in batch.override for i in indices):
        seqs = [batch.get_seq_qual(int(i))[0] for i in indices]
        return _is_phix(seqs)
    import ctypes as ct
    from . import _cdada
    from .utils import _phix_genome, _rc
    genome = _phix_genome()
    n = len(indices)
    starts = np.ascontiguousarray(batch.seq_s[indices], dtype=np.int64)
    ends = np.ascontiguousarray(batch.seq_e[indices], dtype=np.int64)
    flags = np.zeros(n, dtype=bool)
    for ref in (genome, _rc(genome)):
        counts = np.zeros(n, dtype=np.int32)
        _cdada._lib.dada2_match_ref_windows(
            batch.buf,
            starts.ctypes.data_as(ct.POINTER(ct.c_int64)),
            ends.ctypes.data_as(ct.POINTER(ct.c_int64)),
            ct.c_int(n), ref.encode("ascii"), ct.c_int(16), ct.c_int(1),
            counts.ctypes.data_as(ct.POINTER(ct.c_int32)))
        flags |= counts >= 2
    return flags


def _complexity_low(batch: _ReadBatch, indices, threshold) -> np.ndarray:
    seqs = [batch.get_seq_qual(int(i))[0].decode("ascii") for i in indices]
    return _seq_complexity(seqs) < threshold


def _ee_and_minq(batch: _ReadBatch, indices):
    """Per-read expected errors (sequential left-to-right sum, matching
    C_matrixEE) and minimum quality score, for the given reads."""
    if len(indices) == 0:
        return np.zeros(0), np.zeros(0)
    widths = (batch.seq_e - batch.seq_s)[indices]
    max_w = int(widths.max()) if len(widths) else 0
    if max_w == 0:
        return np.zeros(len(indices)), np.full(len(indices), np.inf)
    _, qmat, valid = batch.seq_matrix(indices, max_w)
    ee = np.cumsum(_EE_LUT[qmat], axis=1)[:, -1]
    q_int = np.where(valid, qmat.astype(np.int16) - 33, np.int16(32767))
    minq = q_int.min(axis=1).astype(np.float64)
    minq[widths == 0] = np.inf
    return ee, minq


def _maxn_counts(batch: _ReadBatch, indices):
    if len(indices) == 0:
        return np.zeros(0, dtype=np.int64)
    widths = (batch.seq_e - batch.seq_s)[indices]
    max_w = int(widths.max()) if len(widths) else 0
    if max_w == 0:
        return np.zeros(len(indices), dtype=np.int64)
    smat, _, valid = batch.seq_matrix(indices, max_w)
    other = _IS_OTHER[smat] & valid
    return other.sum(axis=1)


# ---------------------------------------------------------------------------
# fastq_filter  --  single-end filtering
# ---------------------------------------------------------------------------

def fastq_filter(
    fn: str,
    fout: str,
    *,
    trim_left: int = 0,
    trim_right: int = 0,
    trunc_len: int = 0,
    trunc_q: Optional[int] = 2,
    max_len: float = 0,
    min_len: int = 20,
    max_n: int = 0,
    min_q: int = 0,
    max_ee: float = float("inf"),
    rm_phix: bool = True,
    rm_lowcomplex: float = 0.0,
    orient_fwd: Optional[str] = None,
    compress: bool = True,
    n: int = 1_000_000,
    verbose: bool = False,
) -> Tuple[int, int]:
    """Filter and trim a single FASTQ file (mirrors R's fastqFilter).

    Order of operations (identical to R): orient_fwd, max_len, trim_left,
    trim_right, trunc_q, trunc_len, min_len, max_n, min_q, max_ee,
    rm_phix, rm_lowcomplex.  Note that ``trunc_len`` counts bases from the
    *original* 5' end of the read: with ``trim_left=10, trunc_len=240``
    the output reads are 230 bases long.

    Args:
        fn: Input FASTQ path (plain or gzipped).
        fout: Output FASTQ path.
        trim_left: Bases to trim from the 5' end.  Reads shorter than
            ``trim_left + 1`` are discarded.
        trim_right: Bases to trim from the 3' end.  Reads with no bases
            left are discarded.
        trunc_len: Truncate reads at this position (counted from the
            original 5' end).  Shorter reads are discarded.  0 disables.
        trunc_q: Truncate at the first base with quality <= this value.
            ``None`` disables (R has no equivalent; its default is 2).
        max_len: Discard reads longer than this before trimming.  0 or
            ``inf`` disables.
        min_len: Discard reads shorter than this after all trimming.
        max_n: Maximum number of non-ACGT bases allowed.
        min_q: After truncation, discard reads whose minimum quality is
            not strictly greater than this.  Only applied when
            ``min_q > trunc_q`` (R's behaviour).
        max_ee: Maximum expected errors (``sum(10^(-Q/10))``).
        rm_phix: Remove reads matching the PhiX genome.
        rm_lowcomplex: Remove reads with sequence complexity below this.
            0 disables.
        orient_fwd: If set, keep only reads that begin with this sequence
            in forward or reverse-complement orientation, re-orienting the
            latter.  Within each chunk, forward-matching reads are output
            before re-oriented ones (R's behaviour).
        compress: Gzip-compress the output.
        n: Number of records processed per chunk (R's ``n``).  Output
            order within a chunk can depend on this when ``orient_fwd``
            is set.
        verbose: Print a summary line when done.

    Returns:
        ``(reads_in, reads_out)``
    """
    if trunc_q is not None and not (0 <= trunc_q <= 93):
        raise ValueError("trunc_q must be within the Phred+33 range 0..93")

    reads_in = 0
    reads_out = 0
    end_len = _resolve_trunc_ends(int(trunc_len), int(trim_left))
    orient_b = orient_fwd.encode("ascii") if orient_fwd is not None else None

    writer = _GzipRecordWriter(fout, compress)
    try:
        for buf, ls, le in _FastqChunkReader(fn, n):
            batch = _ReadBatch(buf, ls, le)
            reads_in += batch.n
            order = np.arange(batch.n)
            keep = np.ones(batch.n, dtype=bool)

            # --- orient_fwd (fwd matches first, then RC matches) ---
            if orient_b is not None:
                blen = len(orient_b)
                if (batch.widths() < blen).any():
                    # R's narrow() fails on reads shorter than orient.fwd
                    raise ValueError(
                        "orient_fwd is longer than a read in the input "
                        "(R's fastqFilter fails on such input too)")
                fwd_hit = np.zeros(batch.n, dtype=bool)
                rc_hit = np.zeros(batch.n, dtype=bool)
                for i in range(batch.n):
                    s = batch.buf[batch.seq_s[i]:batch.seq_e[i]]
                    if s[:blen] == orient_b:
                        fwd_hit[i] = True
                    elif _rc_bytes(s)[:blen] == orient_b:
                        rc_hit[i] = True
                batch.apply_override_rc(np.flatnonzero(rc_hit))
                keep &= fwd_hit | rc_hit
                order = np.concatenate([np.flatnonzero(fwd_hit),
                                        np.flatnonzero(rc_hit)])
                mask = np.zeros(batch.n, dtype=bool)
                mask[order] = True
                keep &= mask

            widths = batch.widths()

            # --- max_len (before trimming) ---
            if max_len and np.isfinite(max_len) and max_len > 0:
                keep &= widths <= max_len

            # --- trim_left (reads shorter than trim_left+1 are dropped) ---
            if trim_left > 0:
                keep &= widths >= trim_left + 1
                batch.trim_left_(trim_left)
            # --- trim_right ---
            if trim_right > 0:
                keep &= batch.widths() > trim_right
                batch.trim_right_(trim_right)

            # --- trunc_q ---
            if trunc_q is not None:
                batch.truncq_cut_(trunc_q + 33, keep)

            # --- trunc_len ---
            if end_len is not None:
                keep &= batch.widths() >= end_len
                batch.truncate_(end_len)

            # --- min_len ---
            keep &= batch.widths() >= min_len

            # --- max_n ---
            idx = order[keep[order]]
            if len(idx):
                nn = _maxn_counts(batch, idx)
                bad = idx[nn > max_n]
                keep[bad] = False

            # --- min_q / max_ee ---
            idx = order[keep[order]]
            if len(idx) and (max_ee < float("inf") or
                             (min_q > (trunc_q if trunc_q is not None else -1))):
                ee, minq = _ee_and_minq(batch, idx)
                bad = np.zeros(len(idx), dtype=bool)
                if min_q > (trunc_q if trunc_q is not None else -1) and min_q > 0:
                    bad |= ~(minq > min_q)
                if max_ee < float("inf"):
                    bad |= ~(ee <= max_ee)
                keep[idx[bad]] = False

            # --- rm_phix ---
            idx = order[keep[order]]
            if rm_phix and len(idx):
                keep[idx[_phix_flags(batch, idx)]] = False

            # --- rm_lowcomplex ---
            idx = order[keep[order]]
            if rm_lowcomplex > 0 and len(idx):
                keep[idx[_complexity_low(batch, idx, rm_lowcomplex)]] = False

            # --- write passing reads in chunk order ---
            out_idx = order[keep[order]]
            if len(out_idx):
                parts = []
                for i in out_idx:
                    s, q = batch.get_seq_qual(int(i))
                    parts += [batch.get_header(int(i)), b"\n", s, b"\n+\n",
                              q, b"\n"]
                writer.write(b"".join(parts))
                reads_out += len(out_idx)
    finally:
        writer.close()

    if reads_out == 0:
        print(f"The filter removed all reads: {writer.path} not written.")
        writer.remove()

    if verbose:
        pct = (reads_out / reads_in * 100) if reads_in > 0 else 0.0
        print(f"Read in {reads_in} reads, output {reads_out} "
              f"({pct:.1f}%) filtered reads.")

    return reads_in, reads_out


# ---------------------------------------------------------------------------
# fastq_paired_filter  --  paired-end filtering
# ---------------------------------------------------------------------------

def fastq_paired_filter(
    fwd: str,
    filt_fwd: str,
    rev: str,
    filt_rev: str,
    *,
    trim_left: Union[int, Tuple[int, int]] = (0, 0),
    trim_right: Union[int, Tuple[int, int]] = (0, 0),
    trunc_len: Union[int, Tuple[int, int]] = (0, 0),
    trunc_q: Union[None, int, Tuple[Optional[int], Optional[int]]] = 2,
    max_len: Union[float, Tuple[float, float]] = (0, 0),
    min_len: Union[int, Tuple[int, int]] = (20, 20),
    max_n: Union[int, Tuple[int, int]] = (0, 0),
    min_q: Union[int, Tuple[int, int]] = (0, 0),
    max_ee: Union[float, Tuple[float, float]] = (float("inf"), float("inf")),
    rm_phix: Union[bool, Tuple[bool, bool]] = True,
    rm_lowcomplex: Union[float, Tuple[float, float]] = (0.0, 0.0),
    orient_fwd: Optional[str] = None,
    compress: bool = True,
    n: int = 1_000_000,
    verbose: bool = False,
) -> Tuple[int, int]:
    """Filter and trim paired FASTQ files (mirrors R's fastqPairedFilter).

    Both reads of a pair must pass all filters for the pair to be kept.
    Parameters accept scalars (applied to both) or ``(fwd, rev)`` tuples.
    See :func:`fastq_filter` for per-parameter documentation; the order of
    operations and boundary behaviour match R exactly, including:

    - ``trunc_len`` counting from the original 5' end of each read.
    - ``orient_fwd``: pairs whose *reverse* read starts with the primer
      have their forward/reverse reads swapped (not reverse-complemented);
      pairs matching in neither read are discarded.  Within each chunk,
      forward-matching pairs are output before swapped pairs.
    - ``rm_phix``: a pair is discarded when either read matches PhiX
      (when enabled for both reads).

    Returns:
        ``(reads_in, reads_out)`` counted in pairs.
    """
    tl = _as_pair(trim_left)
    tr = _as_pair(trim_right)
    tlen = _as_pair(trunc_len)
    tq = _as_pair(trunc_q)
    mxl = _as_pair(max_len)
    mnl = _as_pair(min_len)
    mn = _as_pair(max_n)
    mq = _as_pair(min_q)
    mee = _as_pair(max_ee)
    rphix = _as_pair(rm_phix)
    rlc = _as_pair(rm_lowcomplex)
    for q in tq:
        if q is not None and not (0 <= q <= 93):
            raise ValueError("trunc_q must be within the Phred+33 range 0..93")

    end_len = (
        _resolve_trunc_ends(int(tlen[0]), int(tl[0])),
        _resolve_trunc_ends(int(tlen[1]), int(tl[1])),
    )
    orient_b = orient_fwd.encode("ascii") if orient_fwd is not None else None

    reads_in = 0
    reads_out = 0

    writer_f = _GzipRecordWriter(filt_fwd, compress)
    writer_r = _GzipRecordWriter(filt_rev, compress)
    try:
        reader_f = _FastqChunkReader(fwd, n)
        reader_r = _FastqChunkReader(rev, n)
        while True:
            chunk_f = next(reader_f, None)
            chunk_r = next(reader_r, None)
            if chunk_f is None and chunk_r is None:
                break
            if chunk_f is None or chunk_r is None:
                raise ValueError(
                    f"Forward and reverse FASTQ files have different numbers "
                    f"of reads ({fwd}, {rev}). Files must be synchronized.")
            bf = _ReadBatch(*chunk_f)
            br = _ReadBatch(*chunk_r)
            if bf.n != br.n:
                raise ValueError(
                    f"Mismatched forward and reverse sequence files: "
                    f"{bf.n}, {br.n}.")
            npairs = bf.n
            reads_in += npairs
            keep = np.ones(npairs, dtype=bool)
            order = np.arange(npairs)
            swapped = np.zeros(npairs, dtype=bool)

            # --- orient_fwd: F matches keep orientation; R matches swap ---
            if orient_b is not None:
                blen = len(orient_b)
                if (bf.widths() < blen).any() or (br.widths() < blen).any():
                    raise ValueError(
                        "orient_fwd is longer than a read in the input "
                        "(R's fastqPairedFilter fails on such input too)")
                keep_f = np.zeros(npairs, dtype=bool)
                keep_r = np.zeros(npairs, dtype=bool)
                for i in range(npairs):
                    sf = bf.buf[bf.seq_s[i]:bf.seq_s[i] + blen]
                    if sf == orient_b:
                        keep_f[i] = True
                    else:
                        sr = br.buf[br.seq_s[i]:br.seq_s[i] + blen]
                        if sr == orient_b:
                            keep_r[i] = True
                keep &= keep_f | keep_r
                swapped = keep_r
                order = np.concatenate([np.flatnonzero(keep_f),
                                        np.flatnonzero(keep_r)])

            def side(idx_batch):
                """Batch serving the forward (0) or reverse (1) role,
                honouring per-pair swaps."""
                return (bf, br) if idx_batch == 0 else (br, bf)

            def role_widths(role):
                base, alt = (bf, br) if role == 0 else (br, bf)
                w = base.widths().copy()
                if swapped.any():
                    w[swapped] = alt.widths()[swapped]
                return w

            # For swapped pairs the F-role read lives in the reverse batch.
            # Implement all window mutations on both batches, selecting
            # per-pair by role.
            class _Role:
                def __init__(self, role):
                    self.role = role  # 0 = forward role, 1 = reverse role

                def _sel(self):
                    # (batch, mask) pairs contributing to this role
                    if self.role == 0:
                        return ((bf, ~swapped), (br, swapped))
                    return ((br, ~swapped), (bf, swapped))

                def widths(self):
                    w = np.empty(npairs, dtype=np.int64)
                    for b, m in self._sel():
                        w[m] = b.widths()[m]
                    return w

                def trim_left_(self, k):
                    for b, m in self._sel():
                        b.seq_s[m] += k
                        b.qual_s[m] += k

                def trim_right_(self, k):
                    for b, m in self._sel():
                        b.seq_e[m] -= k
                        b.qual_e[m] -= k

                def truncq_cut_(self, tq_byte, active):
                    for b, m in self._sel():
                        b.truncq_cut_(tq_byte, active & m)

                def truncate_(self, length):
                    for b, m in self._sel():
                        b.truncate_(length, mask=m)

                def gather(self, indices):
                    """(batch, local_indices) pairs covering ``indices``."""
                    out = []
                    for b, m in self._sel():
                        sel = m[indices]
                        if sel.any():
                            out.append((b, indices[sel], sel))
                    return out

            roleF, roleR = _Role(0), _Role(1)

            # --- max_len ---
            for role, ml in ((roleF, mxl[0]), (roleR, mxl[1])):
                if ml and np.isfinite(ml) and ml > 0:
                    keep &= role.widths() <= ml

            # --- trim_left (pair dropped when either read too short) ---
            keep &= (roleF.widths() >= tl[0] + 1) & (roleR.widths() >= tl[1] + 1)
            if tl[0] > 0:
                roleF.trim_left_(tl[0])
            if tl[1] > 0:
                roleR.trim_left_(tl[1])

            # --- trim_right ---
            if tr[0] > 0:
                keep &= roleF.widths() > tr[0]
                roleF.trim_right_(tr[0])
            if tr[1] > 0:
                keep &= roleR.widths() > tr[1]
                roleR.trim_right_(tr[1])

            # --- trunc_q, then drop zero-width pairs ---
            if tq[0] is not None:
                roleF.truncq_cut_(tq[0] + 33, keep)
            if tq[1] is not None:
                roleR.truncq_cut_(tq[1] + 33, keep)
            keep &= (roleF.widths() > 0) & (roleR.widths() > 0)

            # --- trunc_len ---
            if end_len[0] is not None:
                keep &= roleF.widths() >= end_len[0]
            if end_len[1] is not None:
                keep &= roleR.widths() >= end_len[1]
            if end_len[0] is not None:
                roleF.truncate_(end_len[0])
            if end_len[1] is not None:
                roleR.truncate_(end_len[1])

            # --- min_len ---
            keep &= (roleF.widths() >= mnl[0]) & (roleR.widths() >= mnl[1])

            # --- max_n ---
            idx = order[keep[order]]
            if len(idx):
                bad = np.zeros(len(idx), dtype=bool)
                for role, mx in ((roleF, mn[0]), (roleR, mn[1])):
                    for b, gidx, sel in role.gather(idx):
                        counts = _maxn_counts(b, gidx)
                        sub = np.zeros(len(idx), dtype=bool)
                        sub[np.flatnonzero(sel)] = counts > mx
                        bad |= sub
                keep[idx[bad]] = False

            # --- min_q / max_ee ---
            idx = order[keep[order]]
            if len(idx):
                bad = np.zeros(len(idx), dtype=bool)
                for role, (q_min, q_trunc, ee_max) in (
                        (roleF, (mq[0], tq[0], mee[0])),
                        (roleR, (mq[1], tq[1], mee[1]))):
                    tq_eff = q_trunc if q_trunc is not None else -1
                    need_minq = q_min > tq_eff and q_min > 0
                    need_ee = ee_max < float("inf")
                    if not (need_minq or need_ee):
                        continue
                    for b, gidx, sel in role.gather(idx):
                        ee, minq = _ee_and_minq(b, gidx)
                        sub = np.zeros(len(idx), dtype=bool)
                        flags = np.zeros(len(gidx), dtype=bool)
                        if need_minq:
                            flags |= ~(minq > q_min)
                        if need_ee:
                            flags |= ~(ee <= ee_max)
                        sub[np.flatnonzero(sel)] = flags
                        bad |= sub
                keep[idx[bad]] = False

            # --- rm_phix (pair removed when either enabled read matches) ---
            idx = order[keep[order]]
            if len(idx) and (rphix[0] or rphix[1]):
                is_phi = np.zeros(len(idx), dtype=bool)
                for role, enabled in ((roleF, rphix[0]), (roleR, rphix[1])):
                    if not enabled:
                        continue
                    for b, gidx, sel in role.gather(idx):
                        flags = _phix_flags(b, gidx)
                        sub = np.zeros(len(idx), dtype=bool)
                        sub[np.flatnonzero(sel)] = flags
                        is_phi |= sub
                keep[idx[is_phi]] = False

            # --- rm_lowcomplex ---
            idx = order[keep[order]]
            if len(idx) and (rlc[0] > 0 or rlc[1] > 0):
                is_lowc = np.zeros(len(idx), dtype=bool)
                for role, thr in ((roleF, rlc[0]), (roleR, rlc[1])):
                    if thr <= 0:
                        continue
                    for b, gidx, sel in role.gather(idx):
                        flags = _complexity_low(b, gidx, thr)
                        sub = np.zeros(len(idx), dtype=bool)
                        sub[np.flatnonzero(sel)] = flags
                        is_lowc |= sub
                keep[idx[is_lowc]] = False

            # --- write passing pairs in chunk order ---
            out_idx = order[keep[order]]
            if len(out_idx):
                parts_f, parts_r = [], []
                for i in out_idx:
                    i = int(i)
                    b_f, b_r = (br, bf) if swapped[i] else (bf, br)
                    sf, qf = b_f.get_seq_qual(i)
                    sr, qr = b_r.get_seq_qual(i)
                    parts_f += [b_f.get_header(i), b"\n", sf, b"\n+\n", qf, b"\n"]
                    parts_r += [b_r.get_header(i), b"\n", sr, b"\n+\n", qr, b"\n"]
                writer_f.write(b"".join(parts_f))
                writer_r.write(b"".join(parts_r))
                reads_out += len(out_idx)
    finally:
        writer_f.close()
        writer_r.close()

    if reads_out == 0:
        print(f"The filter removed all reads: {writer_f.path} and "
              f"{writer_r.path} not written.")
        writer_f.remove()
        writer_r.remove()

    if verbose:
        pct = (reads_out / reads_in * 100) if reads_in > 0 else 0.0
        print(f"Read in {reads_in} paired-reads, output {reads_out} "
              f"({pct:.1f}%) filtered paired-reads.")

    return reads_in, reads_out


# ---------------------------------------------------------------------------
# Module-level task runners (picklable for ProcessPoolExecutor)
# ---------------------------------------------------------------------------

def _run_paired_task(args):
    """Run fastq_paired_filter on a single file pair (picklable)."""
    fwd, filt_fwd, rev, filt_rev, kw = args
    return fastq_paired_filter(fwd, filt_fwd, rev, filt_rev, **kw)


def _run_single_task(args):
    """Run fastq_filter on a single file (picklable)."""
    fwd, filt, kw = args
    return fastq_filter(fwd, filt, **kw)


# ---------------------------------------------------------------------------
# filter_and_trim  --  convenience wrapper
# ---------------------------------------------------------------------------

def filter_and_trim(
    fwd: Union[str, Sequence[str]],
    filt: Union[str, Sequence[str]],
    rev: Optional[Union[str, Sequence[str]]] = None,
    filt_rev: Optional[Union[str, Sequence[str]]] = None,
    *,
    trim_left: Union[int, Tuple[int, int]] = 0,
    trim_right: Union[int, Tuple[int, int]] = 0,
    trunc_len: Union[int, Tuple[int, int]] = 0,
    trunc_q: Union[None, int, Tuple[Optional[int], Optional[int]]] = 2,
    max_len: Union[float, Tuple[float, float]] = 0,
    min_len: Union[int, Tuple[int, int]] = 20,
    max_n: Union[int, Tuple[int, int]] = 0,
    min_q: Union[int, Tuple[int, int]] = 0,
    max_ee: Union[float, Tuple[float, float]] = float("inf"),
    rm_phix: Union[bool, Tuple[bool, bool]] = True,
    rm_lowcomplex: Union[float, Tuple[float, float]] = 0.0,
    orient_fwd: Optional[str] = None,
    compress: bool = True,
    n: int = 1_000_000,
    multithread: Union[bool, int] = False,
    verbose: bool = False,
) -> np.ndarray:
    """Filter and trim FASTQ files (single- or paired-end).

    This is a convenience wrapper around :func:`fastq_filter` (single-end)
    and :func:`fastq_paired_filter` (paired-end).  It dispatches to the
    appropriate function based on whether *rev* is provided and optionally
    parallelises across files with
    :class:`~concurrent.futures.ProcessPoolExecutor`.

    Args:
        fwd: Path(s) to forward (or single-end) input FASTQ file(s).
        filt: Path(s) to filtered output FASTQ file(s), same length as
            ``fwd``.
        rev: Path(s) to reverse-read input FASTQ file(s), or ``None`` for
            single-end mode.
        filt_rev: Path(s) to filtered reverse-read output FASTQ file(s).
            Required when ``rev`` is given.
        multithread: ``True`` uses all cores; an int > 1 uses that many
            worker processes; ``False``/1 disables parallelism.
        verbose: Print per-file summaries.

    All other parameters are forwarded — see :func:`fastq_filter` and
    :func:`fastq_paired_filter`.  Defaults match R's ``filterAndTrim``
    (in particular ``trunc_q=2``).

    Returns:
        Integer array of shape ``(n_files, 2)`` with columns
        ``[reads_in, reads_out]``.
    """
    # Normalise to lists
    if isinstance(fwd, str):
        fwd = [fwd]
    if isinstance(filt, str):
        filt = [filt]
    if isinstance(rev, str):
        rev = [rev]
    if isinstance(filt_rev, str):
        filt_rev = [filt_rev]

    n_files = len(fwd)
    if len(filt) != n_files:
        raise ValueError(
            f"fwd ({n_files} files) and filt ({len(filt)} files) must have "
            f"the same length.")

    paired = rev is not None
    if paired:
        if filt_rev is None:
            raise ValueError("filt_rev is required when rev is provided.")
        if len(rev) != n_files:
            raise ValueError(
                f"fwd ({n_files}) and rev ({len(rev)}) must have the same length.")
        if len(filt_rev) != n_files:
            raise ValueError(
                f"fwd ({n_files}) and filt_rev ({len(filt_rev)}) must have "
                f"the same length.")

    # Determine number of workers
    if multithread is True:
        n_workers = os.cpu_count() or 1
    elif isinstance(multithread, int) and multithread > 1:
        n_workers = multithread
    else:
        n_workers = 1

    # Common kwargs (shared by all files)
    common_kw = dict(
        trim_left=trim_left,
        trim_right=trim_right,
        trunc_len=trunc_len,
        trunc_q=trunc_q,
        max_len=max_len,
        min_len=min_len,
        max_n=max_n,
        min_q=min_q,
        max_ee=max_ee,
        rm_phix=rm_phix,
        rm_lowcomplex=rm_lowcomplex,
        orient_fwd=orient_fwd,
        compress=compress,
        n=n,
        verbose=verbose,
    )

    results = np.zeros((n_files, 2), dtype=np.int64)

    # Build per-file argument tuples (picklable for ProcessPoolExecutor)
    if paired:
        tasks = [
            (fwd[i], filt[i], rev[i], filt_rev[i], common_kw)
            for i in range(n_files)
        ]
    else:
        # Normalise tuple params to scalars for single-end
        se_kw = dict(common_kw)
        for key in (
            "trim_left", "trim_right", "trunc_len", "trunc_q",
            "max_len", "min_len", "max_n", "min_q", "max_ee",
            "rm_phix", "rm_lowcomplex",
        ):
            v = se_kw[key]
            if isinstance(v, (list, tuple)):
                se_kw[key] = v[0]
        tasks = [(fwd[i], filt[i], se_kw) for i in range(n_files)]

    if n_workers > 1 and n_files > 1:
        if _loky_executor is not None:
            pool_cm = _loky_executor(max_workers=n_workers)
        else:
            # spawn: forking after an OpenMP region deadlocks the child
            ctx = multiprocessing.get_context("spawn")
            pool_cm = ProcessPoolExecutor(max_workers=n_workers, mp_context=ctx)
        with pool_cm as pool:
            runner = _run_paired_task if paired else _run_single_task
            futures = {pool.submit(runner, t): i for i, t in enumerate(tasks)}
            for fut in futures:
                results[futures[fut]] = fut.result()
    else:
        for i, t in enumerate(tasks):
            results[i] = _run_paired_task(t) if paired else _run_single_task(t)

    return results
