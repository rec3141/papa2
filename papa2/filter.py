"""FASTQ quality filtering and trimming.

Python port of dada2's ``filterAndTrim``, ``fastqFilter``, and
``fastqPairedFilter`` functions.  Reads are streamed four lines at a time
so memory usage is independent of file size.
"""

from __future__ import annotations

import gzip
import math
import os
from concurrent.futures import ProcessPoolExecutor
from typing import Optional, Sequence, Tuple, Union

import numpy as np

from ._cdada import rc as _rc
from .utils import is_phix as _is_phix
from .utils import seq_complexity as _seq_complexity

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

_ACGT = frozenset("ACGTacgt")


def _open_fq(path: str, mode: str = "rt"):
    """Open a FASTQ file, transparently handling gzip."""
    if path.endswith(".gz") or path.endswith(".gzip"):
        return gzip.open(path, mode, encoding="ascii" if "t" in mode else None)
    return open(path, mode)


def _write_opener(path: str, compress: bool):
    """Return a writable file handle, gzipped if requested."""
    if compress:
        if not path.endswith(".gz"):
            path += ".gz"
        return gzip.open(path, "wt", encoding="ascii", compresslevel=6)
    return open(path, "w")


def _fastq_records(fh):
    """Yield (header, seq, sep, qual) tuples from an open FASTQ handle.

    Strips trailing whitespace from each line.  Yields nothing for empty
    files.
    """
    while True:
        header = fh.readline()
        if not header:
            return
        header = header.rstrip("\n\r")
        seq = fh.readline().rstrip("\n\r")
        sep = fh.readline().rstrip("\n\r")
        qual = fh.readline().rstrip("\n\r")
        if not header:
            return
        yield header, seq, sep, qual


def _count_non_acgt(seq: str) -> int:
    """Count bases that are not A, C, G, or T (case-insensitive)."""
    return sum(1 for ch in seq if ch not in _ACGT)


def _decode_quals(qual_str: str) -> list[int]:
    """Decode a Phred+33 quality string to a list of integer Q scores."""
    return [ord(ch) - 33 for ch in qual_str]


def _expected_errors(quals: list[int]) -> float:
    """Compute expected errors from a list of Phred quality scores."""
    return sum(10.0 ** (-q / 10.0) for q in quals)


def _orient_read(seq: str, qual_str: str, orient_fwd: str) -> Optional[Tuple[str, str]]:
    """Orient a read to start with *orient_fwd*.

    Returns ``(seq, qual_str)`` possibly reverse-complemented, or ``None``
    if neither orientation matches.
    """
    orient_upper = orient_fwd.upper()
    n = len(orient_upper)
    if seq[:n].upper() == orient_upper:
        return seq, qual_str
    # Check RC
    seq_rc = _rc(seq)
    if seq_rc[:n].upper() == orient_upper:
        return seq_rc, qual_str[::-1]
    return None


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
    trunc_q: Optional[int] = None,
    max_len: int = 0,
    min_len: int = 20,
    max_n: int = 0,
    min_q: int = 0,
    max_ee: float = float("inf"),
    rm_phix: bool = True,
    rm_lowcomplex: float = 0.0,
    orient_fwd: Optional[str] = None,
    compress: bool = True,
    verbose: bool = False,
) -> Tuple[int, int]:
    """Filter and trim a single FASTQ file.

    Parameters
    ----------
    fn : str
        Path to input FASTQ file (plain or gzipped).
    fout : str
        Path to output FASTQ file.
    trim_left : int
        Number of bases to trim from the 5' end.
    trim_right : int
        Number of bases to trim from the 3' end.
    trunc_len : int
        Truncate reads to exactly this length after trimming.  Reads
        shorter than *trunc_len* after left/right trimming are discarded.
        Set to 0 to disable.
    trunc_q : int or None
        Truncate at the first quality score ``<=`` this value.  ``None``
        disables quality truncation.
    max_len : int
        Discard reads longer than this **before** trimming.  0 disables.
    min_len : int
        Discard reads shorter than this **after** all trimming.
    max_n : int
        Maximum number of ambiguous (non-ACGT) bases allowed.
    min_q : int
        Discard reads containing any quality score below this value.
    max_ee : float
        Maximum expected errors (``sum(10^(-Q/10))``).
    rm_phix : bool
        Remove reads matching the PhiX genome.
    rm_lowcomplex : float
        Remove reads with sequence complexity below this threshold.
        0 disables.
    orient_fwd : str or None
        If set, orient reads so they begin with this primer sequence.
        Reads that match neither in the forward nor reverse-complement
        orientation are discarded.
    compress : bool
        Gzip-compress the output file.
    verbose : bool
        Print a summary line when done.

    Returns
    -------
    tuple of (int, int)
        ``(reads_in, reads_out)``
    """
    reads_in = 0
    reads_out = 0

    # Ensure output directory exists
    out_dir = os.path.dirname(fout)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with _open_fq(fn) as fin, _write_opener(fout, compress) as fout_h:
        for header, seq, sep, qual_str in _fastq_records(fin):
            reads_in += 1

            # --- orient_fwd ---
            if orient_fwd is not None:
                result = _orient_read(seq, qual_str, orient_fwd)
                if result is None:
                    continue
                seq, qual_str = result

            # --- max_len (before trimming) ---
            if max_len > 0 and len(seq) > max_len:
                continue

            # --- trunc_q (truncate at first quality <= threshold) ---
            if trunc_q is not None:
                quals = _decode_quals(qual_str)
                trunc_pos = len(quals)
                for i, q in enumerate(quals):
                    if q <= trunc_q:
                        trunc_pos = i
                        break
                seq = seq[:trunc_pos]
                qual_str = qual_str[:trunc_pos]

            # --- trim_left / trim_right ---
            if trim_left > 0:
                seq = seq[trim_left:]
                qual_str = qual_str[trim_left:]
            if trim_right > 0 and trim_right < len(seq):
                seq = seq[: len(seq) - trim_right]
                qual_str = qual_str[: len(qual_str) - trim_right]
            elif trim_right > 0:
                seq = ""
                qual_str = ""

            # --- trunc_len ---
            if trunc_len > 0:
                if len(seq) < trunc_len:
                    continue
                seq = seq[:trunc_len]
                qual_str = qual_str[:trunc_len]

            # --- min_len (after all trimming) ---
            if len(seq) < min_len:
                continue

            # --- max_n ---
            if max_n >= 0 and _count_non_acgt(seq) > max_n:
                continue

            # Decode quality scores for remaining checks
            quals = _decode_quals(qual_str)

            # --- min_q ---
            if min_q > 0 and quals and min(quals) < min_q:
                continue

            # --- max_ee ---
            if max_ee < float("inf"):
                if _expected_errors(quals) > max_ee:
                    continue

            # --- rm_phix ---
            if rm_phix:
                phix_match = _is_phix([seq])
                if phix_match[0]:
                    continue

            # --- rm_lowcomplex ---
            if rm_lowcomplex > 0:
                complexity = _seq_complexity([seq])
                if complexity[0] < rm_lowcomplex:
                    continue

            # Write passing read
            fout_h.write(f"{header}\n{seq}\n{sep}\n{qual_str}\n")
            reads_out += 1

    if verbose:
        pct = (reads_out / reads_in * 100) if reads_in > 0 else 0.0
        print(
            f"Read in {reads_in} reads, output {reads_out} "
            f"({pct:.1f}%) filtered reads."
        )

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
    trunc_q: Union[None, int, Tuple[Optional[int], Optional[int]]] = (None, None),
    max_len: Union[int, Tuple[int, int]] = (0, 0),
    min_len: Union[int, Tuple[int, int]] = (20, 20),
    max_n: Union[int, Tuple[int, int]] = (0, 0),
    min_q: Union[int, Tuple[int, int]] = (0, 0),
    max_ee: Union[float, Tuple[float, float]] = (float("inf"), float("inf")),
    rm_phix: bool = True,
    rm_lowcomplex: Union[float, Tuple[float, float]] = (0.0, 0.0),
    orient_fwd: Optional[str] = None,
    compress: bool = True,
    verbose: bool = False,
) -> Tuple[int, int]:
    """Filter and trim paired FASTQ files.

    Both reads in a pair must pass **all** filters for either to be kept,
    preserving read-pair synchronisation.

    Parameters accept scalars (applied to both reads) or length-2 tuples
    ``(forward_value, reverse_value)``.

    Parameters
    ----------
    fwd : str
        Path to forward-read input FASTQ.
    filt_fwd : str
        Path to filtered forward-read output FASTQ.
    rev : str
        Path to reverse-read input FASTQ.
    filt_rev : str
        Path to filtered reverse-read output FASTQ.
    trim_left, trim_right, trunc_len, trunc_q, max_len, min_len,
    max_n, min_q, max_ee, rm_lowcomplex
        See :func:`fastq_filter` for per-parameter documentation.  Each
        accepts a scalar or a ``(fwd, rev)`` tuple.
    rm_phix : bool
        Remove reads matching PhiX.
    orient_fwd : str or None
        If set, orient forward reads to start with this sequence (reverse
        reads are reverse-complemented accordingly).
    compress : bool
        Gzip-compress output files.
    verbose : bool
        Print a summary.

    Returns
    -------
    tuple of (int, int)
        ``(reads_in, reads_out)``
    """

    def _as_pair(val, default=None):
        """Normalise a scalar or tuple to a 2-tuple."""
        if isinstance(val, (list, tuple)):
            if len(val) == 2:
                return tuple(val)
            raise ValueError(f"Expected length-2 tuple, got {len(val)}")
        return (val, val)

    tl = _as_pair(trim_left)
    tr = _as_pair(trim_right)
    tlen = _as_pair(trunc_len)
    tq = _as_pair(trunc_q)
    mxl = _as_pair(max_len)
    mnl = _as_pair(min_len)
    mn = _as_pair(max_n)
    mq = _as_pair(min_q)
    mee = _as_pair(max_ee)
    rlc = _as_pair(rm_lowcomplex)

    reads_in = 0
    reads_out = 0

    # Ensure output directories exist
    for p in (filt_fwd, filt_rev):
        d = os.path.dirname(p)
        if d:
            os.makedirs(d, exist_ok=True)

    with (
        _open_fq(fwd) as fin_f,
        _open_fq(rev) as fin_r,
        _write_opener(filt_fwd, compress) as fout_f,
        _write_opener(filt_rev, compress) as fout_r,
    ):
        for (hdr_f, seq_f, sep_f, qual_f), (hdr_r, seq_r, sep_r, qual_r) in zip(
            _fastq_records(fin_f), _fastq_records(fin_r)
        ):
            reads_in += 1

            seqs = [seq_f, seq_r]
            quals_str = [qual_f, qual_r]
            keep = True

            # --- orient_fwd (forward read only; RC the reverse accordingly) ---
            if orient_fwd is not None:
                result = _orient_read(seqs[0], quals_str[0], orient_fwd)
                if result is None:
                    continue
                if result != (seqs[0], quals_str[0]):
                    # Forward was RC'd -- also RC reverse to maintain pairing
                    seqs[0], quals_str[0] = result
                    seqs[1] = _rc(seqs[1])
                    quals_str[1] = quals_str[1][::-1]

            for idx in range(2):
                if not keep:
                    break

                s = seqs[idx]
                qs = quals_str[idx]

                # --- max_len (before trimming) ---
                if mxl[idx] > 0 and len(s) > mxl[idx]:
                    keep = False
                    break

                # --- trunc_q ---
                if tq[idx] is not None:
                    decoded = _decode_quals(qs)
                    trunc_pos = len(decoded)
                    for i, q in enumerate(decoded):
                        if q <= tq[idx]:
                            trunc_pos = i
                            break
                    s = s[:trunc_pos]
                    qs = qs[:trunc_pos]

                # --- trim_left / trim_right ---
                if tl[idx] > 0:
                    s = s[tl[idx]:]
                    qs = qs[tl[idx]:]
                if tr[idx] > 0 and tr[idx] < len(s):
                    s = s[: len(s) - tr[idx]]
                    qs = qs[: len(qs) - tr[idx]]
                elif tr[idx] > 0:
                    s = ""
                    qs = ""

                # --- trunc_len ---
                if tlen[idx] > 0:
                    if len(s) < tlen[idx]:
                        keep = False
                        break
                    s = s[: tlen[idx]]
                    qs = qs[: tlen[idx]]

                # --- min_len ---
                if len(s) < mnl[idx]:
                    keep = False
                    break

                # --- max_n ---
                if mn[idx] >= 0 and _count_non_acgt(s) > mn[idx]:
                    keep = False
                    break

                decoded = _decode_quals(qs)

                # --- min_q ---
                if mq[idx] > 0 and decoded and min(decoded) < mq[idx]:
                    keep = False
                    break

                # --- max_ee ---
                if mee[idx] < float("inf"):
                    if _expected_errors(decoded) > mee[idx]:
                        keep = False
                        break

                # --- rm_lowcomplex ---
                if rlc[idx] > 0:
                    complexity = _seq_complexity([s])
                    if complexity[0] < rlc[idx]:
                        keep = False
                        break

                seqs[idx] = s
                quals_str[idx] = qs

            if not keep:
                continue

            # --- rm_phix (check both reads) ---
            if rm_phix:
                phix_flags = _is_phix(seqs)
                if phix_flags.any():
                    continue

            # Write passing pair
            fout_f.write(f"{hdr_f}\n{seqs[0]}\n{sep_f}\n{quals_str[0]}\n")
            fout_r.write(f"{hdr_r}\n{seqs[1]}\n{sep_r}\n{quals_str[1]}\n")
            reads_out += 1

    if verbose:
        pct = (reads_out / reads_in * 100) if reads_in > 0 else 0.0
        print(
            f"Read in {reads_in} paired-reads, output {reads_out} "
            f"({pct:.1f}%) filtered paired-reads."
        )

    return reads_in, reads_out


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
    trunc_q: Union[None, int, Tuple[Optional[int], Optional[int]]] = None,
    max_len: Union[int, Tuple[int, int]] = 0,
    min_len: Union[int, Tuple[int, int]] = 20,
    max_n: Union[int, Tuple[int, int]] = 0,
    min_q: Union[int, Tuple[int, int]] = 0,
    max_ee: Union[float, Tuple[float, float]] = float("inf"),
    rm_phix: bool = True,
    rm_lowcomplex: Union[float, Tuple[float, float]] = 0.0,
    orient_fwd: Optional[str] = None,
    compress: bool = True,
    multithread: Union[bool, int] = False,
    verbose: bool = False,
) -> np.ndarray:
    """Filter and trim FASTQ files (single- or paired-end).

    This is a convenience wrapper around :func:`fastq_filter` (single-end)
    and :func:`fastq_paired_filter` (paired-end).  It dispatches to the
    appropriate function based on whether *rev* is provided and optionally
    parallelises across files with :class:`~concurrent.futures.ProcessPoolExecutor`.

    Parameters
    ----------
    fwd : str or list of str
        Path(s) to forward (or single-end) input FASTQ file(s).
    filt : str or list of str
        Path(s) to filtered output FASTQ file(s).  Must be the same length
        as *fwd*.
    rev : str, list of str, or None
        Path(s) to reverse-read input FASTQ file(s).  ``None`` for
        single-end mode.
    filt_rev : str, list of str, or None
        Path(s) to filtered reverse-read output FASTQ file(s).  Required
        when *rev* is not ``None``.
    trim_left, trim_right, trunc_len, trunc_q, max_len, min_len,
    max_n, min_q, max_ee, rm_phix, rm_lowcomplex, orient_fwd, compress
        See :func:`fastq_filter` and :func:`fastq_paired_filter`.
    multithread : bool or int
        If ``True``, use all available cores.  If an integer > 1, use
        that many worker processes.  ``False`` or ``1`` disables
        parallelism.
    verbose : bool
        Print per-file summaries.

    Returns
    -------
    numpy.ndarray
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
            f"the same length."
        )

    paired = rev is not None
    if paired:
        if filt_rev is None:
            raise ValueError("filt_rev is required when rev is provided.")
        if len(rev) != n_files:
            raise ValueError(
                f"fwd ({n_files}) and rev ({len(rev)}) must have the same length."
            )
        if len(filt_rev) != n_files:
            raise ValueError(
                f"fwd ({n_files}) and filt_rev ({len(filt_rev)}) must have "
                f"the same length."
            )

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
        verbose=verbose,
    )

    results = np.zeros((n_files, 2), dtype=np.int64)

    if paired:

        def _run_paired(idx):
            return fastq_paired_filter(
                fwd[idx], filt[idx], rev[idx], filt_rev[idx], **common_kw
            )

        if n_workers > 1 and n_files > 1:
            with ProcessPoolExecutor(max_workers=n_workers) as pool:
                futures = {pool.submit(_run_paired, i): i for i in range(n_files)}
                for fut in futures:
                    i = futures[fut]
                    results[i] = fut.result()
        else:
            for i in range(n_files):
                results[i] = _run_paired(i)

    else:

        def _run_single(idx):
            # For single-end, normalise scalar params
            kw = dict(common_kw)
            # Extract single-end values from any tuple params
            for key in (
                "trim_left", "trim_right", "trunc_len", "trunc_q",
                "max_len", "min_len", "max_n", "min_q", "max_ee",
                "rm_lowcomplex",
            ):
                v = kw[key]
                if isinstance(v, (list, tuple)):
                    kw[key] = v[0]
            return fastq_filter(fwd[idx], filt[idx], **kw)

        if n_workers > 1 and n_files > 1:
            with ProcessPoolExecutor(max_workers=n_workers) as pool:
                futures = {pool.submit(_run_single, i): i for i in range(n_files)}
                for fut in futures:
                    i = futures[fut]
                    results[i] = fut.result()
        else:
            for i in range(n_files):
                results[i] = _run_single(i)

    return results
