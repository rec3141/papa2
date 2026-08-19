"""FASTQ reading, filtering, and dereplication."""

import gzip
import os
import ctypes as ct
import numpy as np
from concurrent.futures import ProcessPoolExecutor

# Try to load C dereplication from libpapa2.so
_lib_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "libpapa2.so")
try:
    _lib = ct.CDLL(_lib_path)

    class _DerepResult(ct.Structure):
        _fields_ = [
            ("n_uniques", ct.c_int),
            ("n_reads", ct.c_int),
            ("max_seq_len", ct.c_int),
            ("seqs", ct.POINTER(ct.c_char_p)),
            ("abundances", ct.POINTER(ct.c_int)),
            ("quals", ct.POINTER(ct.c_double)),
            ("map", ct.POINTER(ct.c_int)),
        ]

    _lib.derep_fastq_c.restype = ct.POINTER(_DerepResult)
    _lib.derep_fastq_c.argtypes = [ct.c_char_p, ct.c_int]
    _lib.derep_result_free.restype = None
    _lib.derep_result_free.argtypes = [ct.POINTER(_DerepResult)]
    _HAS_C_DEREP = True
except (OSError, AttributeError):
    _HAS_C_DEREP = False


def _derep_fastq_c(filepath, qual_offset=-1):
    """Fast C dereplication via zlib."""
    path_bytes = filepath.encode('utf-8') if isinstance(filepath, str) else filepath
    res_ptr = _lib.derep_fastq_c(path_bytes, ct.c_int(qual_offset))
    if not res_ptr:
        raise RuntimeError(f"Failed to derep {filepath}")
    res = res_ptr.contents
    nu = res.n_uniques
    ml = res.max_seq_len

    seqs = [res.seqs[i].decode('ascii') for i in range(nu)]
    abunds = np.ctypeslib.as_array(res.abundances, shape=(nu,)).copy()
    quals = np.ctypeslib.as_array(res.quals, shape=(nu * ml,)).copy().reshape(nu, ml)
    nr = res.n_reads
    rmap = np.ctypeslib.as_array(res.map, shape=(nr,)).copy()

    _lib.derep_result_free(res_ptr)
    return {"seqs": seqs, "abundances": abunds, "quals": quals, "map": rmap}


def _derep_one(args):
    filepath, verbose, quality_type = args
    return derep_fastq(filepath, verbose=verbose, quality_type=quality_type)


_QUAL_OFFSETS = {"Auto": -1, "FastqQuality": 33, "SFastqQuality": 64,
                 "SolexaQuality": 64}


def derep_fastq(filepath, verbose=False, with_map=False, multithread=True,
                quality_type="Auto"):
    """Dereplicate a FASTQ file (or a list of them).

    Uses C implementation (zlib) when available for ~2x speedup.
    Always returns the per-read map (read_idx -> unique_idx).
    The with_map parameter is accepted for backward compatibility but ignored.

    A list of paths returns a list of results; with multithread (default)
    the files are processed in a worker pool.

    quality_type follows R: "Auto" (default) detects Phred+33 vs Phred+64
    per file with ShortRead's rule; "FastqQuality" and "SFastqQuality"
    force the offset.

    Returns:
        dict with keys:
            seqs: list[str], unique sequences sorted by abundance (descending)
            abundances: numpy int32 array
            quals: numpy float64 array (n_uniques x max_seqlen), average quality
            map: numpy int32 array, maps each read to its unique index (0-indexed)
    """
    if isinstance(filepath, (list, tuple)):
        files = list(filepath)
        if multithread is True:
            n_workers = min(len(files), os.cpu_count() or 1)
        elif isinstance(multithread, int) and multithread > 1:
            n_workers = min(len(files), multithread)
        else:
            n_workers = 1
        tasks = [(f, verbose, quality_type) for f in files]
        if n_workers > 1:
            try:
                from loky import get_reusable_executor
                with get_reusable_executor(max_workers=n_workers) as pool:
                    return list(pool.map(_derep_one, tasks))
            except ImportError:
                pass
        return [_derep_one(t) for t in tasks]

    if quality_type not in _QUAL_OFFSETS:
        raise ValueError(f"Unknown quality_type: {quality_type!r}")
    offset = _QUAL_OFFSETS[quality_type]

    # Use C implementation if available
    if _HAS_C_DEREP:
        result = _derep_fastq_c(filepath, offset)
        if verbose:
            print(f"Read {result['abundances'].sum()} reads, {len(result['seqs'])} unique sequences")
        return result

    # Fallback: Python implementation
    opener = gzip.open if filepath.endswith(".gz") else open
    with opener(filepath, "rb") as f:
        raw = f.read()

    # Split into lines, extract every 4th (seq) and every 4th+3 (qual)
    lines = raw.split(b'\n')
    # Remove trailing empty line
    if lines and lines[-1] == b'':
        lines.pop()
    n_reads = len(lines) // 4

    # Phase 1: Build dedup index from sequences
    seq_to_idx = {}
    first_seen = []
    counts_list = []
    read_uid = np.empty(n_reads, dtype=np.int32)

    for i in range(n_reads):
        seq = lines[i * 4 + 1].upper()  # bytes
        idx = seq_to_idx.get(seq)
        if idx is None:
            idx = len(first_seen)
            seq_to_idx[seq] = idx
            first_seen.append(seq)
            counts_list.append(0)
        counts_list[idx] += 1
        read_uid[i] = idx

    n_uniques = len(first_seen)
    counts = np.array(counts_list, dtype=np.int32)
    maxlen = max(len(s) for s in first_seen) if first_seen else 0

    # Phase 2: Accumulate quality scores using numpy frombuffer
    _py_offset = float(offset)
    if offset < 0:
        lo, hi = 255, 0
        for i in range(min(n_reads, 10000)):
            ql = lines[i * 4 + 3]
            if ql:
                lo = min(lo, min(ql))
                hi = max(hi, max(ql))
        _py_offset = 64.0 if (hi and lo >= 59 and hi >= 75) else 33.0

    qual_sums = np.zeros((n_uniques, maxlen), dtype=np.float64)

    for i in range(n_reads):
        uid = read_uid[i]
        qline = lines[i * 4 + 3]
        q = np.frombuffer(qline, dtype=np.uint8).astype(np.float64)
        q -= _py_offset
        slen = len(q)
        qual_sums[uid, :slen] += q

    # Phase 3: Sort by abundance descending
    sort_idx = np.argsort(-counts)
    sorted_seqs = [first_seen[i].decode('ascii') for i in sort_idx]
    sorted_counts = counts[sort_idx]
    sorted_quals = qual_sums[sort_idx]

    # Average and NaN-pad
    for i in range(n_uniques):
        slen = len(sorted_seqs[i])
        if sorted_counts[i] > 0:
            sorted_quals[i, :slen] /= sorted_counts[i]
        sorted_quals[i, slen:] = np.nan

    # Remap read indices
    inv_sort = np.empty(n_uniques, dtype=np.int32)
    inv_sort[sort_idx] = np.arange(n_uniques, dtype=np.int32)
    rmap = inv_sort[read_uid]

    if verbose:
        print(f"Read {n_reads} reads, {n_uniques} unique sequences")

    return {
        "seqs": sorted_seqs,
        "abundances": sorted_counts,
        "quals": sorted_quals,
        "map": rmap,
    }


def combine_dereps(dereps):
    """Combine multiple derep results into one pooled derep.

    Exact port of R's dada2:::combineDereps2: unique sequences keep
    first-appearance order across samples, counts are summed, quality
    profiles are abundance-weighted means accumulated in sample order,
    the combined uniques are re-sorted by decreasing abundance (stable,
    like R's order()), and the concatenated per-read map is remapped
    accordingly.
    """
    if isinstance(dereps, dict):
        dereps = [dereps]
    maxlen = max(d["quals"].shape[1] for d in dereps)

    seq_index = {}
    seqs_all = []
    for d in dereps:
        for sq in d["seqs"]:
            if sq not in seq_index:
                seq_index[sq] = len(seqs_all)
                seqs_all.append(sq)
    n = len(seqs_all)

    counts = np.zeros(n, dtype=np.int64)
    qual_sums = np.zeros((n, maxlen), dtype=np.float64)
    maps = []
    for d in dereps:
        idx = np.fromiter((seq_index[sq] for sq in d["seqs"]),
                          dtype=np.int64, count=len(d["seqs"]))
        ab = np.asarray(d["abundances"], dtype=np.int64)
        counts[idx] += ab
        q = d["quals"]
        if q.shape[1] < maxlen:
            pad = np.full((q.shape[0], maxlen - q.shape[1]), np.nan)
            q = np.hstack([q, pad])
        qual_sums[idx] += q * ab[:, None]
        maps.append(idx[np.asarray(d["map"], dtype=np.int64)])

    quals = qual_sums / counts[:, None]

    order = np.argsort(-counts, kind="stable")
    inv = np.empty(n, dtype=np.int64)
    inv[order] = np.arange(n, dtype=np.int64)
    combined_map = inv[np.concatenate(maps)] if maps else np.zeros(0, np.int64)

    return {
        "seqs": [seqs_all[i] for i in order],
        "abundances": counts[order].astype(np.int32),
        "quals": quals[order],
        "map": combined_map.astype(np.int32),
    }
