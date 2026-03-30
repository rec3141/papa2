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
    _lib.derep_fastq_c.argtypes = [ct.c_char_p]
    _lib.derep_result_free.restype = None
    _lib.derep_result_free.argtypes = [ct.POINTER(_DerepResult)]
    _HAS_C_DEREP = True
except (OSError, AttributeError):
    _HAS_C_DEREP = False


def _derep_fastq_c(filepath):
    """Fast C dereplication via zlib."""
    path_bytes = filepath.encode('utf-8') if isinstance(filepath, str) else filepath
    res_ptr = _lib.derep_fastq_c(path_bytes)
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


def derep_fastq(filepath, verbose=False, with_map=False):
    """Dereplicate a FASTQ file.

    Uses C implementation (zlib) when available for ~2x speedup.
    Always returns the per-read map (read_idx -> unique_idx).
    The with_map parameter is accepted for backward compatibility but ignored.

    Returns:
        dict with keys:
            seqs: list[str], unique sequences sorted by abundance (descending)
            abundances: numpy int32 array
            quals: numpy float64 array (n_uniques x max_seqlen), average quality
            map: numpy int32 array, maps each read to its unique index (0-indexed)
    """
    # Use C implementation if available
    if _HAS_C_DEREP:
        result = _derep_fastq_c(filepath)
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
    qual_sums = np.zeros((n_uniques, maxlen), dtype=np.float64)

    for i in range(n_reads):
        uid = read_uid[i]
        qline = lines[i * 4 + 3]
        q = np.frombuffer(qline, dtype=np.uint8).astype(np.float64)
        q -= 33.0
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
