"""ctypes bindings to libpapa2.so C API."""

import ctypes as ct
import numpy as np
import os

# Load the shared library
_lib_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "libpapa2.so")
if not os.path.exists(_lib_path):
    _lib_path = "libpapa2.so"
_lib = ct.CDLL(_lib_path)


class DadaResult(ct.Structure):
    _fields_ = [
        ("nclust", ct.c_int),
        ("nraw", ct.c_int),
        ("maxlen", ct.c_int),
        ("cluster_seqs", ct.POINTER(ct.c_char_p)),
        ("cluster_abunds", ct.POINTER(ct.c_int)),
        ("cluster_n0", ct.POINTER(ct.c_int)),
        ("cluster_n1", ct.POINTER(ct.c_int)),
        ("cluster_nunq", ct.POINTER(ct.c_int)),
        ("cluster_pval", ct.POINTER(ct.c_double)),
        ("ncol_trans", ct.c_int),
        ("trans", ct.POINTER(ct.c_int)),
        ("map", ct.POINTER(ct.c_int)),
        ("pval", ct.POINTER(ct.c_double)),
    ]


_lib.dada2_run.restype = ct.POINTER(DadaResult)
_lib.dada2_run.argtypes = [
    ct.POINTER(ct.c_char_p),  # seqs
    ct.POINTER(ct.c_int),     # abundances
    ct.POINTER(ct.c_int),     # priors
    ct.c_int,                  # nraw
    ct.POINTER(ct.c_double),  # err_mat
    ct.c_int,                  # ncol_err
    ct.POINTER(ct.c_double),  # quals
    ct.c_int,                  # maxlen
    ct.c_int, ct.c_int, ct.c_int,  # match, mismatch, gap_pen
    ct.c_int, ct.c_double, ct.c_int,  # use_kmers, kdist_cutoff, band_size
    ct.c_double, ct.c_double, ct.c_double,  # omegaA, omegaP, omegaC
    ct.c_int, ct.c_int,       # detect_singletons, max_clust
    ct.c_double, ct.c_int, ct.c_int,  # min_fold, min_hamming, min_abund
    ct.c_int, ct.c_int, ct.c_int,  # use_quals, vectorized_alignment, homo_gap_pen
    ct.c_int, ct.c_int,       # multithread, verbose
    ct.c_int, ct.c_int, ct.c_int,  # SSE, gapless, greedy
]

_lib.dada2_result_free.restype = None
_lib.dada2_result_free.argtypes = [ct.POINTER(DadaResult)]

def run_dada(seqs, abundances, err_mat, quals=None, priors=None,
             match=5, mismatch=-4, gap_pen=-8,
             use_kmers=True, kdist_cutoff=0.42, band_size=16,
             omega_a=1e-40, omega_p=1e-4, omega_c=1e-40,
             detect_singletons=False, max_clust=0,
             min_fold=1, min_hamming=1, min_abund=1,
             use_quals=True, vectorized_alignment=True, homo_gap_pen=-8,
             multithread=True, verbose=False,
             sse=2, gapless=True, greedy=True):
    """Call the C++ dada2 algorithm on dereplicated sequences.

    Args:
        seqs: list of str, unique DNA sequences (ACGT only)
        abundances: array-like of int, abundance per unique sequence
        err_mat: numpy array (16, ncol), error rate matrix, row-major
        quals: numpy array (nraw, maxlen) of avg quality scores, or None
        priors: array-like of int (0/1), or None

    Returns:
        dict with keys: cluster_seqs, cluster_abunds, trans, map, pval, etc.
    """
    nraw = len(seqs)
    if nraw == 0:
        return {"cluster_seqs": [], "cluster_abunds": np.array([], dtype=np.int32),
                "trans": np.zeros((16, 0), dtype=np.int32), "map": np.array([], dtype=np.int32),
                "pval": np.array([], dtype=np.float64)}

    # Convert sequences to C strings
    seq_arr = (ct.c_char_p * nraw)()
    for i, s in enumerate(seqs):
        seq_arr[i] = s.encode("ascii") if isinstance(s, str) else s

    # Abundances
    abund_arr = np.asarray(abundances, dtype=np.int32)

    # Priors
    if priors is None:
        prior_arr = np.zeros(nraw, dtype=np.int32)
    else:
        prior_arr = np.asarray(priors, dtype=np.int32)

    # Error matrix (16 x ncol, row-major, C-contiguous)
    err = np.ascontiguousarray(err_mat, dtype=np.float64)
    ncol_err = err.shape[1]

    # Quality matrix (C++ rounds to uint8_t internally via round())
    if quals is not None:
        q = np.ascontiguousarray(quals, dtype=np.float64)
        maxlen = q.shape[1]
        q_ptr = q.ctypes.data_as(ct.POINTER(ct.c_double))
    else:
        maxlen = max(len(s) for s in seqs)
        q_ptr = ct.POINTER(ct.c_double)()

    res_ptr = _lib.dada2_run(
        seq_arr,
        abund_arr.ctypes.data_as(ct.POINTER(ct.c_int)),
        prior_arr.ctypes.data_as(ct.POINTER(ct.c_int)),
        ct.c_int(nraw),
        err.ctypes.data_as(ct.POINTER(ct.c_double)),
        ct.c_int(ncol_err),
        q_ptr,
        ct.c_int(maxlen),
        ct.c_int(match), ct.c_int(mismatch), ct.c_int(gap_pen),
        ct.c_int(int(use_kmers)), ct.c_double(kdist_cutoff), ct.c_int(band_size),
        ct.c_double(omega_a), ct.c_double(omega_p), ct.c_double(omega_c),
        ct.c_int(int(detect_singletons)), ct.c_int(max_clust),
        ct.c_double(min_fold), ct.c_int(min_hamming), ct.c_int(min_abund),
        ct.c_int(int(use_quals)), ct.c_int(int(vectorized_alignment)),
        ct.c_int(homo_gap_pen),
        ct.c_int(int(multithread)), ct.c_int(int(verbose)),
        ct.c_int(sse), ct.c_int(int(gapless)), ct.c_int(int(greedy)),
    )

    if not res_ptr:
        raise RuntimeError("dada2_run returned NULL")

    res = res_ptr.contents
    nc = res.nclust
    nr = res.nraw

    # Extract results using fast buffer copies
    result = {
        "cluster_seqs": [res.cluster_seqs[i].decode("ascii") for i in range(nc)],
        "cluster_abunds": np.ctypeslib.as_array(res.cluster_abunds, shape=(nc,)).copy(),
        "cluster_n0": np.ctypeslib.as_array(res.cluster_n0, shape=(nc,)).copy(),
        "cluster_n1": np.ctypeslib.as_array(res.cluster_n1, shape=(nc,)).copy(),
        "cluster_nunq": np.ctypeslib.as_array(res.cluster_nunq, shape=(nc,)).copy(),
        "trans": np.ctypeslib.as_array(res.trans, shape=(16 * res.ncol_trans,)).copy().reshape(16, res.ncol_trans) if res.ncol_trans > 0 else np.zeros((16, 0), dtype=np.int32),
        "map": np.ctypeslib.as_array(res.map, shape=(nr,)).copy(),
        "pval": np.ctypeslib.as_array(res.pval, shape=(nr,)).copy(),
    }

    _lib.dada2_result_free(res_ptr)
    return result


# =========================================================================
# Taxonomy assignment
# =========================================================================

class TaxResult(ct.Structure):
    _fields_ = [
        ("nseq", ct.c_int),
        ("nlevel", ct.c_int),
        ("rval", ct.POINTER(ct.c_int)),
        ("rboot", ct.POINTER(ct.c_int)),
    ]


_lib.dada2_assign_taxonomy.restype = ct.POINTER(TaxResult)
_lib.dada2_assign_taxonomy.argtypes = [
    ct.POINTER(ct.c_char_p),  # seqs
    ct.c_int,                  # nseq
    ct.POINTER(ct.c_char_p),  # refs
    ct.c_int,                  # nref
    ct.POINTER(ct.c_int),     # ref_to_genus
    ct.POINTER(ct.c_int),     # genusmat
    ct.c_int,                  # ngenus
    ct.c_int,                  # nlevel
    ct.c_int,                  # verbose
]

_lib.dada2_tax_result_free.restype = None
_lib.dada2_tax_result_free.argtypes = [ct.POINTER(TaxResult)]


# =========================================================================
# Paired-read merging functions
# =========================================================================

_lib.dada2_nwalign.restype = ct.c_int
_lib.dada2_nwalign.argtypes = [
    ct.c_char_p,   # s1
    ct.c_char_p,   # s2
    ct.c_int,      # match
    ct.c_int,      # mismatch
    ct.c_int,      # gap_p
    ct.c_int,      # band
    ct.POINTER(ct.c_void_p),  # al1_out
    ct.POINTER(ct.c_void_p),  # al2_out
]

_lib.dada2_eval_pair.restype = None
_lib.dada2_eval_pair.argtypes = [
    ct.c_char_p,              # al1
    ct.c_char_p,              # al2
    ct.POINTER(ct.c_int),     # out_match
    ct.POINTER(ct.c_int),     # out_mismatch
    ct.POINTER(ct.c_int),     # out_indel
]

# Use c_void_p for returned malloc'd strings to prevent ctypes auto-free
_lib.dada2_pair_consensus.restype = ct.c_void_p
_lib.dada2_pair_consensus.argtypes = [
    ct.c_char_p,   # al1
    ct.c_char_p,   # al2
    ct.c_int,      # prefer
    ct.c_int,      # trim_overhang
]

_lib.dada2_rc.restype = ct.c_void_p
_lib.dada2_rc.argtypes = [ct.c_char_p]

_lib.dada2_free_string.restype = None
_lib.dada2_free_string.argtypes = [ct.c_void_p]


def nwalign(s1, s2, match=5, mismatch=-4, gap_p=-8, band=-1):
    """NW ends-free alignment of two ACGT strings.

    Returns:
        (al1, al2): tuple of aligned strings.
    """
    al1_p = ct.c_void_p()
    al2_p = ct.c_void_p()
    s1_b = s1.encode("ascii") if isinstance(s1, str) else s1
    s2_b = s2.encode("ascii") if isinstance(s2, str) else s2

    ret = _lib.dada2_nwalign(s1_b, s2_b, match, mismatch, gap_p, band,
                              ct.byref(al1_p), ct.byref(al2_p))
    if ret != 0:
        raise RuntimeError("dada2_nwalign failed")

    al1 = ct.cast(al1_p, ct.c_char_p).value.decode("ascii")
    al2 = ct.cast(al2_p, ct.c_char_p).value.decode("ascii")
    _lib.dada2_free_string(al1_p)
    _lib.dada2_free_string(al2_p)
    return al1, al2


def eval_pair(al1, al2):
    """Evaluate an alignment: count matches, mismatches, indels (skipping end gaps).

    Returns:
        (nmatch, nmismatch, nindel): tuple of ints.
    """
    m = ct.c_int(0)
    mm = ct.c_int(0)
    ind = ct.c_int(0)
    al1_b = al1.encode("ascii") if isinstance(al1, str) else al1
    al2_b = al2.encode("ascii") if isinstance(al2, str) else al2

    _lib.dada2_eval_pair(al1_b, al2_b, ct.byref(m), ct.byref(mm), ct.byref(ind))
    return m.value, mm.value, ind.value


def pair_consensus(al1, al2, prefer=1, trim_overhang=True):
    """Build consensus from two aligned strings.

    Args:
        prefer: 1 = al1 wins mismatches, 2 = al2 wins.
        trim_overhang: if True, trim overhanging ends.

    Returns:
        consensus string.
    """
    al1_b = al1.encode("ascii") if isinstance(al1, str) else al1
    al2_b = al2.encode("ascii") if isinstance(al2, str) else al2

    result_p = _lib.dada2_pair_consensus(al1_b, al2_b, prefer, int(trim_overhang))
    if not result_p:
        raise RuntimeError("dada2_pair_consensus failed")
    result = ct.cast(result_p, ct.c_char_p).value.decode("ascii")
    _lib.dada2_free_string(result_p)
    return result


def rc(seq):
    """Reverse complement an ACGT string."""
    seq_b = seq.encode("ascii") if isinstance(seq, str) else seq
    result_p = _lib.dada2_rc(seq_b)
    if not result_p:
        raise RuntimeError("dada2_rc failed")
    result = ct.cast(result_p, ct.c_char_p).value.decode("ascii")
    _lib.dada2_free_string(result_p)
    return result


def run_taxonomy(seqs, refs, ref_to_genus, genusmat, ngenus, nlevel, verbose=True):
    """Run dada2 taxonomy assignment via C library.

    Args:
        seqs: list of query sequences (str)
        refs: list of reference sequences (str)
        ref_to_genus: numpy array (nref,) int32, 0-indexed genus ID per ref
        genusmat: numpy array (ngenus, nlevel) int32, genus-to-level mapping
        ngenus: int
        nlevel: int
        verbose: bool

    Returns:
        dict with:
            rval: numpy array (nseq,) int32, 1-indexed best genus per query (0=NA)
            rboot: numpy array (nseq, nlevel) int32, bootstrap counts
    """
    nseq = len(seqs)
    nref = len(refs)

    # Convert sequences to C strings
    seq_arr = (ct.c_char_p * nseq)()
    for i, s in enumerate(seqs):
        seq_arr[i] = s.encode("ascii") if isinstance(s, str) else s

    ref_arr = (ct.c_char_p * nref)()
    for i, s in enumerate(refs):
        ref_arr[i] = s.encode("ascii") if isinstance(s, str) else s

    # Reference-to-genus mapping
    rtg = np.ascontiguousarray(ref_to_genus, dtype=np.int32)

    # Genus matrix (row-major)
    gmat = np.ascontiguousarray(genusmat, dtype=np.int32)

    res_ptr = _lib.dada2_assign_taxonomy(
        seq_arr, nseq,
        ref_arr, nref,
        rtg.ctypes.data_as(ct.POINTER(ct.c_int)),
        gmat.ctypes.data_as(ct.POINTER(ct.c_int)),
        ngenus, nlevel,
        ct.c_int(int(verbose))
    )

    if not res_ptr:
        raise RuntimeError("dada2_assign_taxonomy returned NULL")

    res = res_ptr.contents
    result = {
        "rval": np.ctypeslib.as_array(res.rval, shape=(nseq,)).copy(),
        "rboot": np.ctypeslib.as_array(res.rboot, shape=(nseq * nlevel,)).copy().reshape(nseq, nlevel),
    }

    _lib.dada2_tax_result_free(res_ptr)
    return result


# =========================================================================
# Chimera detection
# =========================================================================

class ChimeraResult(ct.Structure):
    _fields_ = [
        ("n_seqs", ct.c_int),
        ("nflag", ct.POINTER(ct.c_int)),
        ("nsam", ct.POINTER(ct.c_int)),
    ]


_lib.dada2_is_bimera.restype = ct.c_int
_lib.dada2_is_bimera.argtypes = [
    ct.c_char_p,               # seq
    ct.POINTER(ct.c_char_p),   # parents
    ct.c_int,                   # n_parents
    ct.c_int,                   # allow_one_off
    ct.c_int,                   # min_one_off_par_dist
    ct.c_int, ct.c_int, ct.c_int, ct.c_int,  # match, mismatch, gap_p, max_shift
]

_lib.dada2_table_bimera.restype = ct.POINTER(ChimeraResult)
_lib.dada2_table_bimera.argtypes = [
    ct.POINTER(ct.c_int),      # mat
    ct.c_int, ct.c_int,        # nrow, ncol
    ct.POINTER(ct.c_char_p),   # seqs
    ct.c_double,                # min_fold
    ct.c_int,                   # min_abund
    ct.c_int,                   # allow_one_off
    ct.c_int,                   # min_one_off_par_dist
    ct.c_int, ct.c_int, ct.c_int, ct.c_int,  # match, mismatch, gap_p, max_shift
]

_lib.dada2_chimera_result_free.restype = None
_lib.dada2_chimera_result_free.argtypes = [ct.POINTER(ChimeraResult)]


def is_bimera(seq, parents, allow_one_off=False, min_one_off_par_dist=4,
              match=5, mismatch=-4, gap_p=-8, max_shift=16):
    """Check if seq is a bimera of the given parent sequences.

    Args:
        seq: query DNA sequence (str, ACGT).
        parents: list of parent DNA sequences (str, ACGT).
        allow_one_off: allow one mismatch in chimera model.
        min_one_off_par_dist: minimum hamming distance between parents for one-off.
        match, mismatch, gap_p, max_shift: alignment parameters.

    Returns:
        True if seq is a bimera, False otherwise.
    """
    n_parents = len(parents)
    if n_parents == 0:
        return False

    seq_b = seq.encode("ascii") if isinstance(seq, str) else seq
    par_arr = (ct.c_char_p * n_parents)()
    for i, p in enumerate(parents):
        par_arr[i] = p.encode("ascii") if isinstance(p, str) else p

    ret = _lib.dada2_is_bimera(
        seq_b, par_arr, n_parents,
        int(allow_one_off), min_one_off_par_dist,
        match, mismatch, gap_p, max_shift
    )
    return bool(ret)


def table_bimera(mat, seqs, min_fold=1.5, min_abund=2,
                 allow_one_off=False, min_one_off_par_dist=4,
                 match=5, mismatch=-4, gap_p=-8, max_shift=16):
    """Table-level consensus chimera detection.

    Args:
        mat: numpy int32 array (nrow x ncol), column-major (Fortran order).
             Rows = samples, columns = ASVs.
        seqs: list of ASV sequences (str, ACGT), length ncol.
        min_fold: parent must be this-fold more abundant.
        min_abund: parent minimum absolute abundance.
        allow_one_off: allow one mismatch in chimera model.
        min_one_off_par_dist: minimum hamming distance for one-off parents.
        match, mismatch, gap_p, max_shift: alignment parameters.

    Returns:
        dict with:
            nflag: numpy int32 array (ncol,) - per-ASV count of samples flagging chimeric.
            nsam: numpy int32 array (ncol,) - per-ASV count of samples where present.
    """
    mat = np.asfortranarray(mat, dtype=np.int32)
    nrow, ncol = mat.shape

    seq_arr = (ct.c_char_p * ncol)()
    for i, s in enumerate(seqs):
        seq_arr[i] = s.encode("ascii") if isinstance(s, str) else s

    res_ptr = _lib.dada2_table_bimera(
        mat.ctypes.data_as(ct.POINTER(ct.c_int)),
        nrow, ncol,
        seq_arr,
        ct.c_double(min_fold),
        min_abund,
        int(allow_one_off),
        min_one_off_par_dist,
        match, mismatch, gap_p, max_shift
    )

    if not res_ptr:
        raise RuntimeError("dada2_table_bimera returned NULL")

    res = res_ptr.contents
    result = {
        "nflag": np.ctypeslib.as_array(res.nflag, shape=(ncol,)).copy(),
        "nsam": np.ctypeslib.as_array(res.nsam, shape=(ncol,)).copy(),
    }

    _lib.dada2_chimera_result_free(res_ptr)
    return result
