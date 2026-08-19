"""Main DADA2 denoising pipeline."""

import sys
import os
import numpy as np

# Worker pools must not fork a process that has already run an OpenMP
# region (the child inherits a broken libgomp runtime and can deadlock).
# loky provides spawn-safe workers without requiring callers to guard
# their scripts with __main__; plain spawn is the fallback.
try:
    from loky import get_reusable_executor as _loky_executor
except ImportError:
    _loky_executor = None
import multiprocessing
from concurrent.futures import ProcessPoolExecutor
_MP_CTX = multiprocessing.get_context("spawn")


def _make_pool(n_workers):
    if _loky_executor is not None:
        return _loky_executor(max_workers=n_workers)
    return ProcessPoolExecutor(max_workers=n_workers, mp_context=_MP_CTX)
from . import _cdada
from .io import combine_dereps, derep_fastq
from .error import loess_errfun, get_initial_err

DADA_OPTS = {
    "OMEGA_A": 1e-40,
    "OMEGA_P": 1e-4,
    "OMEGA_C": 1e-40,
    "DETECT_SINGLETONS": False,
    "USE_KMERS": True,
    "KDIST_CUTOFF": 0.42,
    "MAX_CONSIST": 10,
    "MATCH": 5,
    "MISMATCH": -4,
    "GAP_PENALTY": -8,
    "BAND_SIZE": 16,
    "VECTORIZED_ALIGNMENT": True,
    "MAX_CLUST": 0,
    "MIN_FOLD": 1,
    "MIN_HAMMING": 1,
    "MIN_ABUNDANCE": 1,
    "USE_QUALS": True,
    "HOMOPOLYMER_GAP_PENALTY": None,
    "SSE": 2,
    "GAPLESS": True,
    "GREEDY": True,
    "PSEUDO_ABUNDANCE": float("inf"),
    "PSEUDO_PREVALENCE": 2,
}


def _run_one_sample(args):
    """Worker function for parallel dada() calls.

    Top-level function so it works with both ThreadPoolExecutor and
    ProcessPoolExecutor. For ProcessPoolExecutor, re-imports _cdada
    in each subprocess since ctypes handles can't be pickled.
    """
    drp, erri, opts, max_clust, verbose, nthreads, prior_seqs = args
    try:
        _cd = _cdada  # ThreadPoolExecutor: module already imported
    except NameError:
        from papa2 import _cdada as _cd  # ProcessPoolExecutor: reimport

    seqs = drp["seqs"]
    if len(seqs) == 0:
        return {"cluster_seqs": [], "cluster_abunds": np.array([]),
                "trans": np.zeros((16, erri.shape[1]), dtype=np.int32),
                "map": np.array([]), "pval": np.array([])}

    homo_gap = opts["GAP_PENALTY"] if opts["HOMOPOLYMER_GAP_PENALTY"] is None else opts["HOMOPOLYMER_GAP_PENALTY"]

    priors = None
    if prior_seqs:
        priors = np.fromiter((s in prior_seqs for s in seqs), dtype=np.int32,
                             count=len(seqs))

    res = _cd.run_dada(
        seqs, drp["abundances"], erri, drp["quals"], priors=priors,
        match=opts["MATCH"], mismatch=opts["MISMATCH"], gap_pen=opts["GAP_PENALTY"],
        use_kmers=opts["USE_KMERS"], kdist_cutoff=opts["KDIST_CUTOFF"],
        band_size=opts["BAND_SIZE"],
        omega_a=opts["OMEGA_A"], omega_p=opts["OMEGA_P"], omega_c=opts["OMEGA_C"],
        detect_singletons=opts["DETECT_SINGLETONS"], max_clust=max_clust,
        min_fold=opts["MIN_FOLD"], min_hamming=opts["MIN_HAMMING"],
        min_abund=opts["MIN_ABUNDANCE"],
        use_quals=opts["USE_QUALS"], vectorized_alignment=opts["VECTORIZED_ALIGNMENT"],
        homo_gap_pen=homo_gap,
        multithread=nthreads, verbose=verbose,
        sse=opts["SSE"], gapless=opts["GAPLESS"], greedy=opts["GREEDY"],
    )

    res["denoised"] = {seq: ab for seq, ab in
                       zip(res["cluster_seqs"], res["cluster_abunds"])}
    return res


def _run_one_file_sample(args):
    """Dereplicate one FASTQ and denoise it in the same worker."""
    filepath, err, opts, max_clust, verbose, nthreads, prior_seqs = args
    drp = derep_fastq(filepath, verbose=verbose)

    max_q = 0
    if drp["quals"].size > 0 and not np.all(np.isnan(drp["quals"])):
        max_q = int(np.nanmax(drp["quals"])) + 1

    erri = err.copy()
    if max_q > erri.shape[1]:
        extra = np.tile(erri[:, -1:], (1, max_q - erri.shape[1]))
        erri = np.hstack([erri, extra])

    return _run_one_sample(
        (drp, erri, opts, max_clust, verbose, nthreads, prior_seqs))


def set_dada_opt(**kwargs):
    for k, v in kwargs.items():
        if k in DADA_OPTS:
            DADA_OPTS[k] = v
        else:
            raise KeyError(f"Unknown DADA option: {k}")


def get_dada_opt(key=None):
    if key is None:
        return dict(DADA_OPTS)
    return DADA_OPTS[key]


def _pseudo_priors_from(results, opts):
    """Sequences qualifying as pseudo-pooling priors (R's rule): present in
    at least PSEUDO_PREVALENCE samples, or with total abundance of at
    least PSEUDO_ABUNDANCE across samples."""
    prevalence = {}
    total = {}
    for res in results:
        for sq, ab in res["denoised"].items():
            if ab > 0:
                prevalence[sq] = prevalence.get(sq, 0) + 1
                total[sq] = total.get(sq, 0) + ab
    return frozenset(
        sq for sq in prevalence
        if prevalence[sq] >= opts["PSEUDO_PREVALENCE"]
        or total[sq] >= opts["PSEUDO_ABUNDANCE"])


def _expand_pooled(pooled, pooled_derep, derep_in):
    """Split a pooled dada result into per-sample results (R's pool=TRUE
    expansion): keep the clusters each sample's uniques map to, renumber
    them, remap the per-unique map, and recompute per-sample abundances.
    trans and pval keep their pooled values, as in R."""
    pooled_map = np.asarray(pooled["map"], dtype=np.int64)
    nclust = len(pooled["cluster_seqs"])
    seq_index = {sq: i for i, sq in enumerate(pooled_derep["seqs"])}

    out = []
    for drp in derep_in:
        uidx = np.fromiter((seq_index[sq] for sq in drp["seqs"]),
                           dtype=np.int64, count=len(drp["seqs"]))
        m = pooled_map[uidx]
        keep = np.zeros(nclust, dtype=bool)
        keep[m[m >= 0]] = True
        new_bi = np.cumsum(keep) - 1  # pooled cluster idx -> sample idx

        map_i = np.where(m >= 0, new_bi[np.maximum(m, 0)], -1).astype(np.int32)
        n_keep = int(keep.sum())
        abunds = np.zeros(n_keep, dtype=np.int64)
        valid = map_i >= 0
        np.add.at(abunds, map_i[valid],
                  np.asarray(drp["abundances"], dtype=np.int64)[valid])

        keep_idx = np.flatnonzero(keep)
        res = {
            "cluster_seqs": [pooled["cluster_seqs"][j] for j in keep_idx],
            "cluster_abunds": abunds.astype(np.int32),
            "map": map_i,
            # Pooled values, pruned per cluster where R prunes them:
            "cluster_n0": np.asarray(pooled["cluster_n0"])[keep_idx]
                if "cluster_n0" in pooled else None,
            "cluster_n1": np.asarray(pooled["cluster_n1"])[keep_idx]
                if "cluster_n1" in pooled else None,
            "cluster_nunq": np.asarray(pooled["cluster_nunq"])[keep_idx]
                if "cluster_nunq" in pooled else None,
            # Pooled (not per-sample) values, as in R:
            "trans": pooled["trans"],
            "pval": pooled["pval"],
        }
        res["denoised"] = {sq: int(ab) for sq, ab in
                           zip(res["cluster_seqs"], res["cluster_abunds"])}
        out.append(res)
    return out


def dada(derep, err=None, error_estimation_function=None, self_consist=False,
         pool=False, priors=None, verbose=True, **opts):
    """Run DADA2 denoising on one or more dereplicated samples.

    Args:
        derep: dict from derep_fastq(), or list of dicts, or FASTQ filepath(s)
        err: numpy array (16, ncol) error matrix, or None for self-consistent learning
        error_estimation_function: callable(trans) -> err_matrix, default loess_errfun
        self_consist: bool, iterate until error model converges
        pool: False (default) processes samples independently; True pools
            all samples into one inference (R's pool=TRUE, including the
            per-sample expansion of the pooled result); "pseudo" runs
            R's pseudo-pooling (a second pass with the first pass's
            consistently-observed sequences as priors, controlled by
            PSEUDO_PREVALENCE / PSEUDO_ABUNDANCE).
        priors: sequences with prior evidence of existence (R's priors);
            they are evaluated against OMEGA_P instead of OMEGA_A
        verbose: bool

    Returns:
        dict (single sample) or list of dicts, each with:
            denoised: dict {seq: abundance}
            cluster_seqs, cluster_abunds, trans, map, pval, err_in, err_out

    Environment variables:
        DADA2_CORES: total cores to plan for (default: os.cpu_count();
            set this under containers/schedulers that allocate fewer)
        DADA2_WORKERS: number of sample-level worker processes
            (0 = auto-detect, default)
        DADA2_OMP_THREADS: OpenMP threads per worker for the within-sample
            comparison loop (0 = auto: cores split evenly across workers)
    """
    if error_estimation_function is None:
        error_estimation_function = loess_errfun

    file_inputs = None

    # Normalize input to list of derep dicts
    if isinstance(derep, str):
        file_inputs = [derep]
        derep = [derep_fastq(derep, verbose=verbose)]
    elif isinstance(derep, dict):
        derep = [derep]
    elif isinstance(derep, list) and len(derep) > 0 and isinstance(derep[0], str):
        file_inputs = list(derep)
        need_dereps = self_consist or err is None or pool is True
        derep = [derep_fastq(f, verbose=verbose) for f in derep] if need_dereps else derep

    single = len(file_inputs) == 1 if file_inputs is not None else len(derep) == 1

    # Merge options
    o = dict(DADA_OPTS)
    o.update(opts)

    # Priors and pooling (R: pool has no effect on a single sample)
    prior_seqs = frozenset(priors) if priors else frozenset()
    pseudo = False
    pseudo_priors = frozenset()
    derep_in = None
    n_inputs = len(file_inputs) if file_inputs is not None else len(derep)
    if n_inputs <= 1:
        pool = False
    if pool is True:
        derep_in = derep
        file_inputs = None  # pooled inference runs on the combined derep
        derep = [combine_dereps(derep_in)]
        if verbose:
            print(f"{len(derep_in)} samples were pooled: "
                  f"{int(derep[0]['abundances'].sum())} reads in "
                  f"{len(derep[0]['seqs'])} unique sequences.")
    elif pool == "pseudo":
        pool = False
        pseudo = True
    elif pool is not False:
        raise ValueError("Invalid pool argument.")

    # Initialize error matrix (matching R: all 1.0)
    initialize_err = False
    if self_consist and err is None:
        max_q = 0
        for d in derep:
            if d["quals"].size > 0 and not np.all(np.isnan(d["quals"])):
                max_q = max(max_q, int(np.nanmax(d["quals"])) + 1)
        err = get_initial_err(max(41, max_q))
        initialize_err = True

    if err is None:
        raise ValueError("Error matrix (err) must be provided unless self_consist=True")

    err_history = []
    nconsist = 0 if initialize_err else 1  # R: init at 0, otherwise start at 1

    # CPU-only standalone work must use processes, not threads:
    # the shared library is not thread-safe across concurrent run_dada
    # invocations.  Cores are split between sample-level workers and the
    # within-sample OpenMP comparison loop, so a lone pooled sample uses
    # every core while many small samples get one core each.
    # DADA2_CORES bounds the total-core estimate (containers and batch
    # schedulers often allocate fewer cores than os.cpu_count reports)
    cores = int(os.environ.get("DADA2_CORES", "0")) or (os.cpu_count() or 1)
    n_workers = int(os.environ.get("DADA2_WORKERS", "0"))
    if n_workers == 0:
        n_workers = min(len(derep), cores)
    omp_threads = int(os.environ.get("DADA2_OMP_THREADS", "0"))
    if omp_threads == 0:
        omp_threads = max(1, cores // max(1, min(n_workers, len(derep))))
    use_parallel = len(derep) > 1 and n_workers > 1

    if file_inputs is not None and not self_consist:
        max_clust_iter = o["MAX_CLUST"]

        def _run_files(active_priors):
            work_args = [(fpath, err, o, max_clust_iter, verbose, omp_threads,
                          active_priors) for fpath in file_inputs]
            if use_parallel:
                with _make_pool(n_workers) as pool_:
                    return list(pool_.map(_run_one_file_sample, work_args))
            return [_run_one_file_sample(a) for a in work_args]

        results = _run_files(prior_seqs)
        if pseudo:
            # R's pseudo-pooling: a second pass with the first pass's
            # consistently-observed sequences added to the priors.
            pseudo_priors = _pseudo_priors_from(results, o)
            results = _run_files(prior_seqs | pseudo_priors)

        for res in results:
            res["err_in"] = err_history[0] if err_history else None
            res["err_out"] = err

        if single:
            return results[0]
        return results

    while True:
        if verbose and self_consist:
            sys.stdout.write(f"   selfConsist step {nconsist}")
            sys.stdout.flush()

        if nconsist > 0:
            err_history.append(err.copy())

        # R uses MAX_CLUST=1 on the initialization pass (nconsist==1 after init)
        max_clust_iter = 1 if initialize_err else o["MAX_CLUST"]

        # Extend error matrix once for all samples
        max_q_all = max((int(np.nanmax(d["quals"])) + 1 if d["quals"].size and not np.all(np.isnan(d["quals"])) else 0) for d in derep)
        erri = err.copy()
        if max_q_all > erri.shape[1]:
            extra = np.tile(erri[:, -1:], (1, max_q_all - erri.shape[1]))
            erri = np.hstack([erri, extra])

        active_priors = prior_seqs | pseudo_priors
        if use_parallel:
            # Parallel multi-sample execution.
            work_args = [(drp, erri, o, max_clust_iter, verbose, omp_threads,
                          active_priors) for drp in derep]
            with _make_pool(n_workers) as pool_:
                results = list(pool_.map(_run_one_sample, work_args))
            if verbose and self_consist:
                sys.stdout.write("." * len(derep))
                sys.stdout.flush()
        else:
            # Sequential path for single-sample or single-worker execution.
            results = []
            for drp in derep:
                if verbose and self_consist:
                    sys.stdout.write(".")
                    sys.stdout.flush()
                res = _run_one_sample(
                    (drp, erri, o, max_clust_iter, verbose, omp_threads,
                     active_priors))
                results.append(res)

        trans_list = [r["trans"] for r in results]

        if verbose and self_consist:
            print()

        # Accumulate transitions
        cur_trans = _accumulate_trans(trans_list)

        # Estimate new error rates
        err = error_estimation_function(cur_trans)

        # After initialization pass: set self-transitions to 1.0 (matching R)
        if initialize_err:
            err[0, :] = 1.0   # A2A
            err[5, :] = 1.0   # C2C
            err[10, :] = 1.0  # G2G
            err[15, :] = 1.0  # T2T
            initialize_err = False

        # Termination: R breaks when done UNLESS a pseudo-pooling second
        # pass is still owed (pseudo requires reaching nconsist >= 2).
        converged = any(np.array_equal(err, h) for h in err_history)
        done = (not self_consist) or converged or nconsist >= o["MAX_CONSIST"]
        if done and (not pseudo or nconsist >= 2):
            if self_consist and verbose:
                if converged:
                    print(f"Convergence after {nconsist} rounds.")
                elif nconsist >= o["MAX_CONSIST"]:
                    print("Self-consistency loop terminated before convergence.")
            break

        # Pseudo-pooling priors for the next pass (R computes these at the
        # end of every full loop iteration)
        if pseudo and nconsist >= 1:
            pseudo_priors = _pseudo_priors_from(results, o)

        nconsist += 1

    # Expand the pooled result into per-sample results (R's pool=TRUE)
    if derep_in is not None:
        results = _expand_pooled(results[0], derep[0], derep_in)
        single = len(derep_in) == 1

    # Attach error info to results
    for res in results:
        res["err_in"] = list(err_history) if self_consist else err_history[0]
        res["err_out"] = err

    if single:
        return results[0]
    return results


def learn_errors(fastq_files, nbases=1e8, error_estimation_function=None,
                 verbose=True, **opts):
    """Learn error rates from FASTQ files.

    Args:
        fastq_files: list of FASTQ file paths
        nbases: target number of bases to use for learning
        error_estimation_function: callable, default loess_errfun
        verbose: bool

    Returns:
        numpy array (16, ncol) of learned error rates
    """
    if isinstance(fastq_files, str):
        fastq_files = [fastq_files]

    # Accumulate samples (in order) until we have enough bases.  Matches R:
    # bases are counted as sum(abundance * sequence length) and accumulation
    # stops after the sample that pushes the total strictly past nbases.
    # Dereplication runs in a worker pool; results are consumed in file
    # order so the selected sample set is identical to a sequential scan.
    dereps = []
    total_bases = 0

    def _n_bases(drp):
        if not drp["seqs"]:
            return 0
        lens = np.fromiter((len(s) for s in drp["seqs"]), dtype=np.int64,
                           count=len(drp["seqs"]))
        return int((np.asarray(drp["abundances"], dtype=np.int64) * lens).sum())

    n_workers = (int(os.environ.get("DADA2_WORKERS", "0"))
                 or int(os.environ.get("DADA2_CORES", "0"))
                 or (os.cpu_count() or 1))
    n_workers = min(n_workers, len(fastq_files))
    if n_workers > 1 and len(fastq_files) > 1:
        with _make_pool(n_workers) as pool:
            for drp in pool.map(derep_fastq, fastq_files):
                if verbose:
                    print(f"Read {int(drp['abundances'].sum())} reads, "
                          f"{len(drp['seqs'])} unique sequences")
                dereps.append(drp)
                total_bases += _n_bases(drp)
                if total_bases > nbases:
                    try:
                        pool.shutdown(wait=False, cancel_futures=True)
                    except TypeError:  # loky's shutdown has no cancel_futures
                        pool.shutdown(wait=False, kill_workers=True)
                    break
    else:
        for fl in fastq_files:
            drp = derep_fastq(fl, verbose=verbose)
            dereps.append(drp)
            total_bases += _n_bases(drp)
            if total_bases > nbases:
                break

    if verbose:
        n_reads_total = sum(d["abundances"].sum() for d in dereps)
        print(f"{int(total_bases)} total bases in {n_reads_total} reads "
              f"from {len(dereps)} samples will be used for learning the error rates.")
        print("Initializing error rates to maximum possible estimate.")

    # R's learnErrors wrapper forces OMEGA_C=0 during learning without
    # changing the general dada() default.
    opts = dict(opts)
    opts.setdefault("OMEGA_C", 0)

    # Run dada with self-consistency
    results = dada(dereps, err=None, error_estimation_function=error_estimation_function,
                   self_consist=True, verbose=verbose, **opts)

    if isinstance(results, dict):
        results = [results]

    return results[0]["err_out"]


def _accumulate_trans(trans_list):
    """Sum transition matrices, zero-padding to max column count."""
    if not trans_list:
        return np.zeros((16, 0), dtype=np.int32)

    max_ncol = max(t.shape[1] for t in trans_list if t.size > 0)
    if max_ncol == 0:
        return np.zeros((16, 0), dtype=np.int32)

    acc = np.zeros((16, max_ncol), dtype=np.int64)
    for t in trans_list:
        if t.size == 0:
            continue
        nc = t.shape[1]
        acc[:, :nc] += t

    return acc.astype(np.int32)
