"""Model-based scheduling of sample-level workers vs. within-sample threads.

DRAFT (issue #5) — the split of cores between worker processes and OpenMP
threads is chosen from the shape of the data instead of a fixed ratio,
and samples are submitted largest-first (LPT) so a heavy sample cannot
straggle at the end of the pool.  Scheduling only: arithmetic, and
therefore byte-parity with R, is untouched; results are always returned
in input order.

Model
-----
Each sample is a malleable task.  Its weight is estimated from the derep
object (``W_i ~ U_i^ALPHA * L_i``) or, for file inputs, from the gzipped
file size.  Runtime on ``t`` threads follows an Amdahl curve with an
efficiency decay::

    T_i(t) = W_i * ((1 - P_PARALLEL) + P_PARALLEL / (EFFICIENCY[t] * t))

For each candidate worker count the predicted makespan is an LPT packing
of ``{T_i(t)}`` onto the workers (a 4/3-approximation), and the width
with the smallest predicted makespan wins.

The constants below are deliberately conservative and machine-agnostic;
the ``benchmarks/real_data`` harness logs the data needed to calibrate
them if that ever becomes worthwhile.  ``DADA2_WORKERS`` /
``DADA2_OMP_THREADS`` override everything.
"""

from __future__ import annotations

import heapq
import os
from typing import List, Optional, Sequence, Tuple

import numpy as np

# Cost-model constants (issue #5).  ALPHA: superlinear growth of per-sample
# cost with unique count (clusters grow with uniques).  P_PARALLEL: fraction
# of per-sample work in the OpenMP comparison loop.  EFFICIENCY: parallel
# efficiency of that loop at a given thread count (memory bandwidth,
# heterogeneous cores).
ALPHA = 1.4
P_PARALLEL = 0.92
EFFICIENCY = {1: 1.00, 2: 0.97, 4: 0.92, 8: 0.85, 16: 0.75, 32: 0.62}


def _efficiency(t: int) -> float:
    """Interpolated parallel efficiency at t threads."""
    if t in EFFICIENCY:
        return EFFICIENCY[t]
    keys = sorted(EFFICIENCY)
    if t <= keys[0]:
        return EFFICIENCY[keys[0]]
    if t >= keys[-1]:
        return EFFICIENCY[keys[-1]]
    lo = max(k for k in keys if k <= t)
    hi = min(k for k in keys if k >= t)
    f = (t - lo) / (hi - lo)
    return EFFICIENCY[lo] * (1 - f) + EFFICIENCY[hi] * f


# ---------------------------------------------------------------------------
# Weights
# ---------------------------------------------------------------------------

def derep_weights(dereps: Sequence[dict]) -> np.ndarray:
    """Per-sample cost estimates from derep objects: U^ALPHA * L."""
    w = np.empty(len(dereps), dtype=np.float64)
    for i, d in enumerate(dereps):
        u = max(1, len(d["seqs"]))
        length = len(d["seqs"][0]) if d["seqs"] else 1
        w[i] = (u ** ALPHA) * length
    return w


def file_weights(paths: Sequence[str]) -> np.ndarray:
    """Per-sample cost proxies for file inputs: gzipped size in bytes.

    Uniques are unknown before the worker dereplicates, and compressed
    size tracks read count x length well enough for load balancing.
    """
    w = np.empty(len(paths), dtype=np.float64)
    for i, p in enumerate(paths):
        try:
            w[i] = max(1, os.path.getsize(p))
        except OSError:
            w[i] = 1.0
    return w


def measured_weights(results: Sequence[dict],
                     fallback: np.ndarray) -> np.ndarray:
    """Weights from a previous pass's measured per-sample wall times.

    Used between selfConsist rounds and for pseudo-pooling's second pass
    (issue #5): measured feedback beats any a-priori model.  Results
    lacking a measurement fall back to the model weight.
    """
    w = fallback.astype(np.float64).copy()
    for i, res in enumerate(results):
        t = res.get("_walltime") if isinstance(res, dict) else None
        if t is not None and t > 0:
            w[i] = t
    return w


# ---------------------------------------------------------------------------
# Makespan model
# ---------------------------------------------------------------------------

def predict_time(weight: float, t: int) -> float:
    """Predicted runtime of one sample on t threads (Amdahl + efficiency)."""
    return weight * ((1.0 - P_PARALLEL) + P_PARALLEL / (_efficiency(t) * t))


def lpt_makespan(times: np.ndarray, n_workers: int) -> float:
    """Makespan of an LPT (longest-processing-time-first) packing."""
    loads = [0.0] * n_workers
    heapq.heapify(loads)
    for t in sorted(times, reverse=True):
        heapq.heappush(loads, heapq.heappop(loads) + float(t))
    return max(loads)


def lpt_order(weights: np.ndarray) -> np.ndarray:
    """Submission order: heaviest sample first (stable for ties)."""
    return np.argsort(-np.asarray(weights, dtype=np.float64), kind="stable")


def choose_split(weights: Sequence[float], cores: int,
                 n_workers_env: int = 0,
                 omp_threads_env: int = 0) -> Tuple[int, int]:
    """Choose (n_workers, omp_threads) for the given sample weights.

    Environment-style overrides win when nonzero (matching the semantics
    of DADA2_WORKERS / DADA2_OMP_THREADS).  Otherwise the candidate
    worker counts {1, C/8, C/4, C/2, C} (bounded by the sample count)
    are compared by predicted LPT makespan.
    """
    weights = np.asarray(list(weights), dtype=np.float64)
    n = len(weights)
    cores = max(1, cores)
    if n == 0:
        return 1, cores

    if n_workers_env > 0:
        w = min(n_workers_env, max(1, n))
        t = omp_threads_env if omp_threads_env > 0 else max(1, cores // w)
        return w, t
    if omp_threads_env > 0:
        w = min(n, max(1, cores // omp_threads_env))
        return w, omp_threads_env

    # Candidate widths, clamped to the sample count (never dropped: with
    # 2 samples the right answer may be exactly 2 workers).
    candidates = sorted({min(max(1, n), c) for c in
                         (1, max(1, cores // 8), max(1, cores // 4),
                          max(1, cores // 2), cores)})

    best = None
    for w in candidates:
        t = max(1, cores // w)
        times = np.array([predict_time(x, t) for x in weights])
        span = lpt_makespan(times, w)
        # Deterministic tie-break: prefer fewer workers (lower memory)
        key = (span, w)
        if best is None or key < best[0]:
            best = (key, w, t)
    return best[1], best[2]
