"""Unit tests for the model-based worker/thread scheduler (issue #5)."""

import numpy as np

from papa2 import _schedule as sch


def test_lpt_order_heaviest_first_and_stable():
    w = np.array([3.0, 10.0, 1.0, 10.0])
    order = list(sch.lpt_order(w))
    assert order == [1, 3, 0, 2]  # ties keep input order


def test_lpt_makespan_balances():
    times = np.array([4.0, 3.0, 3.0, 2.0])
    assert sch.lpt_makespan(times, 2) == 6.0
    assert sch.lpt_makespan(times, 1) == 12.0
    assert sch.lpt_makespan(times, 4) == 4.0


def test_choose_split_env_overrides_win():
    w = np.ones(24)
    assert sch.choose_split(w, 32, n_workers_env=24) == (24, 1)
    assert sch.choose_split(w, 32, omp_threads_env=8) == (4, 8)
    assert sch.choose_split(w, 32, n_workers_env=6, omp_threads_env=5) == (6, 5)


def test_choose_split_single_sample_gets_all_cores():
    w, t = sch.choose_split([100.0], 32)
    assert (w, t) == (1, 32)


def test_choose_split_few_heavy_samples_prefers_narrow_pool():
    """24 equal heavy samples on 32 cores: the sweep measured 4x8 and 8x4
    fastest, 24x1 slowest — the model must not pick the wide pool."""
    w, t = sch.choose_split(np.ones(24) * 1e6, 32)
    assert w <= 8
    assert w * t <= 32
    assert t >= 4


def test_choose_split_never_exceeds_samples_or_cores():
    for n in (1, 2, 5, 24, 120, 500):
        for cores in (1, 4, 32):
            w, t = sch.choose_split(np.ones(n), cores)
            assert 1 <= w <= max(1, n)
            assert 1 <= t
            assert w * t <= max(cores, w)  # threads fill remaining width


def test_choose_split_skewed_weights_narrows_pool():
    """One giant sample among small ones: makespan is bound by the giant,
    so it should get width rather than a 1-thread slot."""
    weights = np.array([1e7] + [1e4] * 31)
    w, t = sch.choose_split(weights, 32)
    assert t > 1


def test_measured_weights_feedback_and_fallback():
    fallback = np.array([5.0, 5.0, 5.0])
    results = [{"_walltime": 1.0}, {}, {"_walltime": 9.0}]
    out = sch.measured_weights(results, fallback)
    assert list(out) == [1.0, 5.0, 9.0]


def test_file_weights_missing_file():
    w = sch.file_weights(["/nonexistent/x.fastq.gz"])
    assert w[0] == 1.0
