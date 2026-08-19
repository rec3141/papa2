"""Regression tests for R-parity behaviours that don't need R installed.

Golden values in this file were generated with R 4.6 / dada2 1.40 and are
frozen here; if these tests fail, byte-parity with R is broken.
"""

import ctypes as ct
import gzip
import os

import numpy as np
import pytest

import papa2
from papa2 import _cdada
from papa2.utils import is_phix, match_ref, _phix_genome, _rc


# ---------------------------------------------------------------------------
# R RNG stream
# ---------------------------------------------------------------------------

def test_r_rng_matches_r_runif():
    """set.seed(100); runif(6) in R, reproduced bit-for-bit."""
    expected = [
        0.307766110869124532, 0.257672501029446721, 0.552322433330118656,
        0.056383150396868587, 0.468549283919855952, 0.483770735096186399,
    ]
    out = np.zeros(6, dtype=np.float64)
    _cdada._lib.r_rng_runif_fill.restype = None
    _cdada._lib.r_rng_runif_fill.argtypes = [
        ct.c_uint32, ct.POINTER(ct.c_double), ct.c_longlong]
    _cdada._lib.r_rng_runif_fill(
        ct.c_uint32(100), out.ctypes.data_as(ct.POINTER(ct.c_double)),
        ct.c_longlong(6))
    assert out.tolist() == expected


# ---------------------------------------------------------------------------
# PhiX / matchRef (R's C_matchRef semantics)
# ---------------------------------------------------------------------------

def test_is_phix_flags_genome_fragments():
    g = _phix_genome()
    assert len(g) == 5386
    frag = g[1000:1250]
    rc_frag = _rc(g[3000:3250])
    random_seq = ("ACGT" * 80)[:250]
    flags = is_phix([frag, rc_frag, random_seq])
    assert flags.tolist() == [True, True, False]


def test_match_ref_circular_and_nonoverlapping():
    # Kmers spanning the circular wrap point must hit.
    ref = "ACGTACGTAAGGTTCCAACCGGTTACGATCGA"  # 32 bp
    wrap = ref[-8:] + ref[:8]
    assert match_ref([wrap], ref, word_size=16)[0] >= 1
    # Non-overlapping counting: R skips word_size positions after each hit.
    q = ref[:30]
    n_over = match_ref([q], ref, word_size=16, non_overlapping=False)[0]
    n_non = match_ref([q], ref, word_size=16, non_overlapping=True)[0]
    assert n_over == 15  # every window of a perfect substring hits
    assert n_non == 1


# ---------------------------------------------------------------------------
# Filter semantics (order of operations frozen from R's fastqFilter)
# ---------------------------------------------------------------------------

def _write_fastq(path, records):
    with gzip.open(path, "wt") as f:
        for name, seq, qual in records:
            f.write(f"@{name}\n{seq}\n+\n{qual}\n")


def _read_fastq(path):
    out = []
    with gzip.open(path, "rt") as f:
        while True:
            h = f.readline().strip()
            if not h:
                break
            s = f.readline().strip()
            f.readline()
            q = f.readline().strip()
            out.append((h, s, q))
    return out


def test_trunclen_counts_from_original_five_prime(tmp_path):
    """R: with trimLeft=10, truncLen=40, output reads are 30 bases."""
    fin = str(tmp_path / "in.fastq.gz")
    fout = str(tmp_path / "out.fastq.gz")
    seq = "ACGT" * 15  # 60 bp
    _write_fastq(fin, [("r1", seq, "I" * 60)])
    n_in, n_out = papa2.fastq_filter(
        fin, fout, trim_left=10, trunc_len=40, trunc_q=None,
        rm_phix=False, min_len=1)
    assert (n_in, n_out) == (1, 1)
    recs = _read_fastq(fout)
    assert len(recs[0][1]) == 30
    assert recs[0][1] == seq[10:40]


def test_truncq_applies_after_trim_left(tmp_path):
    """A low-quality base inside the trimmed-off prefix must not truncate
    the read (R trims left before trimTails)."""
    fin = str(tmp_path / "in.fastq.gz")
    fout = str(tmp_path / "out.fastq.gz")
    seq = "ACGT" * 15
    qual = "#" + "I" * 59  # Q2 at position 0
    _write_fastq(fin, [("r1", seq, qual)])
    n_in, n_out = papa2.fastq_filter(
        fin, fout, trim_left=5, trunc_q=2, rm_phix=False, min_len=20)
    assert (n_in, n_out) == (1, 1)
    assert len(_read_fastq(fout)[0][1]) == 55


def test_min_q_strictly_greater_and_gated_on_truncq(tmp_path):
    """R keeps reads with min(qual) > minQ, and only applies the check
    when minQ > truncQ."""
    fin = str(tmp_path / "in.fastq.gz")
    fout = str(tmp_path / "out.fastq.gz")
    # All qualities exactly Q10 ('+' = 43 -> Q10)
    _write_fastq(fin, [("r1", "ACGT" * 15, "+" * 60)])
    # min_q=10: min == minQ is NOT strictly greater -> dropped
    _, n_out = papa2.fastq_filter(fin, fout, min_q=10, trunc_q=2,
                                  rm_phix=False, min_len=20)
    assert n_out == 0
    # min_q=9: min > minQ -> kept
    _, n_out = papa2.fastq_filter(fin, fout, min_q=9, trunc_q=2,
                                  rm_phix=False, min_len=20)
    assert n_out == 1


def test_paired_orient_swaps_reads(tmp_path):
    """R's fastqPairedFilter swaps F/R (no reverse-complement) when the
    reverse read starts with orient.fwd."""
    finF = str(tmp_path / "F.fastq.gz")
    finR = str(tmp_path / "R.fastq.gz")
    outF = str(tmp_path / "oF.fastq.gz")
    outR = str(tmp_path / "oR.fastq.gz")
    primer = "TCAGC"
    seqA = primer + "A" * 55
    seqB = "G" * 60
    # Pair 1: forward starts with primer (kept as-is).
    # Pair 2: reverse starts with primer (swapped).
    _write_fastq(finF, [("p1 1", seqA, "I" * 60), ("p2 1", seqB, "I" * 60)])
    _write_fastq(finR, [("p1 2", seqB, "I" * 60), ("p2 2", seqA, "I" * 60)])
    n_in, n_out = papa2.fastq_paired_filter(
        finF, outF, finR, outR, orient_fwd=primer, trunc_q=None,
        rm_phix=False, min_len=20)
    assert (n_in, n_out) == (2, 2)
    fw = _read_fastq(outF)
    rv = _read_fastq(outR)
    assert [r[1] for r in fw] == [seqA, seqA]
    assert [r[1] for r in rv] == [seqB, seqB]
    # Swapped pair carries the reverse file's header on the forward side
    assert fw[1][0] == "@p2 2"


def test_zero_output_removes_files(tmp_path):
    fin = str(tmp_path / "in.fastq.gz")
    fout = str(tmp_path / "out.fastq.gz")
    _write_fastq(fin, [("r1", "ACGT" * 15, "#" * 60)])  # all Q2
    _, n_out = papa2.fastq_filter(fin, fout, trunc_q=2, rm_phix=False,
                                  min_len=20)
    assert n_out == 0
    assert not os.path.exists(fout)


# ---------------------------------------------------------------------------
# Chimera defaults (R's per-method defaults)
# ---------------------------------------------------------------------------

def test_remove_bimera_denovo_method_defaults():
    seqs = ["ACGTACGTACGTACGTACGT", "TGCATGCATGCATGCATGCA"]
    seqtab = {"table": np.array([[5, 3]], dtype=np.int32), "seqs": seqs}
    # Small tables have no chimeras either way; this exercises default
    # resolution paths without error.
    for method in ("consensus", "pooled", "per-sample"):
        out = papa2.remove_bimera_denovo(seqtab, method=method)
        assert out["table"].shape[1] == 2
