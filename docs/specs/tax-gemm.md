# Spec: GEMM-based taxonomy scoring (`tax-gemm` branch)

## Status
Speculative — breaks bit-parity with R. Not merged into `parity-perf`.

## Motivation
`assign_taxonomy` spends nearly all its query time in `get_best_genus`:
a sequential float32 sum of `lgk_probability[g, kmer]` over the query's
kmer array, for every genus, repeated 1 (best hit) + 100 (bootstraps)
times per query. The per-element early-break makes it cache-hostile and
SIMD-hostile; the access pattern strides across a 65536-wide row per
genus.

## Idea
Reformulate scoring as dense linear algebra:

1. Per query, gather the used columns once:
   `S = lgk_probability[:, uniq_kmers]` → shape `(ngenus, n_uniq)`.
2. Best hit: `scores = S @ counts` (counts = multiplicity of each unique
   kmer) → argmax.
3. Bootstraps: build the 100 subsample count vectors as a matrix
   `C (n_uniq × 100)` → one SGEMM `S @ C` → column argmaxes.

With SILVA-scale tables (≈7k genera × ≈430 kmers × 100 boots) each query
becomes one ~(7k×430)·(430×100) SGEMM — OpenBLAS territory, expected
10–50× faster than the scalar loop, and it composes with the existing
OMP-over-queries parallelism.

## Why it breaks parity
- Summation order and FMA usage inside GEMM differ from the sequential
  scalar loop → last-ulp differences in float32 scores.
- Argmax on near-tied genera can flip vs R. Frequency: comparable to R's
  own run-to-run nondeterminism (upstream breaks exact ties with an
  unseeded RNG), but systematic rather than random.
- The early-break optimisation is abandoned (GEMM computes everything),
  so worst-case FLOPs rise; wall-clock still wins on throughput.

## Parity containment
- Keep as `method="gemm"` opt-in parameter; default stays the exact
  scalar path.
- Acceptance test: on SILVA + real ASVs, genus assignments must agree
  with the exact path for all queries whose top-2 score gap exceeds an
  epsilon; bootstrap values within ±2 of the exact path in ≥99.9% of
  cells (mirrors R's own noise).

## Sketch
- C side: per-query gather into a thread-local `S` buffer; call
  `cblas_sgemv` / `cblas_sgemm` (link OpenBLAS, optional at build time —
  fall back to scalar when unavailable).
- Bootstrap count matrices reuse the R-RNG unifs stream so the *samples*
  are identical to the exact path; only the arithmetic differs.
