# Spec: GPU denoising path (`gpu` branch)

## Status
Speculative — parity-risky, and this box currently has no working NVIDIA
driver (nvidia-smi cannot reach it), so nothing here is testable locally
today. Recorded for when hardware is available.

## Candidate kernels (by measured CPU cost)
1. `b_compare_omp` inner work: banded vectorized NW alignment + lambda
   computation of every raw against a new cluster center. Embarrassingly
   parallel over raws; the alignment DP is the GPU-friendly part.
   This dominates pooled DENOISE runs (the 1298-sample microscape run
   spent 1.7 h in one DENOISE task).
2. Taxonomy scoring (see tax-gemm — on GPU it is one cuBLAS SGEMM per
   query batch).
3. `final_subs` alignment pass (same kernel as 1).

## Parity concerns
- lambda uses double-precision products of error-matrix entries along
  alignments: GPU FMA/parallel-reduction order changes last-ulp values,
  and OMEGA_A comparisons sit on hard thresholds (1e-40), so cluster
  budding decisions can flip on borderline p-values.
- Mitigation: compute lambda in double on GPU with the same sequential
  order per alignment (feasible: per-alignment loop is short), keep
  reductions per-raw sequential. If that holds, results may remain
  bit-identical — needs empirical verification against the CPU path on
  the 6- and 120-sample benchmarks before any default flip.

## Prior art
rec3141/dada2@gpu-python (the papa2 ancestor) carries earlier GPU
experiments; mine for kernels, but note it predates the -O2 CPU baseline
(the CPU got ~7x faster in papa2, so GPU wins must be re-measured).

## Plan
1. CUDA (or HIP) port of nwalign_vectorized + compute_lambda over raw
   batches, gated behind DADA2_GPU=1 env + runtime device probe.
2. Parity harness: byte-compare denoised outputs vs CPU on both
   benchmarks; if bit-identical, it can graduate to a default;
   otherwise stays opt-in.
