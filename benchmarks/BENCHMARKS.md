# papa2 Benchmarks

This document summarizes the real-data benchmark work done while porting and validating `papa2` against the official R `dada2` package.

It is intentionally a benchmark report, not just a harness README. The harness entrypoints still live in this directory:

- `benchmarks/run_suite.py`
- `benchmarks/run_python_case.py`
- `benchmarks/run_r_case.R`
- `benchmarks/monitor_command.py`
- `benchmarks/plot_results.py`

For the current harness overview, see [README.md](./README.md).

## Baselines

Python baseline:
- `rec3141/papa2`
- local working tree: `/data/papa2`

R baseline:
- untouched upstream `dada2` release `v1.26`
- local checkout: `/tmp/dada2_official`
- commit: `f1e72fb`
- release date: 2022-11-03

This is important: the R comparison was not run against the modified development tree used during the porting work. It was run against a fresh official release.

## Datasets

Two real paired-end datasets were used.

### 1. Curated paired cohort

Source:
- `/matika/seqs/amplicon/cutadapt/515F-806RB/20231105-loutet`

Local copy used for benchmarking:
- `/tmp/papa2_bench_data/20231105-loutet`

Properties:
- 36 complete `R1`/`R2` cutadapt pairs
- 515F-806RB amplicons
- paired-end benchmark on a coherent single run

### 2. Broad 515F pool

Source:
- `/data/scratch/515F-806RB`

Properties:
- large mixed pool of real cutadapt outputs
- user-filtered before benchmarking to remove:
  - 18S samples
  - samples smaller than 100k in size
- used for the wider overnight scaling ladder

## Benchmark design

The goal was to compare real pipeline behavior, not synthetic microbenchmarks.

Pipeline shape:
1. start from cutadapt FASTQs
2. run paired filtering/trimming
3. learn forward/reverse error models
4. denoise forward/reverse reads
5. merge pairs

Important note:
- the cutadapt FASTQs could not be denoised directly in either engine because some reads still contained ambiguous bases (`N`)
- therefore both `papa2` and R `dada2` were benchmarked with filtering included as part of the real workflow

## Parallelism tested

For `papa2`, two parallel regimes were evaluated:

- `core` mode:
  - `DADA2_WORKERS=1`
  - `OMP_NUM_THREADS=N`
- `sample` mode:
  - `DADA2_WORKERS=N`
  - `OMP_NUM_THREADS=1`

For official R `dada2`, the benchmark used numeric `multithread=N`.

The main fair comparison mode turned out to be:
- `papa2 core`
- R `dada2 multithread=N`

`sample` mode was tested, but it never produced a convincing win on the real paired workloads we care about. At higher thread counts it was clearly worse and much heavier on memory.

## Curated paired cohort results

### 4 threads

One paired sample:

| Engine | Wall time | Max RSS |
|---|---:|---:|
| `papa2 core` | `3.84s` | `140 MB` |
| `papa2 sample` | `10.03s` | `138 MB` |
| `dada2` | `34.51s` | `1.11 GB` |

Ten paired samples:

| Engine | Wall time | Max RSS |
|---|---:|---:|
| `papa2 core` | `50.43s` | `484 MB` |
| `papa2 sample` | `50.50s` | `535 MB` |
| `dada2` | `386.75s` | `1.76 GB` |

Interpretation:
- on one large paired sample, `papa2 core` was about `9x` faster than R
- on ten paired samples, `papa2 core` was about `7.7x` faster
- memory use was also much lower
- `sample` mode did not help at 4 threads

### 16 threads

One paired sample:

| Engine | Wall time | Max RSS |
|---|---:|---:|
| `papa2 core` | `2.61s` | `183 MB` |
| `dada2` | `21.24s` | `1.12 GB` |

Ten paired samples:

| Engine | Wall time | Max RSS |
|---|---:|---:|
| `papa2 core` | `30.66s` | `521 MB` |
| `dada2` | `194.3s` | `1.82 GB` |

All 36 paired samples:

| Engine | Wall time | Max RSS |
|---|---:|---:|
| `papa2 core` | `103.08s` | `1.38 GB` |
| `dada2` | `637.0s` | `2.61 GB` |

Interpretation:
- on the coherent 36-sample paired cohort, `papa2 core` was still `6x` to `8x` faster than official R
- memory stayed materially lower as well
- this is the strongest result from the current benchmark set

### `core` vs `sample`

At 16 threads:

One paired sample:
- `core`: `2.90s`
- `sample`: `10.05s`

Ten paired samples:
- `core`: `30.43s`, about `488 MB`
- `sample`: `44.30s`, sampled peak RSS about `8.7 GB`

Conclusion:
- `sample` mode is not the right default for the realistic paired denoising workloads tested here
- `core` mode is the defensible benchmark mode and the one worth keeping

## Broad 515F pool results

The larger overnight ladder was run on `/data/scratch/515F-806RB` after removing 18S samples and very small samples. The goal here was breadth, not strict cohort coherence.

Completed cases at 16 threads:

| Samples | `papa2` wall | `dada2` wall | Speedup | `papa2` RSS | `dada2` RSS |
|---|---:|---:|---:|---:|---:|
| 1 | `2.17s` | `7.25s` | `3.3x` | `96,776 KB` | `860,072 KB` |
| 2 | `2.81s` | `8.29s` | `3.0x` | `97,540 KB` | `923,464 KB` |
| 4 | `37.49s` | `42.97s` | `1.15x` | `207,672 KB` | `1,283,228 KB` |
| 8 | `48.82s` | `58.50s` | `1.20x` | `192,880 KB` | `1,234,348 KB` |
| 16 | `80.91s` | `90.09s` | `1.11x` | `269,836 KB` | `1,355,308 KB` |
| 32 | `195.81s` | `222.15s` | `1.13x` | `433,076 KB` | `1,688,928 KB` |

The overnight run was stopped after this point because the trend was already clear.

Interpretation:
- memory remained clearly better for `papa2`
- the runtime gap compressed dramatically on this mixed pool
- this does not contradict the curated cohort result
- it does indicate that broader real-world throughput is dominated by a different part of the pipeline than the final per-sample `dada(err=...)` pass

## Profiling results

To explain the flattening on the broader benchmark, we profiled a representative `s16` paired case stage by stage.

### Stage breakdown: paired `s16`

`papa2`:
- `filter_and_trim`: `11.0s`
- `learn_errors_fwd`: `66.5s`
- `learn_errors_rev`: `33.6s`
- `derep_fwd`: `1.38s`
- `derep_rev`: `1.01s`
- `dada_fwd`: `10.3s`
- `dada_rev`: `4.02s`
- `merge_pairs`: `1.17s`

official R `dada2`:
- `filter_and_trim`: `10.0s`
- `learn_errors_fwd`: `65.9s`
- `learn_errors_rev`: `31.2s`
- `derep_fwd`: `2.79s`
- `derep_rev`: `1.82s`
- `dada_fwd`: `9.30s`
- `dada_rev`: `3.52s`
- `merge_pairs`: `1.47s`

Main conclusion:
- the wall time is dominated by `learnErrors`
- final denoising with a fixed learned error model is not the main bottleneck
- wrapper work like merge and seqtab handling is not the current limiter

### Internal `learn_errors` profiling

We then profiled `learn_errors` round by round on the forward and reverse reads of the same representative case.

Forward `learn_errors`:
- 6 self-consistency rounds
- total about `41.0s`
- per-round LOESS/error update: about `0.006s`

Reverse `learn_errors`:
- 6 rounds
- total about `18.4s`
- per-round LOESS/error update: again about `0.006s`

So:
- LOESS is negligible
- the time is in repeated denoising work during self-consistency

### Native core profiling inside `learn_errors`

For the stable forward rounds, the native timings were approximately:

| Round | Compare | Shuffle | P-update | Core | Final subs | Total |
|---|---:|---:|---:|---:|---:|---:|
| 2 | `4.56s` | `0.47s` | `0.05s` | `5.18s` | `2.14s` | `7.97s` |
| 3 | `4.60s` | `0.47s` | `0.05s` | `5.21s` | `2.15s` | `8.02s` |
| 4 | `4.66s` | `0.47s` | `0.05s` | `5.28s` | `2.16s` | `8.10s` |
| 5 | `4.67s` | `0.46s` | `0.05s` | `5.27s` | `2.16s` | `8.10s` |

Reverse rounds showed the same structure on a smaller scale.

Main conclusion:
- `compare` is the dominant native hotspot
- `shuffle` matters, but much less
- `p_update` is negligible
- `final_subs` is also expensive

This is an important result: `learn_errors` is not slow because of smoothing or Python overhead. It is slow because it repeatedly runs the heavy denoising core, and the denoising core itself is dominated by comparison work.

## What we tried and what did not help

We attempted a narrow `learn_errors` fast path that returned only transition matrices and skipped some output construction.

That did not produce a real speedup.

Why:
- the transition matrix is currently built from the final alignment objects
- it also depends on the final `correct` flags from the p-value pass
- so the obvious API-level shortcut does not remove the real work

That experiment was reverted.

## Current conclusions

### What `papa2` already demonstrates

- On a coherent paired-end cohort, `papa2` is substantially faster than official R `dada2`.
- The current best real-data result is roughly `6x` to `8x` faster wall time with materially lower memory use.
- On broader mixed real-world pools, `papa2` still wins on memory and remains somewhat faster, but the runtime gap compresses sharply.

### What matters next

If the goal is to push runtime much lower on large real studies, the main target is now clear:
- optimize or redesign `learn_errors`

More specifically:
- the hot path is the repeated denoising work inside self-consistency
- inside that, the main native cost is `compare`
- a secondary but still meaningful cost is the final alignment reconstruction used to build transition counts

### What does not appear to be the next bottleneck

- LOESS
- `merge_pairs`
- derep alone
- taxonomy assignment
- Python-level orchestration by itself

## Benchmark caveats

These results are real and useful, but they are not the last word.

Important caveats:
- the broad 515F pool is less coherent than the curated 36-sample run
- some mixed-cohort runs included zero-pass filtered samples, which complicates strict apples-to-apples interpretation
- the overnight ladder was stopped once the speedup trend was clear, so the widest scales were not fully completed
- chart generation and denser benchmark sweeps are still worth doing once the next round of optimization targets is chosen

## Related materials

- benchmark harness overview: [README.md](./README.md)
- older draft benchmark writeup used as structural reference:
  - `rec3141/microscape` `nextflow-pipeline`
  - [`BENCHMARKS.md`](https://github.com/rec3141/microscape/blob/nextflow-pipeline/BENCHMARKS.md)
