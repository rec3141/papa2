# papa2 Benchmarks

This document summarizes the current real-data benchmark results for `papa2` against the official R `dada2` package.

It complements the harness overview in [README.md](./README.md).

## Baselines

Python:
- `papa2`

R:
- official upstream `dada2` release `v1.26`
- commit `f1e72fb`
- release date `2022-11-03`

The R comparison was run against a clean official release, not a modified development tree.

## Datasets

Two real paired-end 16S datasets were used.

### Curated paired cohort

Properties:
- 36 complete `R1`/`R2` cutadapt pairs
- 515F-806RB amplicons
- coherent single-run dataset

### Broad 515F pool

Properties:
- larger mixed pool of real cutadapt outputs
- filtered before benchmarking to remove:
  - 18S samples
  - very small samples
- used for wider scaling checks

## Benchmark workflow

The comparisons used a real paired-end workflow:

1. filter and trim reads
2. learn forward and reverse error models
3. denoise forward and reverse reads
4. merge pairs

Filtering was included because the cutadapt outputs still contained some ambiguous bases, and neither engine can denoise those reads directly.

## Results

### Curated paired cohort, 4 threads

One paired sample:

| Engine | Wall time | Max RSS |
|---|---:|---:|
| `papa2` | `3.84s` | `140 MB` |
| `dada2` | `34.51s` | `1.11 GB` |

Ten paired samples:

| Engine | Wall time | Max RSS |
|---|---:|---:|
| `papa2` | `50.43s` | `484 MB` |
| `dada2` | `386.75s` | `1.76 GB` |

Takeaways:
- on one large paired sample, `papa2` was about `9x` faster than R
- on ten paired samples, `papa2` was about `7.7x` faster
- memory use was also much lower

### Curated paired cohort, 16 threads

One paired sample:

| Engine | Wall time | Max RSS |
|---|---:|---:|
| `papa2` | `2.61s` | `183 MB` |
| `dada2` | `21.24s` | `1.12 GB` |

Ten paired samples:

| Engine | Wall time | Max RSS |
|---|---:|---:|
| `papa2` | `30.66s` | `521 MB` |
| `dada2` | `194.3s` | `1.82 GB` |

All 36 paired samples:

| Engine | Wall time | Max RSS |
|---|---:|---:|
| `papa2` | `103.08s` | `1.38 GB` |
| `dada2` | `637.0s` | `2.61 GB` |

Takeaways:
- on this coherent paired cohort, `papa2` was roughly `6x` to `8x` faster
- memory stayed materially lower as well
- this is the strongest current benchmark result

### Broad 515F pool, 16 threads

| Samples | `papa2` wall | `dada2` wall | Speedup | `papa2` RSS | `dada2` RSS |
|---|---:|---:|---:|---:|---:|
| 1 | `2.17s` | `7.25s` | `3.3x` | `96,776 KB` | `860,072 KB` |
| 2 | `2.81s` | `8.29s` | `3.0x` | `97,540 KB` | `923,464 KB` |
| 4 | `37.49s` | `42.97s` | `1.15x` | `207,672 KB` | `1,283,228 KB` |
| 8 | `48.82s` | `58.50s` | `1.20x` | `192,880 KB` | `1,234,348 KB` |
| 16 | `80.91s` | `90.09s` | `1.11x` | `269,836 KB` | `1,355,308 KB` |
| 32 | `195.81s` | `222.15s` | `1.13x` | `433,076 KB` | `1,688,928 KB` |

Takeaways:
- memory remained clearly better for `papa2`
- the runtime gap compressed sharply on this broader mixed pool
- that suggests the dominant cost at larger scale is not simply the final fixed-error denoising pass

## Profiling summary

To explain the flattening on the broader benchmark, we profiled a representative paired `s16` case stage by stage.

### Stage breakdown

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
- `learnErrors` dominates the wall time in both engines
- final denoising with a fixed learned error model is not the main bottleneck

### Internal `learn_errors` profiling

Round-by-round profiling showed:
- LOESS is negligible
- the cost is repeated denoising work inside self-consistency

Stable forward rounds were approximately:

| Round | Compare | Shuffle | P-update | Core | Final subs | Total |
|---|---:|---:|---:|---:|---:|---:|
| 2 | `4.56s` | `0.47s` | `0.05s` | `5.18s` | `2.14s` | `7.97s` |
| 3 | `4.60s` | `0.47s` | `0.05s` | `5.21s` | `2.15s` | `8.02s` |
| 4 | `4.66s` | `0.47s` | `0.05s` | `5.28s` | `2.16s` | `8.10s` |
| 5 | `4.67s` | `0.46s` | `0.05s` | `5.27s` | `2.16s` | `8.10s` |

Main conclusion:
- `compare` is the dominant native hotspot
- `shuffle` matters, but much less
- `p_update` is negligible
- `final_subs` is also expensive

This means `learn_errors` is not slow because of smoothing or Python overhead. It is slow because it repeatedly runs the heavy denoising core.

## What this means

- On a coherent paired-end cohort, `papa2` is substantially faster than official R `dada2`.
- The current best real-data result is roughly `6x` to `8x` faster wall time with materially lower memory use.
- On broader mixed real-world pools, `papa2` still wins on memory and remains somewhat faster, but the runtime gap compresses sharply.

If the goal is to reduce runtime further on large real studies, the main target is now clear:
- optimize or redesign `learn_errors`

More specifically:
- the hot path is the repeated denoising work inside self-consistency
- inside that, the main native cost is `compare`
- a secondary cost is the alignment reconstruction used to build transition counts

## Failed experiments

### Python-level sample parallelism

We tested two `papa2` execution modes:

- within-sample threading
- across-sample worker parallelism

The across-sample mode did not produce a convincing win on the paired workloads tested here and was materially worse at higher thread counts.

At 16 threads:
- one paired sample: `2.90s` vs `10.05s`
- ten paired samples: `30.43s` vs `44.30s`
- peak memory also became much worse in the sample-parallel configuration

So the benchmark results in this document use the within-sample threaded mode.

### Transition-only `learn_errors` shortcut

We also tried a narrow fast path that returned only transition matrices during `learn_errors`.

That did not produce a meaningful speedup because the transition matrix still depends on:
- final alignment reconstruction
- final corrected-read flags from the p-value pass

That experiment was reverted.

## Caveats

- the broad 515F pool is less coherent than the curated 36-sample run
- some mixed-cohort runs included zero-pass filtered samples, which complicates strict apples-to-apples interpretation
- the broader overnight ladder was stopped once the trend was clear, so the widest scales were not completed

## Related

- harness overview: [README.md](./README.md)
- older draft benchmark writeup used as structural reference:
  - [`rec3141/microscape` `BENCHMARKS.md`](https://github.com/rec3141/microscape/blob/nextflow-pipeline/BENCHMARKS.md)
