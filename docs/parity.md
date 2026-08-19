# Parity & Performance

papa2's contract is simple: **the same inputs produce the same outputs as
R's dada2**, byte for byte — just faster. This page documents what that
guarantee covers, how it is verified, the few deliberate exceptions, and
the performance you can expect.

The parity target is **dada2 1.40** (Bioconductor 3.21, R 4.6).

---

## What "byte-identical" covers

Every claim below is enforced by tests against a local R installation
(`tests/r_reference.R` + `tests/compare_with_r.py`), by R-free frozen
regression tests (`tests/test_r_parity_regressions.py`), and was verified
on real MiSeq data — a 6-sample and a 120-sample 16S/18S seawater
dataset — before release.

| Stage | Guarantee |
|-------|-----------|
| `filter_and_trim` / `fastq_filter` / `fastq_paired_filter` | Filtered FASTQ files are byte-identical to R's after decompression (252/252 files across both benchmark datasets), including `trunc_q`/`trim_left` ordering, `trunc_len` counted from the original 5' end, `orient_fwd` read swapping, `min_q` strictness, PhiX screening, `match_ids` re-pairing with chunk-boundary carry, `quality_type` auto-detection, and removal of output files when nothing passes |
| `is_phix` / `match_ref` | Exact port of R's `C_matchRef`: full circular PhiX genome, per-strand hit counts, non-overlapping kmer skip |
| `derep_fastq` | Identical unique sequences, abundances, quality means, and read maps (including Phred+64 inputs) |
| `learn_errors` | Identical sample selection (base counting and stop condition); error matrices agree to ~1e-15 (see exceptions) |
| `dada` | Identical ASVs, abundances, and per-unique maps — including `pool=TRUE`, `pool="pseudo"`, and `priors` |
| `merge_pairs` | Identical merged sequences, `prefer` chosen from denoised n0 like R, rejects excluded unless `return_rejects=True` |
| `make_sequence_table` | Identical tables including tied-abundance column order (stable sort, like R's `order()`) |
| `remove_bimera_denovo` | Identical for consensus, pooled, and per-sample methods, with R's per-method defaults (1.5/2 consensus; 2/8 otherwise) |
| `remove_primers` | Byte-identical output: anywhere-in-read matching, IUPAC degenerate codes, reverse-complement re-orientation |
| `assign_taxonomy` | Identical assignments; with `seed=N`, bootstrap values reproduce an R session that calls `set.seed(N)` first (exact port of R's RNG) |

## Deliberate exceptions

1. **LOESS error matrices** agree with R to ~1e-15 rather than bit-for-bit
   (numpy/LAPACK QR/SVD vs R's internal loess solvers). Every downstream
   result is byte-identical on all datasets tested.
2. **Taxonomy bootstrap tie-breaking** uses an unseeded RNG *in R itself*
   — R differs from its own reruns by a handful of bootstrap cells. papa2
   reproduces the behaviour (and therefore matches R only to within R's
   own noise on those cells). Tracked in
   [#4](https://github.com/rec3141/papa2/issues/4), reported upstream as
   [benjjneb/dada2#2211](https://github.com/benjjneb/dada2/issues/2211).
3. **`match_ids` with pre-CASAVA-1.8 ids** matches zero reads, because R
   strips the `#`-suffix from forward ids only. Replicated bug-for-bug;
   tracked in [#3](https://github.com/rec3141/papa2/issues/3), reported
   upstream as
   [benjjneb/dada2#2210](https://github.com/benjjneb/dada2/issues/2210).
4. **`remove_primers(allow_indels=True)`** (R's `with.indels`) is not yet
   implemented and raises `NotImplementedError` rather than silently
   differing.

## Verifying parity yourself

With R and dada2 installed:

```bash
Rscript tests/r_reference.R          # writes tests/r_outputs/
python tests/compare_with_r.py      # 20 function-level comparisons
pytest tests/                        # includes R-free frozen regressions
```

---

## Performance

Measured on a 32-core machine against R dada2 1.40 running multithreaded,
on real 2×151 MiSeq data:

| Workload | R dada2 | papa2 | Speedup |
|----------|---------|-------|---------|
| 6 samples end-to-end (filter → chimera) | 284 s | 156 s | 1.8× |
| 120 samples end-to-end | 1019 s | 512 s | 2.0× |
| `filterAndTrim`, single-threaded | 69 s | 10 s | 7× |
| `assignTaxonomy`, SILVA 138.1 × 2,000 ASVs | 71 s | 15 s | 4.7× |
| `assignTaxonomy`, SILVA 138.1 × 20,697 ASVs | ~35 min | 50 s | ~40× |

Where it comes from:

- an `-O2` C core (the algorithmic heart, shared with R, compiled well)
- vectorized numpy filtering with isal-accelerated gzip
- OpenMP in the comparison loop, chimera detection, taxonomy, and PhiX
  screening
- worker-pool parallelism across samples, composed with the OpenMP
  threading (below)
- parallel list APIs: `derep_fastq([...])`, `merge_pairs` on lists, and
  pooled dereplication in `learn_errors`

## Threading

papa2 composes two levels of parallelism and splits cores between them
automatically:

- **sample-level worker processes** (spawn-safe, via loky)
- **OpenMP threads** inside the C core, per worker

With 32 cores: 120 samples → 32 workers × 1 thread; 6 samples → 6 workers
× 5 threads; one pooled sample → 1 worker × 32 threads.

Environment variables:

| Variable | Meaning |
|----------|---------|
| `DADA2_CORES` | Total cores to plan for. Set this under containers or batch schedulers, where `os.cpu_count()` overreports the real allocation. |
| `DADA2_WORKERS` | Sample-level worker processes (0 = auto). |
| `DADA2_OMP_THREADS` | OpenMP threads per worker (0 = auto: cores split evenly across workers). |
| `OMP_NUM_THREADS` | Caps standalone OpenMP regions (chimera tables, taxonomy, PhiX). Read by the OpenMP runtime **at import time** — set it before importing papa2, or call `papa2.set_num_threads(n)` afterwards. |

!!! tip "Containers and schedulers"
    Inside a task that was allocated `N` CPUs, do this before importing
    papa2:

    ```python
    os.environ.setdefault("DADA2_CORES", str(N))
    os.environ.setdefault("OMP_NUM_THREADS", str(N))
    import papa2
    papa2.set_num_threads(N)
    ```

    The CLI equivalents are `papa2 filter-trim --threads N` and
    `papa2 assign-taxonomy --threads N`.
