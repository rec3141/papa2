# API Reference

This page documents all public functions and objects in papa2. Each module's
docstrings are rendered directly from the source.

---

## papa2.dada — Core Denoising

The main denoising module. Provides the high-level `dada()` and `learn_errors()`
functions as well as the global options dictionary.

### DADA_OPTS

```python
papa2.DADA_OPTS
```

Global dictionary of algorithmic parameters used by `dada()`. Modify with
`set_dada_opt()` or pass keyword arguments directly to `dada()`.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `OMEGA_A` | `1e-40` | Significance threshold for accepting new ASVs |
| `OMEGA_P` | `1e-4` | Significance threshold for prior-guided detection |
| `OMEGA_C` | `1e-40` | Significance threshold for combining ASVs |
| `DETECT_SINGLETONS` | `False` | Detect singleton ASVs |
| `USE_KMERS` | `True` | Use k-mer screen before alignment |
| `KDIST_CUTOFF` | `0.42` | K-mer distance cutoff |
| `MAX_CONSIST` | `10` | Maximum self-consistency iterations |
| `MATCH` | `5` | NW match score |
| `MISMATCH` | `-4` | NW mismatch penalty |
| `GAP_PENALTY` | `-8` | NW gap penalty |
| `BAND_SIZE` | `16` | Banded alignment width |
| `VECTORIZED_ALIGNMENT` | `True` | Use vectorized (SSE) alignment |
| `MAX_CLUST` | `0` | Maximum clusters (0 = unlimited) |
| `MIN_FOLD` | `1` | Minimum fold-abundance for parent |
| `MIN_HAMMING` | `1` | Minimum Hamming distance from parent |
| `MIN_ABUNDANCE` | `1` | Minimum read abundance to consider |
| `USE_QUALS` | `True` | Incorporate quality scores |
| `HOMOPOLYMER_GAP_PENALTY` | `None` | Override gap penalty in homopolymer runs |
| `SSE` | `2` | SSE level (0=off, 1=SSE2, 2=SSE4.1) |
| `GAPLESS` | `True` | Prefer gapless alignments |
| `GREEDY` | `True` | Use greedy clustering |

### dada

::: papa2.dada.dada
    options:
      heading_level: 3

### learn_errors

::: papa2.dada.learn_errors
    options:
      heading_level: 3

### set_dada_opt

::: papa2.dada.set_dada_opt
    options:
      heading_level: 3

### get_dada_opt

::: papa2.dada.get_dada_opt
    options:
      heading_level: 3

---

## papa2.io — FASTQ I/O

FASTQ reading and dereplication.

### derep_fastq

::: papa2.io.derep_fastq
    options:
      heading_level: 3

---

## papa2.filter — Filtering and Trimming

Quality filtering, trimming, and PhiX removal for FASTQ files.

### filter_and_trim

::: papa2.filter.filter_and_trim
    options:
      heading_level: 3

### fastq_filter

::: papa2.filter.fastq_filter
    options:
      heading_level: 3

### fastq_paired_filter

::: papa2.filter.fastq_paired_filter
    options:
      heading_level: 3

---

## papa2.error — Error Models

Functions for estimating and manipulating the DADA2 error rate matrix.

### loess_errfun

::: papa2.error.loess_errfun
    options:
      heading_level: 3

### noqual_errfun

::: papa2.error.noqual_errfun
    options:
      heading_level: 3

### inflate_err

::: papa2.error.inflate_err
    options:
      heading_level: 3

### pacbio_errfun

::: papa2.error.pacbio_errfun
    options:
      heading_level: 3

### make_binned_qual_errfun

::: papa2.error.make_binned_qual_errfun
    options:
      heading_level: 3

### get_initial_err

::: papa2.error.get_initial_err
    options:
      heading_level: 3

---

## papa2.paired — Paired-End Merging

High-level paired-read merging, mirroring R's `mergePairs()`.

### merge_pairs

::: papa2.paired.merge_pairs
    options:
      heading_level: 3

---

## papa2.chimera — Chimera Removal

Chimera detection and removal functions matching R's dada2 interface.

### remove_bimera_denovo

::: papa2.chimera.remove_bimera_denovo
    options:
      heading_level: 3

### is_bimera_denovo

::: papa2.chimera.is_bimera_denovo
    options:
      heading_level: 3

---

## papa2.taxonomy — Taxonomic Classification

Bayesian k-mer classifier for taxonomic assignment.

### assign_taxonomy

::: papa2.taxonomy.assign_taxonomy
    options:
      heading_level: 3

---

## papa2.utils — Utilities

Utility functions for taxonomy assignment, sequence table operations, quality
profiling, FASTA I/O, PhiX detection, and sequence complexity analysis.

### make_sequence_table

::: papa2.utils.make_sequence_table
    options:
      heading_level: 3

### assign_species

::: papa2.utils.assign_species
    options:
      heading_level: 3

### add_species

::: papa2.utils.add_species
    options:
      heading_level: 3

### collapse_no_mismatch

::: papa2.utils.collapse_no_mismatch
    options:
      heading_level: 3

### plot_quality_profile

::: papa2.utils.plot_quality_profile
    options:
      heading_level: 3

### uniquesto_fasta

::: papa2.utils.uniquesto_fasta
    options:
      heading_level: 3

### write_fasta

::: papa2.utils.write_fasta
    options:
      heading_level: 3

### is_phix

::: papa2.utils.is_phix
    options:
      heading_level: 3

### seq_complexity

::: papa2.utils.seq_complexity
    options:
      heading_level: 3

### get_sequences

::: papa2.utils.get_sequences
    options:
      heading_level: 3

### get_uniques

::: papa2.utils.get_uniques
    options:
      heading_level: 3

### get_errors

::: papa2.utils.get_errors
    options:
      heading_level: 3

### merge_sequence_tables

::: papa2.utils.merge_sequence_tables
    options:
      heading_level: 3

### nwhamming

::: papa2.utils.nwhamming
    options:
      heading_level: 3

### is_shift_denovo

::: papa2.utils.is_shift_denovo
    options:
      heading_level: 3

### plot_errors

::: papa2.utils.plot_errors
    options:
      heading_level: 3

### plot_complexity

::: papa2.utils.plot_complexity
    options:
      heading_level: 3

### plot_sankey

::: papa2.utils.plot_sankey
    options:
      heading_level: 3

### track_reads

::: papa2.utils.track_reads
    options:
      heading_level: 3

### remove_primers

::: papa2.utils.remove_primers
    options:
      heading_level: 3

---

## papa2._cdada — C Bindings

Low-level ctypes bindings to the `libpapa2.so` shared library. Most users should
not call these directly — they are used internally by the higher-level functions
above.

### run_dada

::: papa2._cdada.run_dada
    options:
      heading_level: 3

### run_taxonomy

::: papa2._cdada.run_taxonomy
    options:
      heading_level: 3

### nwalign

::: papa2._cdada.nwalign
    options:
      heading_level: 3

### eval_pair

::: papa2._cdada.eval_pair
    options:
      heading_level: 3

### pair_consensus

::: papa2._cdada.pair_consensus
    options:
      heading_level: 3

### rc

::: papa2._cdada.rc
    options:
      heading_level: 3

