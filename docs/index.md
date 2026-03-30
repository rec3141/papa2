# papa2

**Python-first amplicon denoising derived from DADA2.**

papa2 is a Python package that implements the DADA2 (Divisive Amplicon Denoising Algorithm 2) workflow entirely in Python, backed by a standalone C/C++ core (`libpapa2.so`). It produces byte-identical results to R's `dada2` package and is designed for researchers who prefer working in Python ecosystems.

## What is DADA2?

DADA2 is a high-resolution method for correcting and inferring the true sequences in amplicon sequencing data, resolving differences of as little as one nucleotide. Instead of clustering reads into OTUs, it produces **Amplicon Sequence Variants (ASVs)** — exact biological sequences representing true sample diversity.

papa2 brings this proven algorithm to Python with no R dependency.

## Features

- **Byte-identical results** to R's `dada2` package — [20/20 parity tests pass](https://github.com/rec3141/papa2/blob/main/tests/compare_with_r.py)
- **Complete DADA2 port** — all 37 R functions have Python equivalents
- **Full workflow** in Python: filtering, dereplication, error learning, denoising, paired-end merging, chimera removal, sequence tables, taxonomy, and visualization
- **LOESS error modeling** that faithfully mirrors R's `loessErrfun`
- **Bayesian taxonomy classifier** (`assign_taxonomy`) matching R's `assignTaxonomy`
- **Interactive Sankey diagrams** for read tracking through the pipeline
- **Parallel processing** across samples via `ProcessPoolExecutor`
- Drop-in replacement for `dada2` R scripts translated to Python

## Status

The current project is CPU-only. GPU support is planned.

## Provenance

papa2 was split out of development work in:

- Fork: [rec3141/dada2](https://github.com/rec3141/dada2)
- Branch: [gpu-python](https://github.com/rec3141/dada2/tree/gpu-python)

That work branched from upstream DADA2 at commit [`72da7700`](https://github.com/benjjneb/dada2/commit/72da7700b58290e40cdce4b0856314aecf2b9dc4) (2026-02-13).

Upstream: [benjjneb/dada2](https://github.com/benjjneb/dada2)

## Project Layout

```
papa2/
├── papa2/           Python package
│   ├── dada.py      Core denoising: dada(), learn_errors()
│   ├── filter.py    Quality filtering: filter_and_trim()
│   ├── io.py        FASTQ I/O: derep_fastq()
│   ├── error.py     Error models: loess_errfun(), pacbio_errfun()
│   ├── taxonomy.py  Taxonomic classification: assign_taxonomy()
│   ├── paired.py    Paired-end merging: merge_pairs()
│   ├── chimera.py   Chimera removal: remove_bimera_denovo()
│   ├── utils.py     Utilities: sequence tables, plotting, QC
│   └── _cdada.py    ctypes bindings to libpapa2.so
├── src/             Standalone C/C++ core
├── tests/           Parity tests against R dada2
└── libpapa2.so      Compiled shared library
```

## Quick Smoke Test

```python
import papa2
print(papa2.__version__)
```

## Development Setup

```bash
conda env create -f environment.yml
conda activate papa2-dev
make clean libpapa2.so
pytest
```

## License

See [LICENSE](https://github.com/rec3141/papa2/blob/main/LICENSE).
