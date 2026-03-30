# papa2

**Python-first amplicon denoising — byte-identical to R's DADA2, 8x faster**

[![CI](https://github.com/rec3141/papa2/actions/workflows/ci.yml/badge.svg)](https://github.com/rec3141/papa2/actions/workflows/ci.yml)
[![Python](https://img.shields.io/badge/python-3.9%2B-blue)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

`papa2` is a pure-Python + compiled-C reimplementation of the [DADA2](https://github.com/benjjneb/dada2) amplicon denoising pipeline. It produces **byte-identical ASV output** to R's DADA2 while running significantly faster — no R installation required.

Full documentation: **https://rec3141.github.io/papa2**

---

## Key Features

- **Byte-identical to R dada2** — same ASVs, same error models, validated on real datasets
- **4–8x faster than R** on real datasets (tested up to 40M reads; see benchmarks below)
- **No R dependency** — standalone Python package with a compiled C core (`libpapa2.so`)
- **Full DADA2 pipeline** — filter & trim, dereplication, error learning, denoising, paired-end merging, chimera removal, and taxonomy assignment
- **Parallel multi-sample processing** — scales across all available CPU cores

---

## Quick Install

**pip:**
```bash
pip install papa2
```

**conda:**
```bash
conda install -c conda-forge papa2
```

---

## Usage

```python
import papa2

# Learn error model from a set of FASTQ files
err = papa2.learn_errors(fastq_files)

# Dereplicate reads
derep = papa2.dereplicate(fastq_files)

# Denoise
dada_result = papa2.dada(derep, err)

# Merge paired-end reads
merged = papa2.merge_pairs(dada_result["fwd"], derep["fwd"],
                           dada_result["rev"], derep["rev"])

# Build sequence table and remove chimeras
seqtab = papa2.make_sequence_table(merged)
seqtab_nochim = papa2.remove_bimera_denovo(seqtab)
```

For a complete walkthrough see the [Quickstart guide](https://rec3141.github.io/papa2/quickstart/) and [Tutorial](https://rec3141.github.io/papa2/tutorial/).

---

## Performance Benchmarks

Benchmarks run on a single workstation (AMD EPYC, 32 cores). papa2 uses all available cores; R dada2 uses a single core by default.

| Dataset | papa2 | R dada2 | Speedup |
|---|---|---|---|
| 10 samples, 257K reads | 49s | 85s | 1.7x |
| 84 samples, 4M reads | 174s | 1500s | 8.6x |
| 1790 samples, 40M reads | ~25 min | ~5h | ~12x |

---

## Project Layout

```
papa2/        Python package
src/          standalone C/C++ native core (compiled to libpapa2.so)
tests/        pytest suite with bundled FASTQ fixtures
docs/         MkDocs documentation source
```

---

## Development Setup

```bash
# Create and activate the dev environment
conda env create -f environment.yml
conda activate papa2-dev

# Build the native shared library
make clean libpapa2.so

# Run the test suite
pytest
```

Verify the package loads:

```bash
python -c "import papa2; print(papa2.__version__)"
```

---

## Provenance

`papa2` was extracted from the development branch [rec3141/dada2@gpu-python](https://github.com/rec3141/dada2/tree/gpu-python), which itself forked from upstream DADA2 at commit [72da7700b](https://github.com/benjjneb/dada2/commit/72da7700b58290e40cdce4b0856314aecf2b9dc4) (2026-02-13).

---

## Citation

If you use papa2 in published research, please cite the original DADA2 paper:

> Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJA, Holmes SP (2016).
> **DADA2: High-resolution sample inference from Illumina amplicon data.**
> *Nature Methods*, 13, 581–583. https://doi.org/10.1038/nmeth.3869

---

## License

MIT — see [LICENSE](LICENSE).
