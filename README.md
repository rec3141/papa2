# papa2

**Python-first amplicon denoising — byte-identical to R's DADA2**

[![CI](https://github.com/rec3141/papa2/actions/workflows/ci.yml/badge.svg)](https://github.com/rec3141/papa2/actions/workflows/ci.yml)
[![Docs](https://github.com/rec3141/papa2/actions/workflows/docs.yml/badge.svg)](https://rec3141.github.io/papa2/)
[![Container](https://github.com/rec3141/papa2/actions/workflows/container.yml/badge.svg)](https://github.com/rec3141/papa2/pkgs/container/papa2)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue)](https://www.python.org/)
[![License: LGPL v3](https://img.shields.io/badge/License-LGPL_v3-blue.svg)](LICENSE)

`papa2` is a complete Python port of the [DADA2](https://github.com/benjjneb/dada2) amplicon denoising pipeline. All 37 R functions have Python equivalents, producing **byte-identical results** ([20/20 parity tests pass](tests/compare_with_r.py)) with no R dependency.

Full documentation: **https://rec3141.github.io/papa2**

---

## Key Features

- **Byte-identical to R dada2** — same ASVs, same error models, same quality filtering, [validated against R on every function](tests/compare_with_r.py)
- **Complete port** — all 37 user-facing R functions: `filter_and_trim`, `derep_fastq`, `learn_errors`, `dada`, `merge_pairs`, `make_sequence_table`, `remove_bimera_denovo`, `assign_taxonomy`, and more
- **No R dependency** — standalone Python package with a compiled C/C++ core (`libpapa2.so`)
- **Parallel multi-sample processing** — scales across CPU cores via `ProcessPoolExecutor`
- **Interactive visualisation** — Sankey diagrams for read tracking, error rate plots, quality profiles
- **Container-ready** — Docker and Apptainer (HPC) images on GHCR

---

## Quick Install

### Container (recommended for HPC)

```bash
# Docker
docker pull ghcr.io/rec3141/papa2:latest
docker run -v $(pwd):/data ghcr.io/rec3141/papa2 python3 my_script.py

# Apptainer (HPC)
apptainer pull papa2.sif docker://ghcr.io/rec3141/papa2:latest
apptainer exec papa2.sif python3 my_script.py
```

### From source

```bash
git clone https://github.com/rec3141/papa2.git
cd papa2
conda env create -f environment.yml
conda activate papa2-dev
make libpapa2.so
pip install -e .
```

---

## Usage

```python
import papa2

# Forward and reverse FASTQ files
fwd_files = ["sample1_R1.fastq.gz", "sample2_R1.fastq.gz"]
rev_files = ["sample1_R2.fastq.gz", "sample2_R2.fastq.gz"]

# 1. Filter and trim
papa2.filter_and_trim(
    fwd_files, filt_fwd,
    rev=rev_files, filt_rev=filt_rev,
    trunc_len=(240, 200), max_ee=(2, 2), rm_phix=True,
)

# 2. Learn error rates
errF = papa2.learn_errors(filt_fwd)
errR = papa2.learn_errors(filt_rev)

# 3. Dereplicate, denoise, merge
derepFs = [papa2.derep_fastq(f) for f in filt_fwd]
derepRs = [papa2.derep_fastq(f) for f in filt_rev]
dadaFs = papa2.dada(derepFs, err=errF)
dadaRs = papa2.dada(derepRs, err=errR)
mergers = [papa2.merge_pairs(dF, drF, dR, drR)
           for dF, drF, dR, drR in zip(dadaFs, derepFs, dadaRs, derepRs)]

# 4. Sequence table and chimera removal
seqtab = papa2.make_sequence_table(mergers)
seqtab_nochim = papa2.remove_bimera_denovo(seqtab)

# 5. Taxonomy
taxa = papa2.assign_taxonomy(seqtab_nochim["seqs"], "silva_nr99_v138.1_train_set.fa.gz")

# 6. Visualise
track = papa2.track_reads(dereps=derepFs, dadas=dadaFs, mergers=mergers,
                          seqtab=seqtab, seqtab_nochim=seqtab_nochim, taxa=taxa)
papa2.plot_sankey(track, output="read_tracking.html")
```

See the [Quickstart](https://rec3141.github.io/papa2/quickstart/), [Tutorial](https://rec3141.github.io/papa2/tutorial/), and [Big Data Workflow](https://rec3141.github.io/papa2/bigdata/) for complete walkthroughs.

---

## Project Layout

```
papa2/           Python package
  filter.py        filter_and_trim, fastq_filter, fastq_paired_filter
  io.py            derep_fastq
  dada.py          dada, learn_errors
  error.py         loess_errfun, pacbio_errfun, make_binned_qual_errfun
  taxonomy.py      assign_taxonomy
  paired.py        merge_pairs
  chimera.py       remove_bimera_denovo
  utils.py         sequence tables, plotting, QC, taxonomy, export
  _cdada.py        ctypes bindings to libpapa2.so
src/             Standalone C/C++ core
tests/           Parity tests against R dada2
docs/            MkDocs documentation source
```

---

## Provenance

`papa2` was extracted from [rec3141/dada2@gpu-python](https://github.com/rec3141/dada2/tree/gpu-python), which forked from upstream DADA2 at commit [72da7700b](https://github.com/benjjneb/dada2/commit/72da7700b58290e40cdce4b0856314aecf2b9dc4) (2026-02-13).

---

## Citation

If you use papa2 in published research, please cite the original DADA2 paper:

> Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJA, Holmes SP (2016).
> **DADA2: High-resolution sample inference from Illumina amplicon data.**
> *Nature Methods*, 13, 581–583. https://doi.org/10.1038/nmeth.3869

---

## License

LGPL-3.0 — see [LICENSE](LICENSE). This matches the upstream [DADA2 license](https://github.com/benjjneb/dada2/blob/master/LICENSE).
