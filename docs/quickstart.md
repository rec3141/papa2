# Quickstart

## Installation

### Container (Docker / Apptainer)

The fastest way to get started — no compiler or conda environment required.

```bash
# Docker
docker pull ghcr.io/rec3141/papa2:latest
docker run -v $(pwd):/data ghcr.io/rec3141/papa2 python3 my_script.py

# Apptainer (HPC)
apptainer pull papa2.sif docker://ghcr.io/rec3141/papa2:latest
apptainer exec papa2.sif python3 my_script.py
```

---

### Prerequisites

papa2 requires a compiled `libpapa2.so` shared library. Build it from source:

```bash
git clone https://github.com/rec3141/papa2.git
cd papa2
```

### Build the C/C++ core

```bash
# Using the conda environment (recommended)
mamba env create -f environment.yml
conda activate papa2-dev
make libpapa2.so
```

Or build directly with:

```bash
make libpapa2.so
```

### Install the Python package

```bash
pip install -e .
```

### Verify

```python
import papa2
print(papa2.__version__)
# 0.1.0
```

---

## Minimal Example

papa2 ships with two small test FASTQ files in `tests/data/`. This example
runs the core denoising pipeline on them — no external data needed.

```python
import papa2

# The bundled test files (forward reads from two 16S amplicon samples)
fwd_files = ["tests/data/sam1F.fastq.gz", "tests/data/sam2F.fastq.gz"]

# --- 1. Dereplicate ---
dereps = [papa2.derep_fastq(f, verbose=True) for f in fwd_files]
# Read 1500 reads, 896 unique sequences
# Read 1500 reads, 906 unique sequences

# --- 2. Learn error rates ---
err = papa2.learn_errors(fwd_files, verbose=True)
# 750000 total bases in 3000 reads from 2 samples ...

# --- 3. Denoise ---
dadas = papa2.dada(dereps, err=err, verbose=True)

# --- 4. Inspect results ---
for i, dd in enumerate(dadas):
    print(f"Sample {i+1}: {len(dd['denoised'])} ASVs, "
          f"{sum(dd['denoised'].values())} reads")

# --- 5. Build sequence table ---
seqtab = papa2.make_sequence_table(dadas)
print(f"Sequence table: {seqtab.shape[0]} samples x {seqtab.shape[1]} ASVs")
```

!!! note "Test data"
    The bundled files are small subsets (1500 reads each) of the same samples
    used in the [DADA2 tutorial](https://benjjneb.github.io/dada2/tutorial.html).
    For real analyses you'll typically have 10K–200K reads per sample.

---

## Paired-End Example

For paired-end data, dereplication and denoising are done separately on
forward and reverse reads, then merged:

```python
import papa2

# Forward and reverse FASTQ files (replace with your own paths)
fwd_files = ["sample1_R1.fastq.gz", "sample2_R1.fastq.gz"]
rev_files = ["sample1_R2.fastq.gz", "sample2_R2.fastq.gz"]

# Learn errors from forward and reverse reads separately
errF = papa2.learn_errors(fwd_files, verbose=True)
errR = papa2.learn_errors(rev_files, verbose=True)

# Dereplicate
derepFs = [papa2.derep_fastq(f) for f in fwd_files]
derepRs = [papa2.derep_fastq(f) for f in rev_files]

# Denoise
dadaFs = papa2.dada(derepFs, err=errF)
dadaRs = papa2.dada(derepRs, err=errR)

# Merge paired reads
mergers = [
    papa2.merge_pairs(dF, drF, dR, drR)
    for dF, drF, dR, drR in zip(dadaFs, derepFs, dadaRs, derepRs)
]

# Build sequence table and remove chimeras
seqtab = papa2.make_sequence_table(mergers)
seqtab_nochim = papa2.remove_bimera_denovo(seqtab, verbose=True)

print(f"Final: {seqtab_nochim['table'].shape[1]} ASVs across "
      f"{seqtab_nochim['table'].shape[0]} samples")
```

---

## Parallelism

Set `DADA2_WORKERS` to control parallel workers for multi-sample denoising:

```bash
export DADA2_WORKERS=8
```

Set `OMP_NUM_THREADS=1` before importing papa2 for best multi-sample
performance (avoids contention between Python-level and OpenMP parallelism).
