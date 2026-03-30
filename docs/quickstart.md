# Quickstart

## Installation

### Prerequisites

papa2 requires a compiled `libdada2.so` shared library. Build it from source:

```bash
git clone https://github.com/rec3141/dada2.git
cd dada2
git checkout gpu-python
```

### Build the C/C++ core

```bash
# Using the conda environment (recommended)
mamba env create -f environment.yml
conda activate dada2-dev
make libdada2.so
```

Or build directly with:

```bash
make libdada2.so
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

This example walks through the core DADA2 workflow on a single pair of FASTQ files
(forward and reverse reads from a 16S amplicon experiment).

```python
import papa2

# --- 1. Dereplicate FASTQ files ---
# derep_fastq() collapses identical reads and tracks abundances
derepF = papa2.derep_fastq("sample_R1.fastq.gz", verbose=True)
derepR = papa2.derep_fastq("sample_R2.fastq.gz", verbose=True)

# Each derep dict contains:
#   seqs        - unique sequences (sorted by abundance)
#   abundances  - read counts per unique
#   quals       - average quality scores (n_uniques x max_seqlen)
#   map         - per-read index into the unique sequence list

# --- 2. Learn error rates ---
# learn_errors() uses LOESS smoothing on a set of training FASTQ files
# to build a (16 x nqual) error rate matrix (one entry per base transition
# per quality score).
errF = papa2.learn_errors(["sample_R1.fastq.gz"], verbose=True)
errR = papa2.learn_errors(["sample_R2.fastq.gz"], verbose=True)

# --- 3. Denoise with DADA ---
# dada() infers the true biological sequences (ASVs) and returns a
# result dict containing 'denoised' (seq -> abundance mapping).
dadaF = papa2.dada(derepF, err=errF, verbose=True)
dadaR = papa2.dada(derepR, err=errR, verbose=True)

# --- 4. Merge paired-end reads ---
# merge_pairs() aligns denoised forward and reverse reads and returns
# merged amplicon sequences, discarding non-overlapping pairs.
mergers = papa2.merge_pairs(dadaF, derepF, dadaR, derepR, verbose=True)

# mergers is a list of dicts sorted by abundance, each with:
#   sequence   - merged amplicon string
#   abundance  - supporting read count
#   accept     - True if overlap and mismatch criteria were met

# --- 5. Remove chimeras ---
# Build a per-sample sequence table first, then filter bimeras.
seqtab = papa2.make_sequence_table([mergers])

# remove_bimera_denovo() uses consensus detection across samples by default.
seqtab_nochim = papa2.remove_bimera_denovo(seqtab, method="consensus", verbose=True)

print(f"ASVs before chimera removal: {len(seqtab['seqs'])}")
print(f"ASVs after chimera removal:  {len(seqtab_nochim['seqs'])}")

# --- 6. Inspect results ---
# seqtab_nochim['table'] is a numpy array (samples x ASVs)
# seqtab_nochim['seqs']  is the list of ASV sequences
table = seqtab_nochim["table"]
seqs  = seqtab_nochim["seqs"]
print(f"Sequence table shape: {table.shape}")
print(f"First ASV: {seqs[0][:60]}...")
```

---

## Multiple Samples

To process many samples in parallel, pass lists of paths or derep dicts:

```python
import glob
import papa2

fwd_files = sorted(glob.glob("data/*_R1.fastq.gz"))
rev_files  = sorted(glob.glob("data/*_R2.fastq.gz"))

# Learn errors from a representative subset of forward reads
errF = papa2.learn_errors(fwd_files[:10], verbose=True)
errR = papa2.learn_errors(rev_files[:10], verbose=True)

# Dereplicate all samples
derepFs = [papa2.derep_fastq(f) for f in fwd_files]
derepRs = [papa2.derep_fastq(f) for f in rev_files]

# Denoise all forward and reverse reads (runs in parallel via DADA2_WORKERS)
dadaFs = papa2.dada(derepFs, err=errF, verbose=False)
dadaRs = papa2.dada(derepRs, err=errR, verbose=False)

# Merge each sample
mergers = [
    papa2.merge_pairs(dF, drF, dR, drR)
    for dF, drF, dR, drR in zip(dadaFs, derepFs, dadaRs, derepRs)
]

# Build joint sequence table across all samples
seqtab = papa2.make_sequence_table(mergers)

# Remove chimeras
seqtab_nochim = papa2.remove_bimera_denovo(seqtab, verbose=True)
```

!!! tip "Parallelism"
    Set the environment variable `DADA2_WORKERS` to control how many parallel
    workers are used when denoising multiple samples:
    ```bash
    export DADA2_WORKERS=8
    ```
    Set `OMP_NUM_THREADS=1` before importing papa2 for best multi-sample
    performance (avoids contention between Python-level and OpenMP parallelism).
