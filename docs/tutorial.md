# Tutorial

[:material-download: Download tutorial.py](scripts/tutorial.py){ .md-button }

This tutorial walks through the standard DADA2 amplicon processing workflow using
papa2's Python API. The workflow mirrors the canonical
[DADA2 R tutorial](https://benjjneb.github.io/dada2/tutorial.html) but uses Python
syntax and papa2's interface throughout.

**Expected input:** paired-end Illumina reads from a 16S (or ITS) amplicon
experiment, demultiplexed into per-sample FASTQ files.

---

## 0. Setup

```python
import os, glob
import numpy as np
import papa2

# Point to your FASTQ directory
fastq_dir = "data/fastq"

fwd_files = sorted(glob.glob(os.path.join(fastq_dir, "*_R1_001.fastq.gz")))
rev_files = sorted(glob.glob(os.path.join(fastq_dir, "*_R2_001.fastq.gz")))

# Derive sample names from file names
sample_names = [
    os.path.basename(f).replace("_R1_001.fastq.gz", "")
    for f in fwd_files
]
print(sample_names[:5])
```

---

## 1. Inspect Read Quality

Before filtering, visualise per-cycle quality to choose appropriate truncation
lengths. `plot_quality_profile()` produces a heat-map of quality scores across
read positions — the same summary as R's `plotQualityProfile`.

```python
# Plot quality for the first three samples (forward reads)
papa2.plot_quality_profile(fwd_files[:3])
papa2.plot_quality_profile(rev_files[:3])
```

Choose `trunc_len` values where median quality stays above Q30 across all samples.
A typical choice for 2×250 bp V4 16S reads is `trunc_len=(240, 200)`.

---

## 2. Filter and Trim

`filter_and_trim()` applies quality truncation, expected-error filtering,
PhiX removal, and length enforcement — then writes the passing reads to new
FASTQ files.  Both reads in a pair must pass for either to be kept.

```python
filt_dir = "data/filtered"
filt_fwd = [os.path.join(filt_dir, os.path.basename(f)) for f in fwd_files]
filt_rev = [os.path.join(filt_dir, os.path.basename(f)) for f in rev_files]

out = papa2.filter_and_trim(
    fwd_files, filt_fwd,
    rev=rev_files, filt_rev=filt_rev,
    trunc_len=(240, 200),
    max_ee=(2, 2),
    trunc_q=2,
    min_len=50,
    rm_phix=True,
    compress=True,
    multithread=True,
    verbose=True,
)
# out is an (n_files, 2) array: [reads_in, reads_out] per file
print(out)
```

!!! tip "Tuning filter parameters"
    - If too few reads pass, increase `max_ee` or reduce `trunc_len`.
    - If quality drops sharply at read ends, decrease `trunc_len`.
    - Forward and reverse reads can have different truncation lengths.

---

## 3. Dereplicate

`derep_fastq()` collapses identical reads into unique sequences and computes
per-position average quality scores. This step is fast because the C backend
is used when `libpapa2.so` is present.

```python
# Use the filtered files
derepFs = [papa2.derep_fastq(f, verbose=True) for f in filt_fwd]
derepRs = [papa2.derep_fastq(f, verbose=True) for f in filt_rev]

# Each element is a dict:
#   seqs        list[str]   unique sequences, sorted descending by abundance
#   abundances  np.int32    array of read counts
#   quals       np.float64  (n_uniques x max_seqlen) average quality
#   map         np.int32    per-read index into seqs
```

---

## 4. Learn Error Rates

DADA2 builds a parametric error model that describes the probability of observing
each base transition (`A→C`, `A→G`, …) at every quality score. papa2 implements
the same LOESS regression as R's `loessErrfun`.

```python
# Use up to 1e8 bases from the first samples (roughly 10 samples for 1e7 bp each)
errF = papa2.learn_errors(fwd_files, nbases=1e8, verbose=True)
errR = papa2.learn_errors(rev_files, nbases=1e8, verbose=True)

# errF / errR are numpy arrays of shape (16, nqual):
# - 16 rows for the 16 base transitions (A2A, A2C, …, T2T)
# - nqual columns, one per quality-score bin
print("Error matrix shape:", errF.shape)
```

!!! note "Self-consistent learning"
    `learn_errors()` internally runs `dada()` with `self_consist=True`, cycling
    between denoising and error estimation until convergence — exactly as R's
    `learnErrors()` does.

---

## 5. Denoise with DADA

`dada()` applies the core DADA2 algorithm: it partitions the dereplicated reads
into ASVs, evaluating each unique sequence against all more-abundant ones using
the learned error model.

```python
# Denoise forward reads for all samples
dadaFs = papa2.dada(derepFs, err=errF, verbose=False)

# Denoise reverse reads
dadaRs = papa2.dada(derepRs, err=errR, verbose=False)

# Each element of dadaFs / dadaRs is a dict:
#   denoised       dict {seq: abundance}  final ASVs
#   cluster_seqs   list[str]              ASV sequences
#   cluster_abunds np.int32               ASV abundances
#   map            np.int32               unique -> cluster assignment
#   trans          np.int32 (16 x nqual)  observed transition counts

# Report ASV counts per sample
for name, dF in zip(sample_names, dadaFs):
    n_asvs = len(dF["denoised"])
    print(f"  {name}: {n_asvs} ASVs")
```

### Tuning DADA options

All algorithmic parameters from R's `DADA_OPTS` are exposed:

```python
# Tighten the statistical threshold (default is 1e-40)
papa2.set_dada_opt(OMEGA_A=1e-50)

# Or pass options directly to dada()
dadaFs = papa2.dada(derepFs, err=errF, OMEGA_A=1e-50, verbose=False)
```

See `DADA_OPTS` for the full list of tunable parameters.

---

## 6. Merge Paired-End Reads

`merge_pairs()` aligns the denoised forward and reverse ASVs to construct the
full amplicon sequence. Pairs are accepted only when the overlap satisfies both
`min_overlap` and `max_mismatch`.

```python
mergers = []
for name, dF, drF, dR, drR in zip(
        sample_names, dadaFs, derepFs, dadaRs, derepRs):
    m = papa2.merge_pairs(dF, drF, dR, drR,
                          min_overlap=12,
                          max_mismatch=0,
                          verbose=True)
    mergers.append(m)
    n_accept = sum(1 for r in m if r["accept"])
    n_reads  = sum(r["abundance"] for r in m if r["accept"])
    print(f"  {name}: {n_accept} merged pairs, {n_reads} reads")
```

Each element of `mergers` is a list of dicts:

| Key | Type | Description |
|-----|------|-------------|
| `sequence` | `str` | Merged amplicon sequence |
| `abundance` | `int` | Supporting read count |
| `forward` | `int` | Forward cluster index |
| `reverse` | `int` | Reverse cluster index |
| `nmatch` | `int` | Overlap matches |
| `nmismatch` | `int` | Overlap mismatches |
| `nindel` | `int` | Overlap indels |
| `accept` | `bool` | Passed acceptance criteria |

!!! tip "Just concatenate"
    For amplicons where forward and reverse reads do not overlap (e.g. long V1-V3
    16S), use `just_concatenate=True`. The reads are joined with 10 Ns.

---

## 7. Build the Sequence Table

`make_sequence_table()` assembles a (samples × ASVs) abundance matrix from the
per-sample merger lists.

```python
seqtab = papa2.make_sequence_table(mergers)

# seqtab is a dict:
#   table  np.int32 (n_samples x n_asvs)
#   seqs   list[str]

print("Sequence table:", seqtab["table"].shape)
# e.g. (20, 3842) for 20 samples and 3842 unique ASVs

# Inspect the distribution of ASV lengths
import collections
lengths = collections.Counter(len(s) for s in seqtab["seqs"])
print(sorted(lengths.items()))
```

---

## 8. Remove Chimeras

Chimeric sequences are artefacts formed by the spurious joining of two parental
molecules. `remove_bimera_denovo()` identifies and removes bimeras — chimeras with
exactly two parents.

```python
seqtab_nochim = papa2.remove_bimera_denovo(
    seqtab,
    method="consensus",   # flag ASVs flagged as chimeric in most samples
    min_fold=1.5,
    min_abund=2,
    verbose=True
)

# Report read retention
total_before = seqtab["table"].sum()
total_after  = seqtab_nochim["table"].sum()
print(f"Reads retained: {total_after}/{total_before} "
      f"({100*total_after/total_before:.1f}%)")
print(f"ASVs retained: {len(seqtab_nochim['seqs'])}/{len(seqtab['seqs'])}")
```

Available methods:

| Method | Description |
|--------|-------------|
| `"consensus"` | Flag per-sample; remove ASVs chimeric in most samples (default) |
| `"pooled"` | Pool all samples and flag globally |
| `"per-sample"` | Zero only chimeric cells per sample, then drop empty ASVs |

---

## 9. Assign Taxonomy

papa2 includes both a k-mer-based Bayesian classifier (`assign_taxonomy`, matching
R's `assignTaxonomy`) and an exact-match species assigner (`assign_species`).
Download a formatted reference database and point papa2 to it.

### Download reference databases

DADA2-formatted training files are available from
[the DADA2 reference page](https://benjjneb.github.io/dada2/training.html).
Common choices:

| Database | File | Link |
|----------|------|------|
| SILVA v138.1 (genus) | `silva_nr99_v138.1_train_set.fa.gz` | [Download](https://zenodo.org/records/4587955/files/silva_nr99_v138.1_train_set.fa.gz) |
| SILVA v138.1 (species) | `silva_species_assignment_v138.1.fa.gz` | [Download](https://zenodo.org/records/4587955/files/silva_species_assignment_v138.1.fa.gz) |
| GTDB r220 (genus) | `GTDB_bac120_arc53_ssu_r220_genus.fa.gz` | [Download](https://zenodo.org/records/13901193/files/GTDB_bac120_arc53_ssu_r220_genus.fa.gz) |
| GTDB r220 (species) | `GTDB_bac120_arc53_ssu_r220_species.fa.gz` | [Download](https://zenodo.org/records/13901193/files/GTDB_bac120_arc53_ssu_r220_species.fa.gz) |

```bash
# Example: download SILVA v138.1
wget https://zenodo.org/records/4587955/files/silva_nr99_v138.1_train_set.fa.gz
wget https://zenodo.org/records/4587955/files/silva_species_assignment_v138.1.fa.gz
```

### Run the classifier

```python
# Naive Bayesian classifier (genus-level)
taxa = papa2.assign_taxonomy(
    seqtab_nochim["seqs"],
    "silva_nr99_v138.1_train_set.fa.gz",
    min_boot=50,
    verbose=True,
)
# taxa is a DataFrame: rows = ASV sequences, columns = Kingdom..Species

# Exact-match species assignment using the SILVA species database
taxa_with_species = papa2.add_species(
    taxa,
    seqtab_nochim["seqs"],
    "silva_species_assignment_v138.1.fa.gz",
)
```

### Inspect error rates

After learning errors, visualise the fit to check for anomalies:

```python
papa2.plot_errors(errF, output="error_rates_fwd.html")
```

---

## 10. Track Reads Through the Pipeline

It is good practice to track how many reads pass each step. papa2 provides
`track_reads()` and `plot_sankey()` to summarise the entire pipeline —
from raw input through chimera removal and taxonomy — as an interactive
Sankey diagram:

```python
track = papa2.track_reads(
    dereps=derepFs,
    dadas=dadaFs,
    mergers=mergers,
    seqtab=seqtab,
    seqtab_nochim=seqtab_nochim,
    taxa=taxa_with_species,
)
print(track)
# e.g. {'input': 50000, 'denoised': 42000, 'merged': 40000,
#        'tabled': 40000, 'non-chimeric': 38000,
#        'kingdom': 37800, 'phylum': 37200, 'class': 36500,
#        'order': 35800, 'family': 34900, 'genus': 33500, 'species': 28000}

# Save as interactive HTML (or omit output= to display in Jupyter)
papa2.plot_sankey(track, output="read_tracking.html")
```

You can also build a per-sample tracking table for more detail:

```python
import pandas as pd

rows = []
for i, name in enumerate(sample_names):
    n_input    = int(derepFs[i]["abundances"].sum())
    n_denoised = int(dadaFs[i]["cluster_abunds"].sum())
    n_merged   = int(sum(r["abundance"] for r in mergers[i] if r["accept"]))
    rows.append({
        "sample":   name,
        "input":    n_input,
        "denoised": n_denoised,
        "merged":   n_merged,
    })

nochim_table = seqtab_nochim["table"]
for i, row in enumerate(rows):
    row["nonchim"] = int(nochim_table[i].sum())

df = pd.DataFrame(rows)
print(df.to_string(index=False))
```

---

## 11. Export Results

```python
# Write ASVs to a FASTA file
papa2.write_fasta(seqtab_nochim["seqs"], "asvs.fasta")

# Or with a header prefix
papa2.uniquesto_fasta(seqtab_nochim["seqs"], "asvs_with_headers.fasta",
                      prefix="ASV")

# The sequence table and taxonomy can be saved with numpy / pandas:
import numpy as np
np.save("seqtab_nochim.npy", seqtab_nochim["table"])
```

---

## Summary of the Workflow

```
Raw FASTQ files (demultiplexed)
        │
        ▼
filter_and_trim()      # Quality filter, truncate, remove PhiX
        │
        ▼
derep_fastq()          # Collapse identical reads
        │
        ▼
learn_errors()         # Estimate sequencing error rates (LOESS)
        │
        ▼
dada()                 # Infer true sequences (ASVs)
        │
        ▼
merge_pairs()          # Join forward + reverse ASVs
        │
        ▼
make_sequence_table()  # Assemble sample × ASV matrix
        │
        ▼
remove_bimera_denovo() # Remove chimeric ASVs
        │
        ▼
assign_taxonomy()      # Bayesian taxonomic classification
        │
        ▼
add_species()          # Exact-match species assignment
```
