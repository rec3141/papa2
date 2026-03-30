# Big Data Workflow

[:material-download: Download bigdata.py](scripts/bigdata.py){ .md-button }

This workflow is designed for large-scale amplicon studies — hundreds of samples,
HiSeq-scale data (100M+ reads), or multi-run projects. It mirrors the
[DADA2 Big Data tutorial](https://benjjneb.github.io/dada2/bigdata.html) but
uses papa2's Python API.

**Key differences from the [standard tutorial](tutorial.md):**

- Samples are processed **one at a time** in a loop, keeping memory low
- Error rates are learned from a **random subset** of bases
- The workflow is split into **independent stages** that can be run as separate jobs
- Multiple sequencing runs can be **merged** before chimera removal

---

## Strategy

DADA2's sample inference is independent per sample once the error model is
learned. This means the per-sample dereplicate → denoise loop can be streamed
with only one sample in memory at a time, keeping RAM usage around 4–8 GB even
for a full HiSeq lane (~750 samples).

Scaling behaviour:

- **Linear** in number of samples
- **Quadratic** in per-sample depth (fewer samples at higher depth = longer)
- Parallelisable across samples (set `DADA2_WORKERS` or farm out to multiple nodes)

---

## Stage 1: Filter and Trim

Quality-filter and truncate reads with `filter_and_trim()`.  For paired-end
data, both reads in a pair must pass for either to be kept.

```python
import os, glob
import papa2

# Point to your directory of demultiplexed FASTQs
fastq_dir = "data/run1"
fwd_files = sorted(glob.glob(os.path.join(fastq_dir, "*_R1_001.fastq.gz")))
rev_files = sorted(glob.glob(os.path.join(fastq_dir, "*_R2_001.fastq.gz")))

sample_names = [
    os.path.basename(f).replace("_R1_001.fastq.gz", "")
    for f in fwd_files
]
print(f"{len(sample_names)} samples found")

# Filter and trim
filt_dir = "data/run1/filtered"
filt_fwd = [os.path.join(filt_dir, os.path.basename(f)) for f in fwd_files]
filt_rev = [os.path.join(filt_dir, os.path.basename(f)) for f in rev_files]

out = papa2.filter_and_trim(
    fwd_files, filt_fwd,
    rev=rev_files, filt_rev=filt_rev,
    trunc_len=(240, 200),
    max_ee=(2, 2),
    trunc_q=2,
    rm_phix=True,
    multithread=True,
    verbose=True,
)
print(out)
```

---

## Stage 2: Learn Error Rates

Learn the error model from a subset of the data. 100 million bases is usually
sufficient — papa2 will draw from the first N files until the target is reached.

```python
errF = papa2.learn_errors(filt_fwd, nbases=1e8, randomize=True, verbose=True)
errR = papa2.learn_errors(filt_rev, nbases=1e8, randomize=True, verbose=True)
```

---

## Stage 3: Sample Inference (streaming loop)

Process each sample individually — only one dereplicated sample is in memory at
a time. Results are collected into a dict for the sequence table.

### Single-end

```python
import papa2

dds = {}
for name, fwd in zip(sample_names, fwd_files):
    print(f"Processing: {name}")
    derep = papa2.derep_fastq(fwd)
    dd = papa2.dada([derep], err=errF, verbose=False)[0]
    dds[name] = dd["denoised"]

# Build the sequence table
seqtab = papa2.make_sequence_table(dds)
print(f"Sequence table: {seqtab.shape[0]} samples x {seqtab.shape[1]} ASVs")
```

### Paired-end

```python
import papa2

mergers_dict = {}
for name, fwd, rev in zip(sample_names, fwd_files, rev_files):
    print(f"Processing: {name}")
    drF = papa2.derep_fastq(fwd)
    drR = papa2.derep_fastq(rev)
    ddF = papa2.dada([drF], err=errF, verbose=False)[0]
    ddR = papa2.dada([drR], err=errR, verbose=False)[0]
    merger = papa2.merge_pairs(ddF, drF, ddR, drR, verbose=False)
    mergers_dict[name] = merger

# Build the sequence table from the merger results
seqtab = papa2.make_sequence_table(mergers_dict)
print(f"Sequence table: {seqtab.shape[0]} samples x {seqtab.shape[1]} ASVs")
```

!!! tip "Saving intermediate results"
    For very large studies, save the sequence table to disk between stages so
    each stage can be an independent job:

    ```python
    import numpy as np
    np.savez("run1_seqtab.npz", table=seqtab["table"], seqs=seqtab["seqs"])
    ```

---

## Stage 4: Merge Runs (if applicable)

If your study spans multiple sequencing runs, learn errors and infer sequences
per-run, then merge the sequence tables:

```python
import numpy as np
import papa2

# Load per-run tables
run1 = np.load("run1_seqtab.npz", allow_pickle=True)
run2 = np.load("run2_seqtab.npz", allow_pickle=True)

# Merge — sequences are matched by identity across runs
seqtab = papa2.collapse_no_mismatch(
    # Combine into a single table keyed by sequence
    papa2.make_sequence_table({
        **{f"run1_{i}": dict(zip(run1["seqs"], run1["table"][i]))
           for i in range(run1["table"].shape[0])},
        **{f"run2_{i}": dict(zip(run2["seqs"], run2["table"][i]))
           for i in range(run2["table"].shape[0])},
    })
)
```

---

## Stage 5: Remove Chimeras

Use `method="consensus"` for large datasets — it flags chimeras per-sample and
removes ASVs that are chimeric in the majority of samples:

```python
seqtab_nochim = papa2.remove_bimera_denovo(
    seqtab,
    method="consensus",
    verbose=True,
)

total_before = seqtab["table"].sum()
total_after = seqtab_nochim["table"].sum()
print(f"Reads retained: {total_after}/{total_before} "
      f"({100 * total_after / total_before:.1f}%)")
```

---

## Stage 6: Assign Taxonomy

```python
taxa = papa2.assign_taxonomy(
    seqtab_nochim["seqs"],
    "silva_nr99_v138.1_train_set.fa.gz",
    min_boot=50,
    verbose=True,
)

taxa_with_species = papa2.add_species(
    taxa,
    seqtab_nochim["seqs"],
    "silva_species_assignment_v138.1.fa.gz",
)
```

See the [Tutorial](tutorial.md#download-reference-databases) for database
download links.

---

## Stage 7: Track Reads and Export

```python
# Sankey diagram — tracks reads through every stage including taxonomy
track = papa2.track_reads(
    seqtab=seqtab,
    seqtab_nochim=seqtab_nochim,
    taxa=taxa_with_species,
)
papa2.plot_sankey(track, output="read_tracking.html")

# Export ASVs
papa2.write_fasta(seqtab_nochim["seqs"], "asvs.fasta")

# Save final tables
import numpy as np
np.savez(
    "results.npz",
    seqtab=seqtab_nochim["table"],
    seqs=seqtab_nochim["seqs"],
)
```

---

## Performance Tips

- Set `DADA2_WORKERS=N` to parallelise the `dada()` calls across CPU cores.
  Each worker processes one sample, so set this to the number of cores available.
- Set `OMP_NUM_THREADS=1` to avoid contention between Python-level and
  OpenMP-level parallelism.
- For cluster/cloud: each run's Stage 2–3 is independent and can be submitted
  as a separate job. Only Stage 4–6 requires all runs to be complete.
- Memory scales with the largest single sample, not the total dataset size.
  Typical usage is 4–8 GB for HiSeq-depth samples.
