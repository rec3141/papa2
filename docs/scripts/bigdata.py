#!/usr/bin/env python3
"""papa2 big data workflow — streaming per-sample processing.

Full walkthrough: https://rec3141.github.io/papa2/bigdata/
"""

import os, glob
import numpy as np
import papa2

# ── Stage 1: Setup ──────────────────────────────────────────────────────
fastq_dir = "data/run1"
fwd_files = sorted(glob.glob(os.path.join(fastq_dir, "*_R1_001.fastq.gz")))
rev_files = sorted(glob.glob(os.path.join(fastq_dir, "*_R2_001.fastq.gz")))

sample_names = [
    os.path.basename(f).replace("_R1_001.fastq.gz", "") for f in fwd_files
]
print(f"{len(sample_names)} samples found")

# ── Stage 2: Learn Error Rates ──────────────────────────────────────────
errF = papa2.learn_errors(fwd_files, nbases=1e8, randomize=True, verbose=True)
errR = papa2.learn_errors(rev_files, nbases=1e8, randomize=True, verbose=True)

# ── Stage 3: Sample Inference (streaming loop) ─────────────────────────
mergers_dict = {}
for name, fwd, rev in zip(sample_names, fwd_files, rev_files):
    print(f"Processing: {name}")
    drF = papa2.derep_fastq(fwd)
    drR = papa2.derep_fastq(rev)
    ddF = papa2.dada([drF], err=errF, verbose=False)[0]
    ddR = papa2.dada([drR], err=errR, verbose=False)[0]
    merger = papa2.merge_pairs(ddF, drF, ddR, drR, verbose=False)
    mergers_dict[name] = merger

# Build the sequence table
seqtab = papa2.make_sequence_table(mergers_dict)
print(f"Sequence table: {seqtab['table'].shape[0]} samples x {seqtab['table'].shape[1]} ASVs")

# Save intermediate result
np.savez("run1_seqtab.npz", table=seqtab["table"], seqs=seqtab["seqs"])

# ── Stage 5: Remove Chimeras ───────────────────────────────────────────
seqtab_nochim = papa2.remove_bimera_denovo(
    seqtab, method="consensus", verbose=True
)

total_before = seqtab["table"].sum()
total_after = seqtab_nochim["table"].sum()
print(
    f"Reads retained: {total_after}/{total_before} "
    f"({100 * total_after / total_before:.1f}%)"
)

# ── Stage 6: Assign Taxonomy ───────────────────────────────────────────
taxa = papa2.assign_species(
    seqtab_nochim["seqs"],
    "silva_nr99_v138.1_train_set.fa.gz",
    verbose=True,
)

taxa_with_species = papa2.add_species(
    taxa,
    seqtab_nochim["seqs"],
    "silva_species_assignment_v138.1.fa.gz",
)

# ── Stage 7: Track Reads & Export ──────────────────────────────────────
track = papa2.track_reads(
    seqtab=seqtab,
    seqtab_nochim=seqtab_nochim,
    taxa=taxa_with_species,
)
papa2.plot_sankey(track, output="read_tracking.html")

papa2.write_fasta(seqtab_nochim["seqs"], "asvs.fasta")
np.savez(
    "results.npz",
    seqtab=seqtab_nochim["table"],
    seqs=seqtab_nochim["seqs"],
)
print("Done.")
