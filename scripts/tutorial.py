#!/usr/bin/env python3
"""papa2 tutorial — paired-end amplicon denoising pipeline.

Full walkthrough: https://rec3141.github.io/papa2/tutorial/
"""

import os, glob, collections
import numpy as np
import papa2

# ── 0. Setup ────────────────────────────────────────────────────────────────
fastq_dir = "data/fastq"

fwd_files = sorted(glob.glob(os.path.join(fastq_dir, "*_R1_001.fastq.gz")))
rev_files = sorted(glob.glob(os.path.join(fastq_dir, "*_R2_001.fastq.gz")))

sample_names = [
    os.path.basename(f).replace("_R1_001.fastq.gz", "") for f in fwd_files
]
print(f"{len(sample_names)} samples: {sample_names[:5]} ...")

# ── 1. Inspect Read Quality ────────────────────────────────────────────────
papa2.plot_quality_profile(fwd_files[:3])
papa2.plot_quality_profile(rev_files[:3])

# ── 2. Filter and Trim ─────────────────────────────────────────────────────
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
print(out)

# ── 3. Dereplicate ─────────────────────────────────────────────────────────
derepFs = [papa2.derep_fastq(f, verbose=True) for f in filt_fwd]
derepRs = [papa2.derep_fastq(f, verbose=True) for f in filt_rev]

# ── 4. Learn Error Rates ──────────────────────────────────────────────────
errF = papa2.learn_errors(filt_fwd, nbases=1e8, verbose=True)
errR = papa2.learn_errors(filt_rev, nbases=1e8, verbose=True)
print("Error matrix shape:", errF.shape)

# Inspect error rate fit
papa2.plot_errors(errF, output="error_rates_fwd.html")

# ── 5. Denoise ─────────────────────────────────────────────────────────────
dadaFs = papa2.dada(derepFs, err=errF, verbose=False)
dadaRs = papa2.dada(derepRs, err=errR, verbose=False)

for name, dF in zip(sample_names, dadaFs):
    print(f"  {name}: {len(dF['denoised'])} ASVs")

# ── 6. Merge Paired-End Reads ─────────────────────────────────────────────
mergers = []
for name, dF, drF, dR, drR in zip(
    sample_names, dadaFs, derepFs, dadaRs, derepRs
):
    m = papa2.merge_pairs(dF, drF, dR, drR, min_overlap=12, max_mismatch=0, verbose=True)
    mergers.append(m)
    n_accept = sum(1 for r in m if r["accept"])
    n_reads = sum(r["abundance"] for r in m if r["accept"])
    print(f"  {name}: {n_accept} merged pairs, {n_reads} reads")

# ── 7. Build Sequence Table ──────────────────────────────────────────────
seqtab = papa2.make_sequence_table(mergers)
print("Sequence table:", seqtab["table"].shape)

lengths = collections.Counter(len(s) for s in seqtab["seqs"])
print("ASV length distribution:", sorted(lengths.items()))

# ── 8. Remove Chimeras ───────────────────────────────────────────────────
seqtab_nochim = papa2.remove_bimera_denovo(
    seqtab, method="consensus", min_fold=1.5, min_abund=2, verbose=True
)

total_before = seqtab["table"].sum()
total_after = seqtab_nochim["table"].sum()
print(
    f"Reads retained: {total_after}/{total_before} "
    f"({100 * total_after / total_before:.1f}%)"
)
print(f"ASVs retained: {len(seqtab_nochim['seqs'])}/{len(seqtab['seqs'])}")

# ── 9. Assign Taxonomy ──────────────────────────────────────────────────
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

# ── 10. Track Reads & Sankey ─────────────────────────────────────────────
track = papa2.track_reads(
    dereps=derepFs,
    dadas=dadaFs,
    mergers=mergers,
    seqtab=seqtab,
    seqtab_nochim=seqtab_nochim,
    taxa=taxa_with_species,
)
print(track)
papa2.plot_sankey(track, output="read_tracking.html")
print("Sankey diagram saved to read_tracking.html")

# ── 11. Export ───────────────────────────────────────────────────────────
papa2.write_fasta(seqtab_nochim["seqs"], "asvs.fasta")
np.save("seqtab_nochim.npy", seqtab_nochim["table"])
print("Done.")
