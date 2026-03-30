#!/usr/bin/env python3
"""papa2 quickstart — minimal forward-read denoising example.

Full walkthrough: https://rec3141.github.io/papa2/quickstart/
"""

import papa2

# The bundled test files (forward reads from two 16S amplicon samples)
fwd_files = ["tests/data/sam1F.fastq.gz", "tests/data/sam2F.fastq.gz"]

# --- 1. Dereplicate ---
dereps = [papa2.derep_fastq(f, verbose=True) for f in fwd_files]

# --- 2. Learn error rates ---
err = papa2.learn_errors(fwd_files, verbose=True)

# --- 3. Denoise ---
dadas = papa2.dada(dereps, err=err, verbose=True)

# --- 4. Inspect results ---
for i, dd in enumerate(dadas):
    print(
        f"Sample {i+1}: {len(dd['denoised'])} ASVs, "
        f"{sum(dd['denoised'].values())} reads"
    )

# --- 5. Build sequence table ---
seqtab = papa2.make_sequence_table(dadas)
print(f"Sequence table: {seqtab.shape[0]} samples x {seqtab.shape[1]} ASVs")

# --- 6. Visualise read tracking (Sankey diagram) ---
track = papa2.track_reads(dereps=dereps, dadas=dadas, seqtab=seqtab)
papa2.plot_sankey(track, output="read_tracking.html")
print("Sankey diagram saved to read_tracking.html")
