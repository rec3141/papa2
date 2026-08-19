#!/usr/bin/env python3
"""Per-step benchmark: papa2, default parallelism. Same args as bench2_r.R."""
import os, sys, time, json
import numpy as np
import papa2

raw, out = sys.argv[1], sys.argv[2]
tf, tr = int(sys.argv[3]), int(sys.argv[4])
eeF, eeR = float(sys.argv[5]), float(sys.argv[6])
silva = sys.argv[7]
filt = os.path.join(out, "filtered")
os.makedirs(filt, exist_ok=True)

accs = open(os.path.join(raw, "accs.txt")).read().split()
fnFs = [f"{raw}/{a}_1.fastq.gz" for a in accs]
fnRs = [f"{raw}/{a}_2.fastq.gz" for a in accs]
filtFs = [f"{filt}/{a}_F.fastq.gz" for a in accs]
filtRs = [f"{filt}/{a}_R.fastq.gz" for a in accs]

def tt(name, fn):
    t0 = time.perf_counter(); res = fn()
    print(f"TIME {name} {time.perf_counter()-t0:.2f}", flush=True)
    return res

ft = tt("filter", lambda: papa2.filter_and_trim(
    fnFs, filtFs, rev=fnRs, filt_rev=filtRs,
    trunc_len=(tf, tr), max_ee=(eeF, eeR), trunc_q=2, rm_phix=True,
    compress=True, multithread=True))
np.savetxt(os.path.join(out, "filter_result.tsv"), ft, fmt="%d")
keep = ft[:, 1] > 0
accs = [a for a, k in zip(accs, keep) if k]
filtFs = [f for f, k in zip(filtFs, keep) if k]
filtRs = [f for f, k in zip(filtRs, keep) if k]
open(os.path.join(out, "kept.txt"), "w").write("\n".join(accs) + "\n")

dereps = tt("derep", lambda: (papa2.derep_fastq(filtFs), papa2.derep_fastq(filtRs)))
derepFs, derepRs = dereps

errF = tt("learnErrorsF", lambda: papa2.learn_errors(filtFs, verbose=False))
errR = tt("learnErrorsR", lambda: papa2.learn_errors(filtRs, verbose=False))
np.savetxt(os.path.join(out, "err_fwd.tsv"), errF)

dadaFs = tt("dadaF", lambda: papa2.dada(derepFs, err=errF, verbose=False))
dadaRs = tt("dadaR", lambda: papa2.dada(derepRs, err=errR, verbose=False))
if isinstance(dadaFs, dict):
    dadaFs, dadaRs = [dadaFs], [dadaRs]

mergers = tt("merge", lambda: papa2.merge_pairs(dadaFs, derepFs, dadaRs, derepRs))
seqtab = tt("seqtab", lambda: papa2.make_sequence_table({a: m for a, m in zip(accs, mergers)}))
nochim = tt("chimera", lambda: papa2.remove_bimera_denovo(seqtab))
print("dims:", seqtab["table"].shape, "->", nochim["table"].shape)
json.dump({"seqs": seqtab["seqs"], "table": seqtab["table"].tolist()},
          open(os.path.join(out, "seqtab.json"), "w"))
json.dump({"seqs": nochim["seqs"], "table": nochim["table"].tolist()},
          open(os.path.join(out, "seqtab_nochim.json"), "w"))
for a, d in zip(accs, dadaFs):
    open(os.path.join(out, f"dada_{a}_F.txt"), "w").write("\n".join(d["cluster_seqs"]) + "\n")

taxa = tt("taxonomy", lambda: papa2.assign_taxonomy(
    nochim["seqs"], silva, seed=100,
    tax_levels=("Kingdom","Phylum","Class","Order","Family","Genus")))
taxa.to_csv(os.path.join(out, "taxonomy.csv"))
