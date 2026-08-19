#!/usr/bin/env python3
"""Aggregate bench2 logs into LaTeX tables (steps x datasets x engines)."""
import glob, json, os, re, statistics as st, sys

RUNS = sys.argv[1] if len(sys.argv) > 1 else "bench2_runs"
DATASETS = [("mock", "Mock (2$\\times$250)"), ("gut", "Gut (2$\\times$300)"),
            ("mosaic", "MOSAiC (2$\\times$151)")]
STEPS = [("filter", "Filter \\& trim"), ("derep", "Dereplicate"),
         ("learnErrorsF", "Learn errors (F)"), ("learnErrorsR", "Learn errors (R)"),
         ("dadaF", "Denoise (F)"), ("dadaR", "Denoise (R)"),
         ("merge", "Merge pairs"), ("seqtab", "Sequence table"),
         ("chimera", "Chimera removal"), ("taxonomy", "Taxonomy (SILVA)")]

def collect(ds, engine):
    """{step: [times]}, [resmon dicts]"""
    times = {}
    res = []
    for rep in (1, 2, 3):
        log = f"{RUNS}/{ds}_{engine}_rep{rep}.log"
        if not os.path.exists(log):
            continue
        for line in open(log):
            m = re.match(r"TIME (\S+) ([\d.]+)", line)
            if m:
                times.setdefault(m.group(1), []).append(float(m.group(2)))
        rj = f"{RUNS}/{ds}_{engine}_rep{rep}.resmon.json"
        if os.path.exists(rj):
            res.append(json.load(open(rj)))
    return times, res

def fmt(vals):
    if not vals:
        return "--"
    m = st.mean(vals)
    s = st.stdev(vals) if len(vals) > 1 else 0.0
    if m >= 100:
        return f"{m:.0f}$\\pm${s:.0f}"
    if m >= 10:
        return f"{m:.1f}$\\pm${s:.1f}"
    return f"{m:.2f}$\\pm${s:.2f}"

data = {ds: {e: collect(ds, e) for e in ("R", "py")} for ds, _ in DATASETS}

# ---- Table: per-step wall times ----
lines = []
lines.append("\\begin{tabular}{l" + "rrr" * len(DATASETS) + "}")
lines.append("\\toprule")
hdr = " & " + " & ".join(f"\\multicolumn{{3}}{{c}}{{{label}}}" for _, label in DATASETS) + " \\\\"
lines.append(hdr)
sub = "Step" + " & R & papa2 & $\\times$" * len(DATASETS) + " \\\\"
cml = " ".join(f"\\cmidrule(lr){{{2+3*i}-{4+3*i}}}" for i in range(len(DATASETS)))
lines.append(cml)
lines.append(sub)
lines.append("\\midrule")
for key, label in STEPS:
    row = [label]
    for ds, _ in DATASETS:
        tR = data[ds]["R"][0].get(key, [])
        tP = data[ds]["py"][0].get(key, [])
        row.append(fmt(tR)); row.append(fmt(tP))
        if tR and tP and st.mean(tP) > 0:
            row.append(f"{st.mean(tR)/st.mean(tP):.1f}")
        else:
            row.append("--")
    lines.append(" & ".join(row) + " \\\\")
lines.append("\\midrule")
row = ["\\textbf{Total (wall)}"]
tot_line2 = ["Peak memory (GB)"]
tot_line3 = ["CPU time (core-s)"]
for ds, _ in DATASETS:
    wR = [r["wall_s"] for r in data[ds]["R"][1]]
    wP = [r["wall_s"] for r in data[ds]["py"][1]]
    row.append("\\textbf{" + fmt(wR) + "}"); row.append("\\textbf{" + fmt(wP) + "}")
    row.append("\\textbf{" + (f"{st.mean(wR)/st.mean(wP):.1f}" if wR and wP else "--") + "}")
    mR = [r["peak_rss_gb"] for r in data[ds]["R"][1]]
    mP = [r["peak_rss_gb"] for r in data[ds]["py"][1]]
    tot_line2 += [fmt(mR), fmt(mP), f"{st.mean(mR)/st.mean(mP):.1f}" if mR and mP else "--"]
    cR = [r["cpu_user_s"] + r["cpu_sys_s"] for r in data[ds]["R"][1]]
    cP = [r["cpu_user_s"] + r["cpu_sys_s"] for r in data[ds]["py"][1]]
    tot_line3 += [fmt(cR), fmt(cP), f"{st.mean(cR)/st.mean(cP):.1f}" if cR and cP else "--"]
lines.append(" & ".join(row) + " \\\\")
lines.append(" & ".join(tot_line2) + " \\\\")
lines.append(" & ".join(tot_line3) + " \\\\")
lines.append("\\bottomrule")
lines.append("\\end{tabular}")
open(os.path.join(RUNS, "table_steps.tex"), "w").write("\n".join(lines) + "\n")
print("wrote table_steps.tex")
for ds, _ in DATASETS:
    for e in ("R", "py"):
        n = len(data[ds][e][1])
        print(f"  {ds} {e}: {n} replicates")
