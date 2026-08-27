#!/usr/bin/env python3
"""
LINE 2 figure: MHC-I complex expression by mutation silencing strategy.

Grouped bar chart with SEM error bars, one panel per cohort, matching the
classic Line 2 presentation figure.
"""
import sys, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RES_DIR = sys.argv[1]
FIG_DIR = os.path.join(RES_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

MHC = ["B2M", "HLA-A", "HLA-B", "HLA-C"]
BLUE = "#4C93A8"
RED = "#C8553D"
INK = "#16213E"
GREY = "#5A6785"

print("[1/3] loading ...", flush=True)
frames = {}
for c in ["brca", "ov"]:
    f = "%s/%s_line2_stats.csv" % (RES_DIR, c)
    if os.path.exists(f):
        frames[c] = pd.read_csv(f)
    else:
        print("      missing %s - skipping" % f, flush=True)

if not frames:
    raise SystemExit("no line2 stats found; run line2_mhc.py first")

print("[2/3] drawing ...", flush=True)
fig, axes = plt.subplots(1, len(frames), figsize=(7.2 * len(frames), 6.4),
                         squeeze=False)
axes = axes[0]

for ax, (cohort, df) in zip(axes, frames.items()):
    d = df[df.Gene.isin(MHC)].set_index("Gene").loc[
        [g for g in MHC if g in df.Gene.values]]
    x = np.arange(len(d))
    w = 0.38

    ax.bar(x - w / 2, d.mean_high_log2, w, yerr=d.sem_high, capsize=4,
           color=BLUE, edgecolor="black", linewidth=0.8,
           error_kw=dict(lw=1.2, ecolor="black"), label="High Silencing")
    ax.bar(x + w / 2, d.mean_low_log2, w, yerr=d.sem_low, capsize=4,
           color=RED, edgecolor="black", linewidth=0.8,
           error_kw=dict(lw=1.2, ecolor="black"), label="Low Silencing")

    top = max((d.mean_high_log2 + d.sem_high).max(),
              (d.mean_low_log2 + d.sem_low).max())
    for i, (_, r) in enumerate(d.iterrows()):
        star = "n.s." if r.p_bonferroni >= 0.05 else (
            "***" if r.p_bonferroni < 0.001 else
            "**" if r.p_bonferroni < 0.01 else "*")
        ax.text(i, top * 1.035, star, ha="center", fontsize=10,
                color=GREY if star == "n.s." else INK,
                fontweight="normal" if star == "n.s." else "bold")
        ax.text(i, top * 1.085, "d=%.2f" % r.cohens_d, ha="center",
                fontsize=8.5, color=GREY)

    ax.set_xticks(x)
    ax.set_xticklabels(d.index, fontsize=11)
    ax.set_ylim(0, top * 1.16)
    ax.set_xlabel("Gene", fontsize=11.5, fontweight="bold")
    ax.set_ylabel("Expression, log2(RSEM + 1)", fontsize=11.5, fontweight="bold")
    n_hi = int(d.n_high.iloc[0])
    n_lo = int(d.n_low.iloc[0])
    ax.set_title("%s   (High n=%d,  Low n=%d)" % (cohort.upper(), n_hi, n_lo),
                 fontsize=13, fontweight="bold", color=INK, pad=12)
    ax.legend(title="Patient Strategy", frameon=True, fontsize=10,
              title_fontsize=10, loc="lower left")
    ax.grid(axis="y", ls="-", lw=0.7, c="#E6EBF2")
    ax.set_axisbelow(True)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)

fig.suptitle("Line 2  |  MHC-I Complex Expression by Mutation Silencing Strategy",
             fontsize=15, fontweight="bold", color=INK, x=0.012, ha="left", y=0.985)
fig.tight_layout(rect=[0, 0, 1, 0.945], w_pad=3.0)

out = "%s/line2_mhc_expression.png" % FIG_DIR
fig.savefig(out, dpi=200, facecolor="white")
print("[3/3] wrote %s" % out, flush=True)
