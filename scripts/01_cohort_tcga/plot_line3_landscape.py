#!/usr/bin/env python3
"""
LINE 3 figure: the immune escape landscape.

Waterfall of patients ranked by mutation silencing ratio (Line 1), coloured by
B2M copy number status (Line 3), with an oncoprint strip underneath showing the
full antigen-presentation machinery.

Reading the figure: if red bars cluster on the left, the two escape routes
reinforce each other. If red and blue interleave evenly, they are independent.
"""
import sys, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.colors import ListedColormap, BoundaryNorm

RES_DIR = sys.argv[1]
FIG_DIR = os.path.join(RES_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

RED = "#C8553D"
BLUE = "#6FB3D0"
DARK = "#5C1A12"
CREAM = "#FBF0EC"
INK = "#16213E"
GREY = "#5A6785"

STRIP = ["B2M", "HLA-A", "HLA-B", "HLA-C", "TAP1", "TAP2", "NLRC5"]

cohorts = []
for c in ["brca", "ov"]:
    f = "%s/%s_line3_landscape.tsv" % (RES_DIR, c)
    if os.path.exists(f):
        cohorts.append((c, pd.read_csv(f, sep="\t")))
    else:
        print("      missing %s - skipping" % f, flush=True)

if not cohorts:
    raise SystemExit("no line3 landscape files; run line3_cna.py first")

print("[1/3] loaded %d cohort(s)" % len(cohorts), flush=True)

n = len(cohorts)
fig = plt.figure(figsize=(15, 5.6 * n))
gs = fig.add_gridspec(2 * n, 1,
                      height_ratios=[4, 1] * n,
                      hspace=0.12 if n == 1 else 0.45)

for k, (cohort, df) in enumerate(cohorts):
    df = df.sort_values("silencing_ratio", ascending=False).reset_index(drop=True)
    x = np.arange(len(df))
    colours = np.where(df.B2M_deleted, RED, BLUE)

    ax = fig.add_subplot(gs[2 * k, 0])
    ax.bar(x, df.silencing_ratio, width=1.0, color=colours, linewidth=0)
    med = df.silencing_ratio.median()
    ax.axhline(med, ls="--", lw=1.1, c=GREY)
    ax.text(len(df) * 0.995, med, "  median %.3f" % med,
            va="bottom", ha="right", fontsize=9, color=GREY)

    pct = 100.0 * df.B2M_deleted.mean()
    ax.set_ylabel("Mutation Silencing Ratio", fontsize=11, fontweight="bold")
    ax.set_title("%s Immune Escape Landscape  (N=%d,  B2M lost in %.1f%%)"
                 % (cohort.upper(), len(df), pct),
                 fontsize=13, fontweight="bold", color=INK, pad=10)
    ax.set_xlim(-0.5, len(df) - 0.5)
    ax.set_ylim(0, max(df.silencing_ratio.max() * 1.08, 0.05))
    ax.set_xticks([])
    ax.legend(handles=[Patch(facecolor=RED, label="B2M Deleted"),
                       Patch(facecolor=BLUE, label="B2M Intact")],
              loc="upper right", frameon=True, fontsize=10)
    ax.grid(axis="y", ls="-", lw=0.7, c="#E6EBF2")
    ax.set_axisbelow(True)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)

    axs = fig.add_subplot(gs[2 * k + 1, 0])
    rows = [g for g in STRIP if "%s_cna" % g in df.columns]
    M = np.vstack([(df["%s_cna" % g].values < 0).astype(int) for g in rows])
    axs.imshow(M, aspect="auto", interpolation="nearest",
               cmap=ListedColormap([CREAM, DARK]),
               norm=BoundaryNorm([-0.5, 0.5, 1.5], 2))
    axs.set_yticks(np.arange(len(rows)))
    axs.set_yticklabels(rows, fontsize=9)
    axs.set_xlabel("Individual Patients  (ranked by silencing ratio, high to low)",
                   fontsize=11, fontweight="bold")
    axs.set_xticks(np.linspace(0, len(df) - 1, 6).astype(int))
    axs.set_xticklabels(np.linspace(0, len(df) - 1, 6).astype(int), fontsize=9)
    axs.tick_params(axis="y", length=0)
    for s in axs.spines.values():
        s.set_edgecolor("#C4CCD8")
    axs.text(1.005, 0.5, "copy\nnumber\nloss", transform=axs.transAxes,
             va="center", fontsize=8.5, color=GREY)

fig.suptitle("Line 3  |  Genomic Disruption of Antigen Presentation (Hardware Kill-Switch)",
             fontsize=15, fontweight="bold", color=INK, x=0.012, ha="left",
             y=0.995 if n == 1 else 0.997)
fig.tight_layout(rect=[0, 0, 0.985, 0.97 if n == 1 else 0.985])

out = "%s/line3_escape_landscape.png" % FIG_DIR
fig.savefig(out, dpi=200, facecolor="white")
print("[2/3] wrote %s" % out, flush=True)
print("[3/3] done", flush=True)
