#!/usr/bin/env python3
"""
LINE 1 figure 2: the genes the tumour must express versus the genes it hides.

Four horizontal bar panels:
  top row    - BRCA low silencing  |  BRCA high silencing
  bottom row - OV   low silencing  |  OV   high silencing

Only genes with at least MIN_MUT mutations are shown, so every bar rests on
a defensible number of observations.
"""
import sys, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RES_DIR = sys.argv[1]
MIN_B = int(sys.argv[2]) if len(sys.argv) > 2 else 10
MIN_O = int(sys.argv[3]) if len(sys.argv) > 3 else 5
TOP_N = 15

FIG_DIR = os.path.join(RES_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

TEAL = "#1C7293"
AMBER = "#D98E28"
INK = "#16213E"

print("[1/3] loading ...", flush=True)
brca = pd.read_csv("%s/brca_all_genes.csv" % RES_DIR)
ov = pd.read_csv("%s/ov_all_genes.csv" % RES_DIR)
brca = brca[brca.Total_Mutations >= MIN_B]
ov = ov[ov.Total_Mutations >= MIN_O]
print("      BRCA genes: %d   OV genes: %d" % (len(brca), len(ov)), flush=True)


def pick(df, high, n):
    return (df.sort_values(["Silenced_Percent", "Total_Mutations"],
                           ascending=[not high, False])
              .head(n).sort_values("Silenced_Percent"))


def panel(ax, df, colour, title, subtitle):
    y = np.arange(len(df))
    ax.barh(y, df.Silenced_Percent, color=colour, height=0.68,
            edgecolor="white", linewidth=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels(df.Hugo_Symbol, fontsize=10)
    ax.set_xlim(0, 118)
    ax.set_xticks([0, 25, 50, 75, 100])
    for yi, (pct, n) in enumerate(zip(df.Silenced_Percent, df.Total_Mutations)):
        ax.text(pct + 2.5, yi, "%.0f%%  (n=%d)" % (pct, n),
                va="center", fontsize=8.5, color=INK)
    ax.set_title(title, fontsize=12, fontweight="bold", color=INK,
                 loc="left", pad=8)
    ax.text(0, 1.015, subtitle, transform=ax.transAxes,
            fontsize=9.5, color="#5A6785")
    ax.set_xlabel("Silenced mutations (%)", fontsize=10.5)
    ax.grid(axis="x", ls="-", lw=0.7, c="#E6EBF2")
    ax.set_axisbelow(True)
    for s in ["top", "right", "left"]:
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="y", length=0)


print("[2/3] drawing ...", flush=True)
fig, axes = plt.subplots(2, 2, figsize=(15, 11))

panel(axes[0, 0], pick(brca, False, TOP_N), TEAL,
      "BREAST  -  LOW silencing",
      "drivers the tumour is forced to keep expressing (n >= %d mutations)" % MIN_B)
panel(axes[0, 1], pick(brca, True, TOP_N), AMBER,
      "BREAST  -  HIGH silencing",
      "passengers the tumour hides from the immune system")
panel(axes[1, 0], pick(ov, False, TOP_N), TEAL,
      "OVARIAN  -  LOW silencing",
      "drivers the tumour is forced to keep expressing (n >= %d mutations)" % MIN_O)
panel(axes[1, 1], pick(ov, True, TOP_N), AMBER,
      "OVARIAN  -  HIGH silencing",
      "passengers the tumour hides from the immune system")

fig.suptitle("Line 1  |  Stealth Mode:  what the tumour must show and what it hides",
             fontsize=15, fontweight="bold", color=INK, x=0.012, ha="left", y=0.985)
fig.tight_layout(rect=[0, 0, 1, 0.955], h_pad=3.5, w_pad=3.0)

out = "%s/line1_top_genes.png" % FIG_DIR
fig.savefig(out, dpi=200, facecolor="white")
print("[3/3] wrote %s" % out, flush=True)
