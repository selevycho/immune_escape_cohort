#!/usr/bin/env python3
"""
LINE 1 figure 1: convergence of transcriptomic silencing between
breast and ovarian cancer.

Left  : scatter of Silenced_Percent, BRCA vs OV, one dot per gene
Right : threshold sensitivity, showing the ratio holds at every cutoff
"""
import sys, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

RES_DIR = sys.argv[1]
FIG_DIR = os.path.join(RES_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

TEAL = "#1C7293"
AMBER = "#D98E28"
INK = "#16213E"
GREY = "#9AA5B8"

print("[1/3] loading tables ...", flush=True)
brca = pd.read_csv("%s/brca_all_genes.csv" % RES_DIR)
ov = pd.read_csv("%s/ov_all_genes.csv" % RES_DIR)
m = brca.merge(ov, on="Hugo_Symbol", suffixes=("_brca", "_ov"))
print("      genes in both cohorts: %d" % len(m), flush=True)

sw_b = pd.read_csv("%s/brca_threshold_sweep.csv" % RES_DIR)
sw_o = pd.read_csv("%s/ov_threshold_sweep.csv" % RES_DIR)

shared_low = set(pd.read_csv(
    "%s/convergence_shared_LOW_silencing.csv" % RES_DIR).Hugo_Symbol)
shared_high = set(pd.read_csv(
    "%s/convergence_shared_HIGH_silencing.csv" % RES_DIR).Hugo_Symbol)

rho, p_rho = stats.spearmanr(m.Silenced_Percent_brca, m.Silenced_Percent_ov)

print("[2/3] drawing ...", flush=True)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6.2),
                               gridspec_kw={"width_ratios": [1.45, 1]})

# ---- panel A : scatter ----
rng = np.random.default_rng(7)
jx = rng.normal(0, 1.1, len(m))
jy = rng.normal(0, 1.1, len(m))

is_low = m.Hugo_Symbol.isin(shared_low)
is_high = m.Hugo_Symbol.isin(shared_high)
other = ~(is_low | is_high)

ax1.scatter(m.Silenced_Percent_brca[other] + jx[other.values],
            m.Silenced_Percent_ov[other] + jy[other.values],
            s=18, c=GREY, alpha=0.45, linewidths=0, label="other genes")
ax1.scatter(m.Silenced_Percent_brca[is_low] + jx[is_low.values],
            m.Silenced_Percent_ov[is_low] + jy[is_low.values],
            s=70, c=TEAL, edgecolors="white", linewidths=0.8,
            label="shared LOW silencing", zorder=3)
ax1.scatter(m.Silenced_Percent_brca[is_high] + jx[is_high.values],
            m.Silenced_Percent_ov[is_high] + jy[is_high.values],
            s=70, c=AMBER, edgecolors="white", linewidths=0.8,
            label="shared HIGH silencing", zorder=3)

for _, r in m[is_low | is_high].iterrows():
    ax1.annotate(r.Hugo_Symbol,
                 (r.Silenced_Percent_brca, r.Silenced_Percent_ov),
                 fontsize=8.5, color=INK,
                 xytext=(6, 5), textcoords="offset points")

ax1.plot([-3, 103], [-3, 103], ls="--", lw=1, c="#C4CCD8", zorder=0)
ax1.set_xlim(-5, 108)
ax1.set_ylim(-5, 108)
ax1.set_xlabel("Silenced mutations in BREAST (%)", fontsize=11.5)
ax1.set_ylabel("Silenced mutations in OVARIAN (%)", fontsize=11.5)
ax1.set_title("A   Silencing converges between the two cancers",
              fontsize=13, fontweight="bold", color=INK, loc="left", pad=12)
ax1.text(0.03, 0.95,
         "Spearman rho = %.3f\np = %.2g\nn = %d genes" % (rho, p_rho, len(m)),
         transform=ax1.transAxes, va="top", fontsize=10.5, color=INK,
         bbox=dict(boxstyle="round,pad=0.45", fc="#EEF3F7", ec="none"))
ax1.legend(loc="lower right", frameon=False, fontsize=9.5)
for s in ["top", "right"]:
    ax1.spines[s].set_visible(False)

# ---- panel B : threshold sweep ----
ax2.plot(sw_b.RSEM_threshold, sw_b.silenced_pct_of_all_mutations,
         "o-", c=TEAL, lw=2.2, ms=8, label="BRCA")
ax2.plot(sw_o.RSEM_threshold, sw_o.silenced_pct_of_all_mutations,
         "s-", c=AMBER, lw=2.2, ms=8, label="OV")
for _, r in sw_b.iterrows():
    ax2.annotate("%.1f%%" % r.silenced_pct_of_all_mutations,
                 (r.RSEM_threshold, r.silenced_pct_of_all_mutations),
                 fontsize=9, color=TEAL, xytext=(0, 9),
                 textcoords="offset points", ha="center")
for _, r in sw_o.iterrows():
    ax2.annotate("%.1f%%" % r.silenced_pct_of_all_mutations,
                 (r.RSEM_threshold, r.silenced_pct_of_all_mutations),
                 fontsize=9, color=AMBER, xytext=(0, -16),
                 textcoords="offset points", ha="center")
ax2.axvline(5.0, ls=":", lw=1.4, c="#C4CCD8")
ax2.text(5.4, ax2.get_ylim()[0] + 1, "chosen\ncutoff",
         fontsize=9, color=GREY)
ax2.set_xscale("log")
ax2.set_xticks(sw_b.RSEM_threshold)
ax2.set_xticklabels([str(int(t)) for t in sw_b.RSEM_threshold])
ax2.set_xlabel("RSEM threshold for \"not expressed\"", fontsize=11.5)
ax2.set_ylabel("Silenced fraction of all mutations (%)", fontsize=11.5)
ax2.set_title("B   The gap between cancers is threshold-independent",
              fontsize=13, fontweight="bold", color=INK, loc="left", pad=12)
ax2.legend(frameon=False, fontsize=10.5)
ax2.grid(axis="y", ls="-", lw=0.7, c="#E6EBF2")
ax2.set_axisbelow(True)
for s in ["top", "right"]:
    ax2.spines[s].set_visible(False)

fig.suptitle("Line 1  |  Transcriptomic Silencing of Mutations (Stealth Mode)",
             fontsize=15, fontweight="bold", color=INK, x=0.012, ha="left", y=0.985)
fig.tight_layout(rect=[0, 0, 1, 0.955])

out = "%s/line1_convergence.png" % FIG_DIR
fig.savefig(out, dpi=200, facecolor="white")
print("[3/3] wrote %s" % out, flush=True)
