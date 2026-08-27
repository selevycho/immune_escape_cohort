#!/usr/bin/env python3
"""
Clean redraw of the two injection figures.

The first version of injection_losses.png computed segment widths
proportional to their share and placed labels under them unconditionally.
At 2.6-3.6% per segment the labels are wider than the segments and
collided. This version uses fixed-width cards instead - width carries no
information here, the printed numbers do, and cards cannot overlap
regardless of how small a share is.

The VAF figure was structurally fine; its titles sat outside the axes and
got clipped by bbox_inches='tight'. This version reserves headroom for
them explicitly.

Usage:
  python plot_injection_clean.py [outdir]
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.patches import FancyBboxPatch

RES = os.path.expanduser("~/immune_escape_project/results")
OUT = sys.argv[1] if len(sys.argv) > 1 else f"{RES}/figures"
SNV = f"{RES}/verify_step2_per_mutation.tsv"

EMBER, GOLD, FLAME = "#E8402A", "#F2A623", "#F27127"
BONE, ASH, DUSK = "#F4EFEA", "#A29591", "#7A6C68"
GREEN, RUST, SLATE = "#6F9E44", "#8E3020", "#5A4744"

have = {f.name for f in fm.fontManager.ttflist}
HEAD = next((f for f in ["TeX Gyre Heros Cn", "Liberation Sans Narrow",
                         "Liberation Sans", "DejaVu Sans"] if f in have),
            "DejaVu Sans")
BODY = next((f for f in ["Carlito", "Open Sans", "Liberation Sans",
                         "DejaVu Sans"] if f in have), "DejaVu Sans")

d = pd.read_csv(SNV, sep="\t")
d["alt"] = d.alt.astype(str)
snv = d[(d.alt.str.len() == 1) & (d.alt != "-")].copy()

n = len(snv)
landed = int(snv.landed.sum())
lost = snv[~snv.landed]

cats = [
    ("MHC", int(lost.in_mhc.sum()), RUST),
    ("Low coverage", int(((~lost.in_mhc) & (lost.tumour_depth < 5)).sum()), SLATE),
    ("Covered, no alt read",
     int(((~lost.in_mhc) & (lost.tumour_depth >= 5)).sum()), EMBER),
]

os.makedirs(OUT, exist_ok=True)

# ===================================================================== 1
fig = plt.figure(figsize=(12, 3.6), dpi=300)
fig.patch.set_alpha(0)
ax = fig.add_axes([0.03, 0.06, 0.94, 0.80])
ax.axis("off")
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)

ax.text(0, 1.0, f"{n} SUBSTITUTIONS ATTEMPTED", family=HEAD, fontsize=15,
        color=BONE, va="top", fontweight="bold")

# one wide card for what landed
LW, LH, LY = 0.985, 0.30, 0.52
ax.add_patch(FancyBboxPatch((0, LY), LW, LH,
                            boxstyle="round,pad=0,rounding_size=0.02",
                            facecolor=GREEN, edgecolor="none", zorder=3))
ax.text(0.02, LY + LH / 2, "LANDED", family=HEAD, fontsize=16,
        color="#14200A", va="center", fontweight="bold", zorder=5)
ax.text(0.985, LY + LH / 2, f"{landed}  ({100*landed/n:.1f}%)",
        family=HEAD, fontsize=19, color="#14200A", ha="right", va="center",
        fontweight="bold", zorder=5)

# three fixed-width cards for the losses - width is decorative, not
# proportional, so three thin true segments never have to share a label
CW, GAP = 0.315, 0.028
CY, CH = 0.02, 0.36
x = 0
for name, k, col in cats:
    ax.add_patch(FancyBboxPatch((x, CY), CW, CH,
                                boxstyle="round,pad=0,rounding_size=0.02",
                                facecolor="#241A18", edgecolor=col, lw=1.6,
                                zorder=3))
    ax.add_patch(FancyBboxPatch((x, CY), 0.012, CH,
                                boxstyle="round,pad=0,rounding_size=0.006",
                                facecolor=col, edgecolor="none", zorder=4))
    ax.text(x + 0.03, CY + CH - 0.10, name.upper(), family=HEAD,
            fontsize=11.5, color=ASH, va="center", fontweight="bold", zorder=5)
    ax.text(x + 0.03, CY + 0.10,
            f"{k}   ({100*k/len(lost):.0f}% of losses)", family=BODY,
            fontsize=12.5, color=BONE, va="center", zorder=5)
    x += CW + GAP

f1 = f"{OUT}/injection_losses.png"
fig.savefig(f1, dpi=300, transparent=True, bbox_inches="tight",
            pad_inches=0.14)
plt.close(fig)
print(f"written to {f1}")

# ===================================================================== 2
ok = snv[snv.landed & (snv.target_vaf > 0)]
r = np.corrcoef(ok.target_vaf, ok.observed_vaf)[0, 1]

bins = [0, .05, .10, .15, .20, .30, 1.01]
labs = ["5-10%", "10-15%", "15-20%", "20-30%", ">30%"]
ok = ok.copy()
ok["bin"] = pd.cut(ok.target_vaf, bins=bins, labels=["<5%"] + labs,
                   right=False)
devs = [(ok[ok.bin == l].observed_vaf - ok[ok.bin == l].target_vaf).mean()
        for l in labs if len(ok[ok.bin == l]) > 3]
labs = [l for l in labs if len(ok[ok.bin == l]) > 3]

fig = plt.figure(figsize=(12, 4.6), dpi=300)
fig.patch.set_alpha(0)

# left: scatter, with headroom reserved above the axes for the title
ax = fig.add_axes([0.07, 0.14, 0.40, 0.68])
ax.set_facecolor((0, 0, 0, 0))
for s in ax.spines.values():
    s.set_visible(False)
ax.scatter(ok.target_vaf, ok.observed_vaf, s=8, color=EMBER, alpha=0.35,
           edgecolors="none", zorder=3)
lim = max(ok.target_vaf.max(), ok.observed_vaf.max()) * 1.05
ax.plot([0, lim], [0, lim], color=GOLD, lw=1.3, ls=(0, (5, 4)), zorder=5)
ax.set_xlim(0, lim)
ax.set_ylim(0, lim)
ax.tick_params(colors=DUSK, labelsize=11, length=0)
for t in ax.get_xticklabels() + ax.get_yticklabels():
    t.set_fontfamily(BODY)
ax.grid(color="#3A2F2D", lw=0.6, zorder=0)
ax.set_axisbelow(True)
ax.set_xlabel("requested", color=ASH, fontsize=12, family=BODY, labelpad=8)
ax.set_ylabel("observed in the reads", color=ASH, fontsize=12, family=BODY,
              labelpad=8)
fig.text(0.07, 0.90, f"ALLELE FRACTION  ·  r = {r:.3f}", family=HEAD,
         fontsize=14, color=BONE, fontweight="bold")

# right: deviation by bin
ax2 = fig.add_axes([0.57, 0.14, 0.40, 0.68])
ax2.set_facecolor((0, 0, 0, 0))
for s in ax2.spines.values():
    s.set_visible(False)
xs = np.arange(len(labs))
ax2.bar(xs, devs, width=0.6, color=EMBER, zorder=3, edgecolor="none")
ax2.axhline(0, color=SLATE, lw=1.0, zorder=4)
ax2.set_xticks(xs)
ax2.set_xticklabels(labs)
ax2.tick_params(axis="x", colors=ASH, labelsize=12, length=0, pad=8)
ax2.tick_params(axis="y", colors=DUSK, labelsize=11, length=0)
for t in ax2.get_xticklabels() + ax2.get_yticklabels():
    t.set_fontfamily(BODY)
ax2.grid(axis="y", color="#3A2F2D", lw=0.6, zorder=0)
ax2.set_axisbelow(True)
ax2.set_xlabel("requested fraction", color=ASH, fontsize=12, family=BODY,
               labelpad=8)
fig.text(0.57, 0.90, "DEVIATION FROM WHAT WAS ASKED", family=HEAD,
         fontsize=14, color=BONE, fontweight="bold")

f2 = f"{OUT}/injection_vaf.png"
fig.savefig(f2, dpi=300, transparent=True, bbox_inches="tight",
            pad_inches=0.14)
plt.close(fig)
print(f"written to {f2}")
