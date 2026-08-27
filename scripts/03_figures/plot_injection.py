#!/usr/bin/env python3
"""
The two injection figures, drawn from the verification tables.

Usage:
  python plot_injection.py [outdir]
"""
import os
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from deckplot import *
from matplotlib.patches import FancyBboxPatch

RES = os.path.expanduser("~/immune_escape_project/results")
OUT = sys.argv[1] if len(sys.argv) > 1 else f"{RES}/figures"
SRC = f"{RES}/verify_step2_per_mutation.tsv"

if not os.path.exists(SRC):
    sys.exit(f"missing {SRC}\nrun verify_step2_injection.py first")

d = pd.read_csv(SRC, sep="\t")
d["alt"] = d.alt.astype(str)
# indels went through a separate routine in step 8 and are counted there
snv = d[(d.alt.str.len() == 1) & (d.alt != "-")].copy()

n = len(snv)
landed = int(snv.landed.sum())
lost = snv[~snv.landed]
CATS = [
    ("MHC", int(lost.in_mhc.sum()), RUST),
    ("Low coverage", int(((~lost.in_mhc) & (lost.tumour_depth < 5)).sum()), SLATE),
    ("Covered, no alt read",
     int(((~lost.in_mhc) & (lost.tumour_depth >= 5)).sum()), EMBER),
]

print(f"  attempted {n}, landed {landed} ({100*landed/n:.1f}%)")
for name, k, _ in CATS:
    print(f"    {name:<24}{k:>5}   {100*k/len(lost):.0f}% of losses")

# ---------------------------------------------------------------- losses
fig = figure(12, 3.4)
ax = blank(fig, (0.03, 0.05, 0.94, 0.78))
fig.text(0.03, 0.90, f"{n} SUBSTITUTIONS ATTEMPTED", family=HEAD,
         fontsize=15, color=BONE, fontweight="bold", va="bottom")

LH, LY = 0.34, 0.58
ax.add_patch(FancyBboxPatch((0, LY), 1.0, LH,
                            boxstyle="round,pad=0,rounding_size=0.022",
                            facecolor=GREEN, edgecolor="none", zorder=3))
ax.text(0.022, LY + LH / 2, "LANDED", family=HEAD, fontsize=16,
        color="#14200A", va="center", fontweight="bold", zorder=5)
ax.text(0.978, LY + LH / 2, f"{landed}   {100*landed/n:.1f}%", family=HEAD,
        fontsize=19, color="#14200A", ha="right", va="center",
        fontweight="bold", zorder=5)

CW, GAP, CH, CY = 0.315, 0.0275, 0.40, 0.04
x = 0
for name, k, col in CATS:
    card(ax, x, CY, CW, CH, name, f"{k}   ({100*k/len(lost):.0f}% of losses)",
         accent=col)
    x += CW + GAP
print(save(fig, f"{OUT}/injection_losses.png"))

# ------------------------------------------------------------------- vaf
ok = snv[snv.landed & (snv.target_vaf > 0)]
r = np.corrcoef(ok.target_vaf, ok.observed_vaf)[0, 1]

BINS = [0, .05, .10, .15, .20, .30, 1.01]
NAMES = ["<5%", "5-10%", "10-15%", "15-20%", "20-30%", ">30%"]
ok = ok.copy()
ok["bin"] = pd.cut(ok.target_vaf, bins=BINS, labels=NAMES, right=False)
labs, devs = [], []
for nm in NAMES:
    g = ok[ok.bin == nm]
    if len(g) > 3:
        labs.append(nm)
        devs.append(float((g.observed_vaf - g.target_vaf).mean()))

print(f"  correlation r = {r:.3f}")
for nm, v in zip(labs, devs):
    print(f"    {nm:<10}{v:+.4f}")

fig = figure(12, 4.6)
ax = panel(fig, 0.07, 0.38, bottom=0.17, title="allele fraction")
ax.scatter(ok.target_vaf, ok.observed_vaf, s=8, color=EMBER, alpha=0.32,
           edgecolors="none", zorder=3)
lim = max(ok.target_vaf.max(), ok.observed_vaf.max()) * 1.06
ax.plot([0, lim], [0, lim], color=GOLD, lw=1.3, ls=(0, (5, 4)), zorder=5)
ax.set_xlim(0, lim)
ax.set_ylim(0, lim)
note(ax, f"r = {r:.3f}", loc="upper left")
finish(ax, "requested", "observed in the reads", grid="both")

ax2 = panel(fig, 0.57, 0.39, bottom=0.17,
            title="deviation from what was asked")
xs = np.arange(len(labs))
ax2.bar(xs, devs, width=0.58, color=EMBER, zorder=3, edgecolor="none")
ax2.axhline(0, color=SLATE, lw=1.1, zorder=4)
headroom(ax2, devs, pad=0.22)
bar_labels(ax2, xs, devs)
ax2.set_xticks(xs)
ax2.set_xticklabels(labs)
ax2.tick_params(axis="x", colors=ASH, labelsize=12, length=0, pad=7)
ax2.set_xlim(-0.7, len(labs) - 0.3)
finish(ax2, "requested fraction")
print(save(fig, f"{OUT}/injection_vaf.png"))
