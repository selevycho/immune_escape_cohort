#!/usr/bin/env python3
"""
Do the two predictors agree?

They correlate at r = 0.975 across samples while returning totals in a
ratio of 0.69. Both facts matter and they are not in tension: the two rank
samples identically and count different numbers of things.

The difference has two causes and neither is a disagreement about
binding. mhcflurry sees every placed mutation; pVACseq sees only what
Mutect2 recovered, which at 74.5% recall is a quarter fewer to work from.
And they count different objects — one row per peptide-allele pair against
one per epitope per transcript.

Agreement on ranking from two routes that share no code is the useful
result. Neither is carrying a systematic error the other would expose.

Usage:
  python plot_predictor_agreement.py [outdir]
"""
import os
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from deckplot import *

RES = os.path.expanduser("~/immune_escape_project/results")
OUT = sys.argv[1] if len(sys.argv) > 1 else f"{RES}/figures"

mf = pd.read_csv(f"{RES}/mhcflurry_summary.tsv", sep="\t")
pv = pd.read_csv(f"{RES}/pvacseq_summary.tsv", sep="\t")

j = pv.merge(mf[["sample", "strong", "mutations"]], on="sample",
             suffixes=("_pv", "_mf"))
j["ratio"] = j.strong_pv / j.strong_mf

r = float(np.corrcoef(j.strong_mf, j.strong_pv)[0, 1])

W = 60
print("=" * W)
print(" THE TWO ROUTES")
print("=" * W)
print(f"\n  samples compared   {len(j)}")
print(f"  mhcflurry total    {int(j.strong_mf.sum())}")
print(f"  NetMHCpan total    {int(j.strong_pv.sum())}")
print(f"  ratio              {j.strong_pv.sum()/j.strong_mf.sum():.2f}")
print(f"  correlation        r = {r:.3f}")
print(f"  median per-sample ratio  {j.ratio.median():.2f}")

os.makedirs(OUT, exist_ok=True)
fig = figure(12, 4.6)

# -------------------------------------------- left: one against the other
axl = panel(fig, 0.065, 0.375, bottom=0.20,
            title="strong binders, one route against the other")

axl.scatter(j.strong_mf, j.strong_pv, s=52, color=EMBER, alpha=0.8,
            edgecolors="none", zorder=4)

lim = max(j.strong_mf.max(), j.strong_pv.max()) * 1.08
axl.plot([0, lim], [0, lim], color=BONE, lw=1.2, ls=(0, (5, 4)),
         alpha=0.5, zorder=3)
axl.set_xlim(0, lim)
axl.set_ylim(0, lim)
note(axl, f"r = {r:.3f}", loc="upper left")
finish(axl, "mhcflurry", "NetMHCpan", grid="both")

# ----------------------------------------------- right: the ratio spread
axr = panel(fig, 0.585, 0.375, bottom=0.20,
            title="NetMHCpan as a fraction of mhcflurry")

ratios = np.array(sorted(j.ratio.values))
xs = np.arange(len(ratios))
axr.bar(xs, ratios, width=0.85, zorder=3, edgecolor="none",
        color=[GOLD if v > 1 else SLATE for v in ratios])

axr.axhline(1.0, color=BONE, lw=1.3, ls=(0, (5, 4)), alpha=0.6, zorder=5)
med = np.median(ratios)
axr.axhline(med, color=EMBER, lw=1.6, zorder=5)
axr.text(0.5, med + 0.05, f"median {med:.2f}", family=BODY, fontsize=11.5,
         color=EMBER, va="bottom", zorder=6)

axr.set_xticks([])
axr.set_ylim(0, max(ratios) * 1.15)
axr.set_xlim(-1, len(ratios))
finish(axr, "one bar per sample")

print("\n" + save(fig, f"{OUT}/predictor_agreement.png"))
