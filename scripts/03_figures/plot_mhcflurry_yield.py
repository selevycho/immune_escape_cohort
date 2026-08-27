#!/usr/bin/env python3
"""
What came out of the binding prediction.

The funnel first: fifty thousand peptides tested, seven and a half
thousand binding at all, eighteen hundred binding strongly. Those numbers
mean little on their own — what matters is that 89% of the missense
mutations produced at least one binder, so the input to the escape
question is not thin.

The length panel is the sanity check. Nine residues is the canonical
class I length and the predictor was trained accordingly; if nine-mers
did not dominate, something would be wrong with how the peptides were
cut. They do, by a factor of four.

Usage:
  python plot_mhcflurry_yield.py [outdir]
"""
import os
import sys
import glob
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from deckplot import *
from matplotlib.patches import FancyBboxPatch

WS = os.environ.get("WS", "/gpfs/bwfor/work/ws/fr_os136-immune_escape")
RES = os.path.expanduser("~/immune_escape_project/results")
OUT = sys.argv[1] if len(sys.argv) > 1 else f"{RES}/figures"
COHORT = f"{WS}/simulation/cohort"

d = pd.read_csv(f"{RES}/mhcflurry_summary.tsv", sep="\t")

frames = []
for p in glob.glob(f"{COHORT}/*/neoantigens/neoantigens_all.tsv"):
    frames.append(pd.read_csv(p, sep="\t",
                              usecols=lambda c: c in ("length", "binder")))
P = pd.concat(frames) if frames else pd.DataFrame()

n_mut = int(d.mutations.sum())
n_pep = int(d.peptides.sum())
n_strong = int(d.strong.sum())
n_weak = int(d.weak.sum())
n_yield = int(d.with_binder.sum())

W = 60
print("=" * W)
print(" YIELD")
print("=" * W)
print(f"\n  missense mutations   {n_mut}")
print(f"  peptides tested      {n_pep:,}")
print(f"  weak binders         {n_weak}")
print(f"  strong binders       {n_strong}")
print(f"\n  mutations with a binder  {n_yield} of {n_mut}"
      f"   ({100*n_yield/n_mut:.1f}%)")

os.makedirs(OUT, exist_ok=True)
fig = figure(12, 4.6)

# --------------------------------------------------- left: the funnel
axl = blank(fig, (0.045, 0.10, 0.42, 0.78))
fig.text(0.045, 0.905, "FROM PEPTIDES TO BINDERS", family=HEAD,
         fontsize=14, color=BONE, fontweight="bold", va="bottom")

stages = [("peptides tested", n_pep, SLATE),
          ("bind at all", n_weak + n_strong, GOLD),
          ("bind strongly", n_strong, EMBER)]
widest = stages[0][1]

BY, BH, GAP = 0.66, 0.19, 0.09
for i, (label, k, colour) in enumerate(stages):
    y = BY - i * (BH + GAP)
    w = max(k / widest, 0.06)
    axl.add_patch(FancyBboxPatch((0, y), w, BH,
                                 boxstyle="round,pad=0,rounding_size=0.025",
                                 facecolor=colour, edgecolor="none",
                                 zorder=3))
    # the label goes inside a wide bar and outside a narrow one, so it
    # never has to compete with the number beside it
    if w > 0.42:
        axl.text(0.022, y + BH / 2, label.upper(), family=HEAD,
                 fontsize=12, color="#1A1412", va="center",
                 fontweight="bold", zorder=6)
        axl.text(w - 0.022, y + BH / 2, f"{k:,}", family=HEAD,
                 fontsize=16, color="#1A1412", ha="right", va="center",
                 fontweight="bold", zorder=6)
    else:
        axl.text(w + 0.028, y + BH / 2, f"{k:,}", family=HEAD,
                 fontsize=16, color=BONE, va="center",
                 fontweight="bold", zorder=6)
        axl.text(w + 0.028, y - 0.045, label, family=BODY, fontsize=11.5,
                 color=DUSK, va="center", zorder=6)

# the lowest bar bottoms out at 0.10 and its caption sits at 0.055, so
# this line has to clear both
axl.text(0, -0.03, f"{n_yield} of {n_mut} mutations produced at least one",
         family=BODY, fontsize=12, color=GOLD, va="center")

# ------------------------------------------------- right: by length
axr = panel(fig, 0.585, 0.375, bottom=0.20,
            title="strong binders, by peptide length")

if len(P) and "length" in P.columns:
    rows = []
    for L, g in P.groupby("length"):
        ns = int((g.binder == "STRONG").sum())
        rows.append({"length": int(L), "tested": len(g), "strong": ns,
                     "rate": 100 * ns / len(g)})
    t = pd.DataFrame(rows).sort_values("length")

    print(f"\n  {'length':<9}{'tested':>10}{'strong':>9}{'rate':>9}")
    for _, r in t.iterrows():
        print(f"  {int(r.length):<9}{int(r.tested):>10}"
              f"{int(r.strong):>9}{r.rate:>8.2f}%")

    xs = np.arange(len(t))
    axr.bar(xs, t.strong, width=0.55, zorder=3, edgecolor="none",
            color=[EMBER if L == 9 else SLATE for L in t.length])

    top = t.strong.max()
    for x, r in zip(xs, t.itertuples()):
        axr.text(x, r.strong + top * 0.045, f"{int(r.strong)}",
                 family=HEAD, fontsize=14, color=BONE, ha="center",
                 va="bottom", zorder=6)
        axr.text(x, top * 0.045, f"{r.rate:.2f}%", family=BODY,
                 fontsize=11,
                 color="#2A1008" if r.length == 9 else BONE,
                 ha="center", va="bottom", zorder=6)

    axr.set_xticks(xs)
    axr.set_xticklabels([f"{int(L)}-mer" for L in t.length])
    axr.tick_params(axis="x", colors=ASH, labelsize=12, length=0, pad=8)
    axr.set_ylim(0, top * 1.22)
    axr.set_xlim(-0.65, len(t) - 0.35)
    finish(axr)

print("\n" + save(fig, f"{OUT}/mhcflurry_yield.png"))
