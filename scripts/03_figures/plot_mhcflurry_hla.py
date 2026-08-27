#!/usr/bin/env python3
"""
How much of the presented repertoire depends on the genotype.

Two samples with the same mutations do not present the same neoantigens.
A patient carrying six distinct class I alleles has six shapes of binding
groove to work with; one carrying four has four. The left panel measures
that: 0.91 strong binders per mutation at four alleles, 1.74 at six.

Nearly double, from the genotype alone. That is the mechanism the project
is about, arriving here as a side effect of counting binders — and it is
also the reason losing one allele matters.

The right panel asks how concentrated each patient's repertoire is. If
the binders were spread evenly across six alleles, losing one would cost
a sixth. They are not: the single most productive allele carries a third
of them, and in one sample more than half.

Usage:
  python plot_mhcflurry_hla.py [outdir]
"""
import os
import sys
import glob
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from deckplot import *

WS = os.environ.get("WS", "/gpfs/bwfor/work/ws/fr_os136-immune_escape")
RES = os.path.expanduser("~/immune_escape_project/results")
OUT = sys.argv[1] if len(sys.argv) > 1 else f"{RES}/figures"
COHORT = f"{WS}/simulation/cohort"

d = pd.read_csv(f"{RES}/mhcflurry_summary.tsv", sep="\t")
d["per_mut"] = d.strong / d.mutations

frames = []
for p in glob.glob(f"{COHORT}/*/neoantigens/neoantigens_binders.tsv"):
    sid = p.split("/cohort/")[1].split("/")[0]
    t = pd.read_csv(p, sep="\t", usecols=lambda c: c in ("allele", "binder"))
    t["sample"] = sid
    frames.append(t)
B = pd.concat(frames) if frames else pd.DataFrame()
strong = B[B.binder == "STRONG"] if len(B) and "binder" in B else B

W = 62
print("=" * W)
print(" BINDERS PER MUTATION, BY ALLELE COUNT")
print("=" * W)
grp = d.groupby("n_alleles").agg(
    samples=("sample", "size"),
    per_mut=("per_mut", "median")).reset_index()
print(f"\n  {'alleles':<10}{'samples':>9}{'per mutation':>15}")
for _, r in grp.iterrows():
    print(f"  {int(r.n_alleles):<10}{int(r.samples):>9}{r.per_mut:>15.2f}")

conc = []
if len(strong):
    for sid, g in strong.groupby("sample"):
        vc = g.allele.value_counts()
        if len(vc):
            conc.append(100 * vc.iloc[0] / len(g))
    print(f"\n  share on each sample's best allele:")
    print(f"    median {np.median(conc):.0f}%, "
          f"range {min(conc):.0f}-{max(conc):.0f}%")

os.makedirs(OUT, exist_ok=True)
fig = figure(12, 4.6)

# ------------------------------- left: yield against genotype breadth
axl = panel(fig, 0.065, 0.375, bottom=0.20,
            title="strong binders per mutation")
xs = np.arange(len(grp))
vals = grp.per_mut.values

axl.bar(xs, vals, width=0.55, zorder=3, edgecolor="none",
        color=[SLATE, GOLD, EMBER][:len(grp)])

top = vals.max()
for x, r in zip(xs, grp.itertuples()):
    axl.text(x, r.per_mut + top * 0.05, f"{r.per_mut:.2f}", family=HEAD,
             fontsize=15, color=BONE, ha="center", va="bottom", zorder=6)
    axl.text(x, top * 0.05, f"{int(r.samples)} samples", family=BODY,
             fontsize=11, color="#2A1008" if r.per_mut > 1.5 else BONE,
             ha="center", va="bottom", zorder=6)

axl.set_xticks(xs)
axl.set_xticklabels([f"{int(n)} alleles" for n in grp.n_alleles])
axl.tick_params(axis="x", colors=ASH, labelsize=12, length=0, pad=8)
axl.set_ylim(0, top * 1.24)
axl.set_xlim(-0.65, len(grp) - 0.35)
finish(axl)

# ---------------------------- right: how concentrated the repertoire is
axr = panel(fig, 0.585, 0.375, bottom=0.20,
            title="share carried by the best allele")

if conc:
    conc = np.array(sorted(conc))
    xs2 = np.arange(len(conc))
    axr.bar(xs2, conc, width=0.85, zorder=3, edgecolor="none",
            color=[EMBER if v >= 40 else SLATE for v in conc])

    med = np.median(conc)
    axr.axhline(med, color=BONE, lw=1.3, ls=(0, (5, 4)), alpha=0.6,
                zorder=5)
    axr.text(0.5, med + 1.8, f"median {med:.0f}%", family=BODY,
             fontsize=11.5, color=ASH, va="bottom", zorder=6)

    axr.set_xticks([])
    axr.set_ylim(0, max(conc) * 1.20)
    axr.set_xlim(-1, len(conc))
    axr.yaxis.set_major_formatter(lambda v, p: f"{v:.0f}%")
    finish(axr, "one bar per sample")

print("\n" + save(fig, f"{OUT}/mhcflurry_hla.png"))
