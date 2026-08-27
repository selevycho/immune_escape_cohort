#!/usr/bin/env python3
"""
The second route, and what it recovers.

pVACseq starts where mhcflurry does not: from the variants Mutect2
actually called, annotated by VEP, rather than from the mutations that
were placed. It scores 193 304 epitopes and finds 1 255 strong binders,
against mhcflurry's 1 830.

The gene count is the more useful number of the two. A patient's clinical
question is not how many peptides bind but how many distinct proteins
present something, and the median sample here has seven such genes.

Usage:
  python plot_netmhcpan.py [outdir]
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
PVAC = f"{WS}/simulation/pvacseq/cohort"

pv = pd.read_csv(f"{RES}/pvacseq_summary.tsv", sep="\t")

frames = []
for p in glob.glob(f"{PVAC}/*/pvacseq_out/MHC_Class_I/*.all_epitopes.tsv"):
    t = pd.read_csv(p, sep="\t", low_memory=False)
    col = "Best MT Percentile"
    if col not in t.columns:
        col = next((c for c in t.columns if "Percentile" in c), None)
    if col is None or "HLA Allele" not in t.columns:
        continue
    v = pd.to_numeric(t[col], errors="coerce")
    frames.append(pd.DataFrame({"allele": t["HLA Allele"], "p": v}))
A = pd.concat(frames) if frames else pd.DataFrame()

W = 60
print("=" * W)
print(" NetMHCpan")
print("=" * W)
print(f"\n  samples          {len(pv)}")
print(f"  epitopes         {int(pv.epitopes.sum()):,}")
print(f"  strong binders   {int(pv.strong.sum())}")
print(f"  weak binders     {int(pv.weak.sum())}")
if "genes" in pv.columns:
    print(f"\n  genes with a strong binder")
    print(f"    median per sample   {pv.genes.median():.0f}")
    print(f"    range               {int(pv.genes.min())}-{int(pv.genes.max())}")

os.makedirs(OUT, exist_ok=True)
fig = figure(12, 4.6)

# --------------------------------------------- left: genes per sample
axl = panel(fig, 0.065, 0.375, bottom=0.20,
            title="genes presenting a strong binder")

g = np.array(sorted(pv.genes.values)) if "genes" in pv.columns else np.array([])
if len(g):
    xs = np.arange(len(g))
    axl.bar(xs, g, width=0.85, zorder=3, edgecolor="none",
            color=[EMBER if v >= 20 else SLATE for v in g])
    med = np.median(g)
    axl.axhline(med, color=BONE, lw=1.3, ls=(0, (5, 4)), alpha=0.6, zorder=5)
    axl.text(0.5, med + max(g) * 0.03, f"median {med:.0f}", family=BODY,
             fontsize=11.5, color=ASH, va="bottom", zorder=6)
    axl.set_xticks([])
    axl.set_ylim(0, max(g) * 1.18)
    axl.set_xlim(-1, len(g))
    finish(axl, "one bar per sample")

# ------------------------------------------------- right: by locus
axr = panel(fig, 0.585, 0.375, bottom=0.20,
            title="strong binders, by locus")

if len(A):
    s = A[A.p < 0.5]
    by_locus = s.allele.str.extract(r"HLA-([ABC])")[0].value_counts()
    by_locus = by_locus.reindex(["A", "B", "C"]).fillna(0)

    print(f"\n  by locus:")
    for loc, k in by_locus.items():
        print(f"    HLA-{loc}   {int(k):>5}   ({100*k/by_locus.sum():.0f}%)")

    xs2 = np.arange(3)
    vals = by_locus.values
    axr.bar(xs2, vals, width=0.5, zorder=3, edgecolor="none",
            color=[SLATE, EMBER, GOLD])
    top = vals.max()
    for x, v in zip(xs2, vals):
        axr.text(x, v + top * 0.04, f"{int(v)}", family=HEAD, fontsize=15,
                 color=BONE, ha="center", va="bottom", zorder=6)
        axr.text(x, top * 0.04, f"{100*v/vals.sum():.0f}%", family=BODY,
                 fontsize=11.5, color=BONE, ha="center", va="bottom",
                 zorder=6)
    axr.set_xticks(xs2)
    axr.set_xticklabels([f"HLA-{l}" for l in "ABC"])
    axr.tick_params(axis="x", colors=ASH, labelsize=12.5, length=0, pad=8)
    axr.set_ylim(0, top * 1.20)
    axr.set_xlim(-0.65, 2.65)
    finish(axr)

print("\n" + save(fig, f"{OUT}/netmhcpan_yield.png"))
