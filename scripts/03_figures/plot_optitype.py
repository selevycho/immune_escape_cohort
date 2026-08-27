#!/usr/bin/env python3
"""
What OptiType produced across the cohort.

Three questions the results have to answer.

Are the genotypes distinct? Forty backbones were chosen so that none is
reused, so forty different genotypes should come back. Identical calls
across samples would mean the tool defaulted rather than typed — which is
exactly what happened on an earlier run, silently, when it was given fewer
threads than it wanted.

How much diversity is there? A cohort drawn from five superpopulations
should show alleles that are common in some and absent in others. A narrow
allele set would mean the panel is not capturing the variation.

Where is it thin? Three samples came back without HLA-C. Whether that is
coverage or read count decides whether it is fixable.

Usage:
  python plot_optitype.py [outdir]
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

man = pd.read_csv(f"{COHORT}/manifest.tsv", sep="\t")

rows = []
for _, m in man.iterrows():
    sid = m.sample_id
    hits = glob.glob(f"{COHORT}/{sid}/optitype/*result.tsv") + \
           glob.glob(f"{COHORT}/{sid}/optitype/*/*result.tsv")
    if not hits:
        rows.append({"sample": sid, "cohort": m.cohort,
                     "pop": m.get("superpopulation"), "ok": False})
        continue
    t = pd.read_csv(hits[0], sep="\t")
    if not len(t):
        rows.append({"sample": sid, "cohort": m.cohort,
                     "pop": m.get("superpopulation"), "ok": False})
        continue
    x = t.iloc[0]
    r = {"sample": sid, "cohort": m.cohort, "pop": m.get("superpopulation"),
         "ok": True,
         "reads": float(x.get("Reads", np.nan)),
         "objective": float(x.get("Objective", np.nan))}
    for loc in "ABC":
        for k in (1, 2):
            # an uncalled locus comes back as NaN, and str(NaN) is the
            # string "nan", which counts as an allele unless caught here
            val = x.get(f"{loc}{k}")
            r[f"{loc}{k}"] = (None if pd.isna(val)
                              or str(val).strip().lower() in ("", "nan")
                              else str(val).strip())
    rows.append(r)

d = pd.DataFrame(rows)
good = d[d.ok].copy()

ALLELE_COLS = [f"{l}{k}" for l in "ABC" for k in (1, 2)]
good["genotype"] = good[ALLELE_COLS].astype(str).agg("|".join, axis=1)

W = 68
print("=" * W)
print(" OPTITYPE ACROSS THE COHORT")
print("=" * W)
print(f"\n  samples typed              {len(good)} of {len(d)}")
print(f"  distinct genotypes         {good.genotype.nunique()}")
if good.genotype.nunique() < len(good):
    dup = good.genotype.value_counts()
    dup = dup[dup > 1]
    print(f"  repeated genotypes         {len(dup)}  <-- investigate")

print(f"\n  distinct alleles per locus:")
counts = {}
for loc in "ABC":
    a = pd.concat([good[f"{loc}1"], good[f"{loc}2"]]).dropna()
    a = a[a.notna()]
    counts[loc] = a.nunique()
    print(f"    HLA-{loc}   {a.nunique():>3} distinct   "
          f"{len(a)} calls   most common {a.value_counts().index[0]}"
          f" ({a.value_counts().iloc[0]}x)")

print(f"\n  heterozygous at:")
for loc in "ABC":
    het = (good[f"{loc}1"] != good[f"{loc}2"]) & good[f"{loc}1"].notna() \
          & good[f"{loc}2"].notna()
    print(f"    HLA-{loc}   {int(het.sum())} of {len(good)}")

missing = {}
for loc in "ABC":
    m_ = good[good[f"{loc}1"].isna()]
    missing[loc] = list(m_["sample"])
    if len(m_):
        print(f"\n  HLA-{loc} missing in {len(m_)}: {' '.join(m_['sample'])}")
        print(f"    reads reaching OptiType: "
              f"{', '.join(f'{v:.0f}' for v in m_.reads if np.isfinite(v))}")
        print(f"    cohort median: {good.reads.median():.0f}")

if "pop" in good.columns and good["pop"].notna().any():
    print(f"\n  by superpopulation:")
    for p, g in good.groupby("pop"):
        a = pd.concat([g.A1, g.A2, g.B1, g.B2, g.C1, g.C2]).dropna()
        print(f"    {str(p):<6}{len(g):>3} samples, {a.nunique():>3} distinct alleles")

# ------------------------------------------------------------------ figure
os.makedirs(OUT, exist_ok=True)
fig = figure(12, 4.6)

ax = panel(fig, 0.06, 0.36, bottom=0.19, title="distinct alleles seen")
xs = np.arange(3)
vals = [counts[l] for l in "ABC"]
ax.bar(xs, vals, width=0.5, color=[EMBER, FLAME, GOLD], zorder=3,
       edgecolor="none")
for x, v in zip(xs, vals):
    ax.text(x, v + max(vals) * 0.035, str(v), family=HEAD, fontsize=19,
            color=BONE, ha="center", va="bottom", zorder=6)
ax.set_xticks(xs)
ax.set_xticklabels([f"HLA-{l}" for l in "ABC"])
ax.tick_params(axis="x", colors=ASH, labelsize=13, length=0, pad=8)
headroom(ax, vals + [0], pad=0.20)
ax.set_xlim(-0.65, 2.65)
finish(ax)
note(ax, f"in {len(good)} samples", loc="upper right", color=DUSK, size=11)

ax2 = blank(fig, (0.55, 0.19, 0.43, 0.62))
fig.text(0.55, 0.845, "WHAT CAME BACK", family=HEAD, fontsize=14,
         color=BONE, fontweight="bold", va="bottom")

facts = [
    (f"{good.genotype.nunique()} of {len(good)}", "genotypes distinct",
     GREEN if good.genotype.nunique() == len(good) else EMBER),
    (f"{len(good)} of {len(d)}", "samples typed at all", GREEN),
    (f"{sum(len(v) for v in missing.values())}", "loci left uncalled",
     GOLD if any(missing.values()) else GREEN),
]
y = 0.74
for big, lab, col in facts:
    ax2.add_patch(FancyBboxPatch((0, y), 1.0, 0.21,
                                 boxstyle="round,pad=0,rounding_size=0.03",
                                 facecolor=INK, edgecolor=col, lw=1.5,
                                 zorder=3))
    ax2.text(0.035, y + 0.105, big, family=HEAD, fontsize=20, color=col,
             va="center", fontweight="bold", zorder=6)
    ax2.text(0.42, y + 0.105, lab, family=BODY, fontsize=13, color=ASH,
             va="center", zorder=6)
    y -= 0.28

print("\n" + save(fig, f"{OUT}/optitype_results.png"))
d.to_csv(f"{RES}/verify_step4.tsv", sep="\t", index=False)
print(f"written to {RES}/verify_step4.tsv")
