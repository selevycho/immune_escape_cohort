#!/usr/bin/env python3
"""
OptiType accuracy against the published 1000 Genomes types.

Two panels, because accuracy alone does not say what went wrong.

The left one is per-locus sensitivity: of the alleles the reference
records, how many did OptiType recover. B is much better than A and C,
which is the shape of the result.

The right one is the reason. Eight of the alleles OptiType reported appear
in none of the 2693 typed individuals in the panel — they cannot be right,
and they are what the wrong calls look like. Plotting the count of
impossible calls beside the accuracy turns "78%" from a bare number into
an explained one.

Every label is placed relative to the bar it belongs to, and axis limits
carry explicit headroom, so nothing can collide however the numbers move.

Usage:
  python plot_optitype_accuracy.py [outdir]
"""
import os
import re
import sys
from collections import Counter
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from deckplot import *

WS = os.environ.get("WS", "/gpfs/bwfor/work/ws/fr_os136-immune_escape")
RES = os.path.expanduser("~/immune_escape_project/results")
OUT = sys.argv[1] if len(sys.argv) > 1 else f"{RES}/figures"
ACC = f"{RES}/verify_step4_accuracy.tsv"
CACHE = f"{WS}/ref/1000G_HLA_types.txt"
COHORT = f"{WS}/simulation/cohort"

for p in (ACC, CACHE):
    if not os.path.exists(p):
        sys.exit(f"missing {p}")

d = pd.read_csv(ACC, sep="\t")
d = d[d.in_reference]
ref = pd.read_csv(CACHE, sep="\t", dtype=str)
man = pd.read_csv(f"{COHORT}/manifest.tsv", sep="\t")


def tf(a):
    if a is None or (isinstance(a, float) and pd.isna(a)):
        return None
    s = str(a).strip().split("/")[0].replace("HLA-", "")
    if not s or s.lower() in ("nan", "na", "-"):
        return None
    rest = s.split("*", 1)[1] if "*" in s else s
    p = rest.split(":")
    return f"{p[0].strip()}:{p[1].strip()}" if len(p) >= 2 \
        and p[0].strip().isdigit() else None


def cols_for(loc):
    return [c for c in ref.columns
            if re.fullmatch(rf"HLA[-_ ]?{loc}[ _-]*[12]", c.strip(), re.I)]


sens, impossible, uncalled = {}, {}, {}
for loc in "ABC":
    g = d[d[f"{loc}_n_truth"].notna()]
    sens[loc] = 100 * g[f"{loc}_correct"].sum() / max(1, g[f"{loc}_n_truth"].sum())
    uncalled[loc] = int((d[f"{loc}_bucket"] == "not called").sum())

    pool = Counter()
    for cc in cols_for(loc):
        for v in ref[cc]:
            a = tf(v)
            if a:
                pool[a] += 1
    n_bad = 0
    for _, m in man.iterrows():
        p = f"{COHORT}/{m.sample_id}/optitype/{m.sample_id}_result.tsv"
        if not os.path.exists(p):
            continue
        t = pd.read_csv(p, sep="\t")
        if not len(t):
            continue
        for k in (1, 2):
            a = tf(t.iloc[0].get(f"{loc}{k}"))
            if a and pool[a] == 0:
                n_bad += 1
    impossible[loc] = n_bad

tot_cor = int(sum(d[f"{l}_correct"].sum() for l in "ABC"))
tot_pub = int(sum(d[f"{l}_n_truth"].sum() for l in "ABC"))

print("=" * 62)
print(" OPTITYPE ACCURACY")
print("=" * 62)
print(f"\n  overall   {tot_cor} of {tot_pub}   ({100*tot_cor/tot_pub:.1f}%)")
for loc in "ABC":
    print(f"  HLA-{loc}     {sens[loc]:.1f}%   "
          f"{impossible[loc]} impossible calls   "
          f"{uncalled[loc]} uncalled")

os.makedirs(OUT, exist_ok=True)
fig = figure(12, 4.6)

# ---------------------------------------------------- left: sensitivity
ax = panel(fig, 0.065, 0.375, bottom=0.20,
           title="alleles recovered, of those published")
xs = np.arange(3)
vals = [sens[l] for l in "ABC"]
ax.bar(xs, vals, width=0.5, zorder=3, edgecolor="none",
       color=[GREEN if v >= 85 else (GOLD if v >= 75 else EMBER)
              for v in vals])
for x, v, loc in zip(xs, vals, "ABC"):
    g = d[d[f"{loc}_n_truth"].notna()]
    ax.text(x, v + 2.6, f"{v:.1f}%", family=HEAD, fontsize=17, color=BONE,
            ha="center", va="bottom", zorder=6)
    ax.text(x, 4, f"{int(g[f'{loc}_correct'].sum())} of "
                  f"{int(g[f'{loc}_n_truth'].sum())}", family=BODY,
            fontsize=11.5, color="#2A1008" if v >= 85 else BONE,
            ha="center", va="bottom", zorder=6)
ax.set_xticks(xs)
ax.set_xticklabels([f"HLA-{l}" for l in "ABC"])
ax.tick_params(axis="x", colors=ASH, labelsize=13, length=0, pad=8)
ax.set_ylim(0, 112)
ax.set_xlim(-0.65, 2.65)
ax.yaxis.set_major_formatter(lambda v, p: f"{v:.0f}%")
finish(ax)
note(ax, f"overall {100*tot_cor/tot_pub:.1f}%", loc="upper right", size=12)

# ------------------------------------------- right: impossible calls
ax2 = panel(fig, 0.575, 0.385, bottom=0.20,
            title="calls seen in no typed individual")
vals2 = [impossible[l] for l in "ABC"]
ax2.bar(xs, vals2, width=0.5, color=EMBER, zorder=3, edgecolor="none")
for x, v in zip(xs, vals2):
    ax2.text(x, v + max(vals2) * 0.06, str(v), family=HEAD, fontsize=18,
             color=BONE, ha="center", va="bottom", zorder=6)
ax2.set_xticks(xs)
ax2.set_xticklabels([f"HLA-{l}" for l in "ABC"])
ax2.tick_params(axis="x", colors=ASH, labelsize=13, length=0, pad=8)
headroom(ax2, vals2 + [0], pad=0.28)
ax2.set_xlim(-0.65, 2.65)
finish(ax2)
note(ax2, f"out of {len(ref)} people", loc="upper right", color=DUSK,
     size=11)

print("\n" + save(fig, f"{OUT}/optitype_accuracy.png"))
