#!/usr/bin/env python3
"""
The sensitivity curve, and the ceiling above it.

Two things are being said at once. Recall climbs with allele fraction, as
it must — a mutation in fewer reads is harder to distinguish from a
sequencing error. And it stops climbing at 94%, which is the more
interesting half: six per cent of mutations sitting in thirteen reads at
full depth are still not recovered.

Binning is by observed fraction rather than by requested. The injection
lands at about 0.94 of what was asked, so a curve drawn against the
request would be shifted by that much and every bin would understate its
own label.

The reads-per-bin axis is the point of the whole figure. Below three
supporting reads nothing is recoverable by any caller, and at 34x that
floor arrives at 9% allele fraction. It is arithmetic, not a property of
Mutect2.

Usage:
  python plot_mutect_recall.py [outdir]
"""
import os
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from deckplot import *

RES = os.path.expanduser("~/immune_escape_project/results")
OUT = sys.argv[1] if len(sys.argv) > 1 else f"{RES}/figures"
SRC = f"{RES}/step3_rescored_full.tsv"

if not os.path.exists(SRC):
    sys.exit(f"missing {SRC}")

d = pd.read_csv(SRC, sep="\t")
BINS = [0, .05, .10, .15, .20, .30, 1.01]
LABELS = ["<5%", "5-10%", "10-15%", "15-20%", "20-30%", ">30%"]
d["bin"] = pd.cut(d.observed_vaf, bins=BINS, labels=LABELS, right=False)

rows = []
for lab in LABELS:
    g = d[d.bin == lab]
    if not len(g):
        continue
    rows.append({"bin": lab, "n": len(g),
                 "emitted": int(g.emitted.sum()),
                 "passed": int(g["pass"].sum()),
                 "recall": 100 * g["pass"].mean(),
                 "reads": g.alt_reads.median()})
t = pd.DataFrame(rows)

print("=" * 62)
print(" SENSITIVITY BY OBSERVED ALLELE FRACTION")
print("=" * 62)
print(f"\n  {'fraction':<11}{'n':>6}{'emitted':>10}{'PASS':>7}"
      f"{'recall':>9}{'reads':>8}")
for _, r in t.iterrows():
    print(f"  {r['bin']:<11}{int(r.n):>6}{int(r.emitted):>10}"
          f"{int(r.passed):>7}{r.recall:>8.1f}%{r.reads:>8.0f}")
print(f"\n  overall  {int(d['pass'].sum())} of {len(d)}"
      f"   ({100*d['pass'].mean():.1f}%)")
print(f"  ceiling  {t.recall.max():.1f}%")

os.makedirs(OUT, exist_ok=True)
fig = figure(12, 4.6)

# ------------------------------------------------- left: the curve
ax = panel(fig, 0.065, 0.40, bottom=0.19,
           title="recovered by Mutect2")
xs = np.arange(len(t))

ax.bar(xs, t.recall, width=0.55, zorder=3, edgecolor="none",
       color=[GREEN if v >= 90 else (GOLD if v >= 60 else EMBER)
              for v in t.recall])

# the ceiling: where the curve stops rising
# the ceiling line is labelled on the left, where no bar reaches it, so
# it cannot collide with the value printed above the tallest bar
ceil = t.recall.max()
ax.axhline(ceil, color=BONE, lw=1.2, ls=(0, (5, 4)), alpha=0.55, zorder=5)
ax.text(-0.5, ceil + 2, f"ceiling {ceil:.0f}%", family=BODY, fontsize=11.5,
        color=ASH, ha="left", va="bottom", zorder=6)

for x, r in zip(xs, t.itertuples()):
    ax.text(x, r.recall + 2.5, f"{r.recall:.0f}%", family=HEAD,
            fontsize=14, color=BONE, ha="center", va="bottom", zorder=6)
    # a bar of no height has nowhere to put its count without landing on
    # the percentage above it; that one goes on the slide by hand
    if r.recall >= 12:
        ax.text(x, 3.5, f"{int(r.n)}", family=BODY, fontsize=11,
                color="#2A1008" if r.recall >= 60 else BONE,
                ha="center", va="bottom", zorder=6)

ax.set_xticks(xs)
ax.set_xticklabels(t["bin"])
ax.tick_params(axis="x", colors=ASH, labelsize=12, length=0, pad=8)
ax.set_ylim(0, 112)
ax.set_xlim(-0.65, len(t) - 0.35)
ax.yaxis.set_major_formatter(lambda v, p: f"{v:.0f}%")
finish(ax, "allele fraction in the reads")
ax.set_title("")

# --------------------------------------- right: reads, and the floor
ax2 = panel(fig, 0.575, 0.385, bottom=0.19,
            title="supporting reads at that fraction")
ax2.bar(xs, t.reads, width=0.55, color=SLATE, zorder=3, edgecolor="none")

# three reads is where a variant stops being separable from an error
ax2.axhline(3, color=EMBER, lw=1.6, ls=(0, (5, 4)), zorder=5)
ax2.text(len(t) - 0.4, 3.6, "three reads", family=BODY, fontsize=11.5,
         color=EMBER, ha="right", va="bottom", zorder=6)

for x, r in zip(xs, t.itertuples()):
    ax2.text(x, r.reads + 0.35, f"{r.reads:.0f}", family=BODY,
             fontsize=12.5, color=BONE, ha="center", va="bottom", zorder=6)

ax2.set_xticks(xs)
ax2.set_xticklabels(t["bin"])
ax2.tick_params(axis="x", colors=ASH, labelsize=12, length=0, pad=8)
headroom(ax2, list(t.reads) + [0], pad=0.22)
ax2.set_xlim(-0.65, len(t) - 0.35)
finish(ax2, "allele fraction in the reads")
ax2.set_title("")

print("\n" + save(fig, f"{OUT}/mutect_recall.png"))
