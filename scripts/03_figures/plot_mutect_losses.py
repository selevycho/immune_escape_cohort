#!/usr/bin/env python3
"""
Where the missing quarter goes.

Three hundred and eighty-one verified substitutions were not recovered,
and they divide into two groups that mean different things. Two hundred
and eighty-seven were called and then rejected by a filter: the caller
found them and decided against them, which is a threshold and can be
moved. Ninety-four never appeared in the output at all: the caller did not
consider them, which is a detection limit and cannot.

A single recall figure hides that distinction entirely, and the
distinction is the whole reason the next section — the filter sweep — has
anything to say.

The filter tags are worth showing rather than summarising. Two account for
almost everything: weak_evidence, which is a statement about read count,
and strand_bias, which is a statement about where in the read the
mutation sits. Neither is wrong in general.

Usage:
  python plot_mutect_losses.py [outdir]
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
SRC = f"{RES}/step3_rescored_full.tsv"

if not os.path.exists(SRC):
    sys.exit(f"missing {SRC}")

d = pd.read_csv(SRC, sep="\t")
missed = d[~d["pass"]]
filtered = missed[missed.emitted]
absent = missed[~missed.emitted]

# GATK's tag names are identifiers, not English. On a slide they should
# read as the reasons they are.
READABLE = {
    "weak_evidence":   "too few reads support it",
    "strand_bias":     "all support on one strand",
    "map_qual":        "reads map ambiguously",
    "haplotype":       "conflicts with a nearby variant",
    "normal_artifact": "the normal carries it too",
    "base_qual":       "low base quality",
    "clustered_events": "too many variants close together",
    "slippage":        "sits in a repeat",
    "position":        "only near read ends",
    "contamination":   "attributed to contamination",
    "germline":        "looks germline",
    "panel_of_normals": "seen in a panel of normals",
    "multiallelic":    "more than one alternate here",
    "fragment":        "mate pairs disagree",
    "duplicate":       "duplicate reads only",
    "low_allele_frac": "allele fraction below threshold",
    "n_ratio":         "too many ambiguous bases",
    "strict_strand":   "strand support too uneven",
}

tags = {}
for f in filtered.filters.dropna():
    for t in str(f).split(";"):
        label = READABLE.get(t, t.replace("_", " "))
        tags[label] = tags.get(label, 0) + 1
tags = dict(sorted(tags.items(), key=lambda x: -x[1]))

W = 64
print("=" * W)
print(" WHERE THE MISSES GO")
print("=" * W)
print(f"\n  verified              {len(d)}")
print(f"  recovered             {int(d['pass'].sum())}")
print(f"  missed                {len(missed)}")
print(f"    called, filtered    {len(filtered)}   "
      f"({100*len(filtered)/len(missed):.0f}%)")
print(f"    never called        {len(absent)}   "
      f"({100*len(absent)/len(missed):.0f}%)")
print(f"\n  filter tags:")
for t, k in tags.items():
    print(f"    {t:<22}{k}")
print(f"\n  median supporting reads")
print(f"    recovered           {d[d['pass']].alt_reads.median():.0f}")
print(f"    filtered            {filtered.alt_reads.median():.0f}")
print(f"    never called        {absent.alt_reads.median():.0f}")

os.makedirs(OUT, exist_ok=True)
fig = figure(12, 4.4)
ax = blank(fig, (0.03, 0.05, 0.94, 0.84))

fig.text(0.03, 0.925, f"{len(missed)} SUBSTITUTIONS NOT RECOVERED",
         family=HEAD, fontsize=15, color=BONE, fontweight="bold",
         va="bottom")

# ------------------------------------------------ the split, as one bar
BAR_Y, BAR_H = 0.60, 0.22
frac_filtered = len(filtered) / len(missed)

ax.add_patch(FancyBboxPatch((0, BAR_Y), frac_filtered, BAR_H,
                            boxstyle="round,pad=0,rounding_size=0.02",
                            facecolor=GOLD, edgecolor="none", zorder=3))
ax.add_patch(FancyBboxPatch((frac_filtered + 0.006, BAR_Y),
                            1 - frac_filtered - 0.006, BAR_H,
                            boxstyle="round,pad=0,rounding_size=0.02",
                            facecolor=RUST, edgecolor="none", zorder=3))

ax.text(0.022, BAR_Y + BAR_H / 2, "CALLED, THEN FILTERED", family=HEAD,
        fontsize=13.5, color="#2A1D05", va="center", fontweight="bold",
        zorder=6)
ax.text(frac_filtered - 0.022, BAR_Y + BAR_H / 2,
        f"{len(filtered)}", family=HEAD, fontsize=17, color="#2A1D05",
        ha="right", va="center", fontweight="bold", zorder=6)

ax.text(frac_filtered + 0.028, BAR_Y + BAR_H / 2, "NEVER CALLED",
        family=HEAD, fontsize=13.5, color=BONE, va="center",
        fontweight="bold", zorder=6)
ax.text(0.978, BAR_Y + BAR_H / 2, f"{len(absent)}", family=HEAD,
        fontsize=17, color=BONE, ha="right", va="center",
        fontweight="bold", zorder=6)

ax.text(0.0, BAR_Y - 0.075, "a threshold, and it can be moved",
        family=BODY, fontsize=12, color=GOLD, va="center")
ax.text(frac_filtered + 0.028, BAR_Y - 0.075,
        "a detection limit, and it cannot", family=BODY, fontsize=12,
        color="#C9756A", va="center")

# ------------------------------------------------------- the filter tags
fig.text(0.03, 0.40, "WHY THE FILTERED ONES WERE REJECTED", family=HEAD,
         fontsize=13, color=ASH, fontweight="bold", va="bottom")

top = list(tags.items())[:4]
widest = max(k for _, k in top)
y = 0.28
for t, k in top:
    w = 0.44 * k / widest
    ax.add_patch(FancyBboxPatch((0.42, y), max(w, 0.02), 0.055,
                                boxstyle="round,pad=0,rounding_size=0.012",
                                facecolor=GOLD if k > 50 else SLATE,
                                edgecolor="none", zorder=3))
    ax.text(0.405, y + 0.028, t, family=BODY, fontsize=12.5, color=BONE,
            ha="right", va="center", zorder=6)
    ax.text(0.42 + max(w, 0.02) + 0.014, y + 0.028, str(k), family=BODY,
            fontsize=12.5, color=BONE, va="center", zorder=6)
    y -= 0.085

print("\n" + save(fig, f"{OUT}/mutect_losses.png"))
