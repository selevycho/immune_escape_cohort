#!/usr/bin/env python3
"""
Indels, against the substitutions that share the file with them.

Both kinds went into the same BAM and came out of the same Mutect2 run, so
the comparison is between two kinds of mutation rather than between two
experiments. Indels come out ahead: 82.7% against 74.5%.

That ordering surprises people who expect indels to be the harder case. A
frameshift is several bases of disagreement with the reference in a row,
which is a stronger signal than a single substituted base — a substitution
has to be told apart from a sequencing error at that one position, and an
indel does not.

The thirteen that were missed divide the same way the substitutions did,
into filtered and never-called, and seven of them sit in homopolymers
where an inserted base is indistinguishable from polymerase slippage.

Numbers are small — 23 indels against 1 497 substitutions — so the bars
carry Wilson intervals. Two of two insertions landing is consistent with a
true rate anywhere above about 30%, and pretending otherwise would be
dishonest at this sample size.

Usage:
  python plot_indel_recall.py [outdir]
"""
import os
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from deckplot import *

RES = os.path.expanduser("~/immune_escape_project/results")
OUT = sys.argv[1] if len(sys.argv) > 1 else f"{RES}/figures"
IND = f"{RES}/indel_recall.tsv"
SNV = f"{RES}/step3_rescored_full.tsv"

for p in (IND, SNV):
    if not os.path.exists(p):
        sys.exit(f"missing {p}")

i = pd.read_csv(IND, sep="\t")
s = pd.read_csv(SNV, sep="\t")


def wilson(k, n):
    if n == 0:
        return (np.nan, np.nan)
    z, p = 1.96, k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return (100 * max(0, c - h), 100 * min(1, c + h))


groups = [("substitutions", int(s["pass"].sum()), len(s))]
for ty in ("DEL", "INS"):
    g = i[i.type == ty]
    if len(g):
        label = "deletions" if ty == "DEL" else "insertions"
        groups.append((label, int(g.found.sum()), len(g)))
groups.append(("indels, both", int(i.found.sum()), len(i)))

W = 62
print("=" * W)
print(" RECOVERY BY MUTATION TYPE")
print("=" * W)
print(f"\n  {'':<18}{'found':>8}{'of':>6}{'rate':>9}{'95% CI':>16}")
for label, k, n in groups:
    lo, hi = wilson(k, n)
    print(f"  {label:<18}{k:>8}{n:>6}{100*k/n:>8.1f}%"
          f"{f'{lo:.0f}-{hi:.0f}%':>16}")

miss = i[~i.found]
print(f"\n  missed {len(miss)} of {len(i)}")
if len(miss):
    print(f"    median supporting reads   {miss.support.median():.0f}")
    print(f"    median allele fraction    {miss.observed_vaf.median():.3f}")

os.makedirs(OUT, exist_ok=True)
fig = figure(12, 4.4)

# ------------------------------------------------- left: the comparison
ax = panel(fig, 0.065, 0.36, bottom=0.21, title="recovered")
xs = np.arange(len(groups))
vals = [100 * k / n for _, k, n in groups]
errs = []
for _, k, n in groups:
    lo, hi = wilson(k, n)
    errs.append([100 * k / n - lo, hi - 100 * k / n])

ax.bar(xs, vals, width=0.55, zorder=3, edgecolor="none",
       color=[SLATE] + [EMBER] * (len(groups) - 2) + [GOLD])

# capped error bars read as letters at this size; a plain vertical line
# says the same thing without competing with the numbers
ax.errorbar(xs, vals, yerr=np.array(errs).T, fmt="none", ecolor=BONE,
            elinewidth=1.6, capsize=0, alpha=0.55, zorder=5)

# each label sits above its own interval, not above the tallest one —
# the insertions bar carries a 93% upper bound and everything placed at a
# shared height landed inside it
for x, (label, k, n), v, e in zip(xs, groups, vals, errs):
    ax.text(x, v + e[1] + 3.5, f"{v:.1f}%",
            family=HEAD, fontsize=14, color=BONE, ha="center",
            va="bottom", zorder=6)
    ax.text(x, 4, f"{k} of {n}", family=BODY, fontsize=11, color=BONE,
            ha="center", va="bottom", zorder=6)

ax.set_xticks(xs)
ax.set_xticklabels([g[0].replace(", both", "") for g in groups],
                   fontsize=11)
ax.tick_params(axis="x", colors=ASH, labelsize=11, length=0, pad=8)
ax.set_ylim(0, 108)
ax.set_xlim(-0.65, len(groups) - 0.35)
ax.yaxis.set_major_formatter(lambda v, p: f"{v:.0f}%")
finish(ax)

# ------------------------------- right: detection against read support
ax2 = panel(fig, 0.575, 0.375, bottom=0.21,
            title="found and missed, by supporting reads")

found = i[i.found]
lost = i[~i.found]

rng = np.random.default_rng(7)
if len(found):
    ax2.scatter(found.support + rng.uniform(-.15, .15, len(found)),
                found.observed_vaf, s=52, color=GREEN, alpha=0.85,
                edgecolors="none", zorder=4, label="found")
if len(lost):
    ax2.scatter(lost.support + rng.uniform(-.15, .15, len(lost)),
                lost.observed_vaf, s=52, color=EMBER, alpha=0.9,
                edgecolors="none", zorder=5, label="missed")

leg = ax2.legend(frameon=False, loc="lower right", fontsize=11.5,
                 markerscale=1.1)
for t in leg.get_texts():
    t.set_color(ASH)
    t.set_fontfamily(BODY)

finish(ax2, "reads carrying the indel", "allele fraction", grid="both")

print("\n" + save(fig, f"{OUT}/indel_recall.png"))
