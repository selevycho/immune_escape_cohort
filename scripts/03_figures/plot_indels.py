#!/usr/bin/env python3
"""
The two indel figures, matching the substitution pair.

Indels get their own slides rather than a line on the substitution ones
because they were injected by a different BAMSurgeon entry point, into a
separate BAM, and because their failure mode is different: a substitution
that fails is usually a coverage problem, while every indel that failed
here sits in a repeat where realignment cannot decide which copy a read
belongs to.

The sample is small — 75 events against 1 528 — so the first figure prints
Wilson intervals rather than bare rates. A 100% success on sixteen
insertions is compatible with a true rate near 80%, and saying so is part
of reporting it.

Usage:
  python plot_indels.py [outdir]
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
SRC = f"{RES}/verify_step8_per_indel.tsv"

if not os.path.exists(SRC):
    sys.exit(f"missing {SRC}\nrun verify_step8_indels.py first")

d = pd.read_csv(SRC, sep="\t")
n, landed = len(d), int(d.landed.sum())
lost = d[~d.landed]


def wilson(k, m):
    if m == 0:
        return (np.nan, np.nan)
    z, p = 1.96, k / m
    den = 1 + z * z / m
    c = (p + z * z / (2 * m)) / den
    h = z * np.sqrt(p * (1 - p) / m + z * z / (4 * m * m)) / den
    return (100 * max(0, c - h), 100 * min(1, c + h))


print(f"  indels attempted {n}, landed {landed} ({100*landed/n:.1f}%)")
for ty, g in d.groupby("type"):
    k = int(g.landed.sum())
    lo, hi = wilson(k, len(g))
    print(f"    {ty:<6}{k:>4} of {len(g):<4}  {100*k/len(g):>5.1f}%   "
          f"CI {lo:.0f}-{hi:.0f}%")

# ============================================================ figure one
fig = figure(12, 4.4)

TYPES = [t for t in ["INS", "DEL"] if (d.type == t).any()]
rates, ns, errs = [], [], []
for ty in TYPES:
    g = d[d.type == ty]
    k = int(g.landed.sum())
    rates.append(100 * k / len(g))
    ns.append((k, len(g)))
    lo, hi = wilson(k, len(g))
    errs.append([rates[-1] - lo, hi - rates[-1]])

ax = panel(fig, 0.07, 0.36, bottom=0.20, title="confirmed in the reads")
xs = np.arange(len(TYPES))
ax.bar(xs, rates, width=0.5, zorder=3, edgecolor="none",
       color=[GREEN if r >= 95 else EMBER for r in rates])
ax.errorbar(xs, rates, yerr=np.array(errs).T, fmt="none", ecolor=GOLD,
            elinewidth=1.5, capsize=7, capthick=1.5, zorder=5)
for x, (k, tot) in zip(xs, ns):
    ax.text(x, 4, f"{k} of {tot}", family=BODY, fontsize=12.5, color=BONE,
            ha="center", va="bottom", zorder=6)
ax.set_xticks(xs)
ax.set_xticklabels(["insertions", "deletions"])
ax.tick_params(axis="x", colors=ASH, labelsize=12.5, length=0, pad=8)
ax.set_ylim(0, 120)
ax.yaxis.set_major_formatter(lambda v, p: f"{v:.0f}%")
finish(ax)
note(ax, "bars show 95% intervals", loc="upper right", color=DUSK, size=11)

# the genes where every failure sits
ax2 = fig.add_axes([0.56, 0.20, 0.40, 0.63])
ax2.axis("off")
ax2.set_xlim(0, 1)
ax2.set_ylim(0, 1)
fig.text(0.56, 0.885, "WHERE THE FAILURES SIT", family=HEAD, fontsize=14,
         color=BONE, fontweight="bold", va="bottom")

if len(lost) and "gene" in lost.columns:
    gc = lost.gene.value_counts()
    y, step = 0.76, 0.26
    for g, k in gc.items():
        ax2.add_patch(FancyBboxPatch((0, y - 0.02), 0.34, 0.19,
                                     boxstyle="round,pad=0,rounding_size=0.03",
                                     facecolor="#3A1A13", edgecolor=EMBER,
                                     lw=1.4, zorder=3))
        ax2.text(0.17, y + 0.075, g, family=HEAD, fontsize=15, color=BONE,
                 ha="center", va="center", fontweight="bold", zorder=6)
        ax2.text(0.39, y + 0.075, f"{k} lost", family=BODY, fontsize=13,
                 color=GOLD, va="center")
        y -= step
    ax2.text(0, 0.03, "segmental duplications and tandem repeats",
             family=BODY, fontsize=12, color=DUSK, va="center")
print(save(fig, f"{OUT}/indel_landing.png"))

# ============================================================ figure two
ok = d[d.landed & (d.target_vaf > 0) & (d.observed_vaf > 0)]
r = np.corrcoef(ok.target_vaf, ok.observed_vaf)[0, 1]
ratio = float((ok.observed_vaf / ok.target_vaf).median())
print(f"  r = {r:.3f}, observed/requested median {ratio:.3f}")

fig = figure(12, 4.4)
ax = panel(fig, 0.07, 0.38, bottom=0.18, title="allele fraction")
for ty, col in [("DEL", EMBER), ("INS", GOLD)]:
    g = ok[ok.type == ty]
    if len(g):
        ax.scatter(g.target_vaf, g.observed_vaf, s=34, color=col, alpha=0.75,
                   edgecolors="none", zorder=3, label=ty)
lim = max(ok.target_vaf.max(), ok.observed_vaf.max()) * 1.08
ax.plot([0, lim], [0, lim], color=BONE, lw=1.1, ls=(0, (5, 4)), zorder=5,
        alpha=0.5)
ax.set_xlim(0, lim)
ax.set_ylim(0, lim)
note(ax, f"r = {r:.3f}", loc="upper left")
leg = ax.legend(frameon=False, loc="lower right", fontsize=11.5)
for t in leg.get_texts():
    t.set_color(ASH)
    t.set_fontfamily(BODY)
finish(ax, "requested", "observed in the reads", grid="both")

# the same shortfall the substitutions show, side by side
ax2 = panel(fig, 0.58, 0.36, bottom=0.18, title="observed / requested")
labels = ["indels", "substitutions"]
vals = [ratio, 0.83]
xs = np.arange(2)
ax2.bar(xs, vals, width=0.5, color=[EMBER, SLATE], zorder=3, edgecolor="none")
ax2.axhline(1.0, color=GOLD, lw=1.3, ls=(0, (5, 4)), zorder=5)
for x, v in zip(xs, vals):
    ax2.text(x, v + 0.03, f"{v:.2f}", family=BODY, fontsize=13, color=BONE,
             ha="center", va="bottom", zorder=6)
ax2.set_xticks(xs)
ax2.set_xticklabels(labels)
ax2.tick_params(axis="x", colors=ASH, labelsize=12.5, length=0, pad=8)
ax2.set_ylim(0, 1.22)
ax2.set_xlim(-0.7, 1.7)
finish(ax2)
note(ax2, "1.00 = exactly as asked", loc="upper right", color=DUSK, size=11)
print(save(fig, f"{OUT}/indel_vaf.png"))
