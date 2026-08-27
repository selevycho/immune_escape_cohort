#!/usr/bin/env python3
"""
Three figures for the sweep slides.

One — filter stringency. Every setting sees the same variants; they differ
only in what they keep. Loosening the filters returns 222 more true
mutations and not one false call, which turns the missing quarter from a
limit into a choice.

Two — sequencing depth. Both BAMs downsampled and recalled at four levels.
Recall falls almost linearly, about twelve points for every five-fold
step down, and at ten-fold three quarters of the truth set is invisible.

Three — what does not matter. Mutation burden, peptide length beyond the
canonical nine, and the number of IMGT subtypes handed to LOHHLA. Negative
results, collected in one place so they are stated rather than assumed.

Usage:
  python plot_sweeps.py [outdir]
"""
import os
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from deckplot import *
from matplotlib.patches import FancyBboxPatch

WS = os.environ.get("WS", "/gpfs/bwfor/work/ws/fr_os136-immune_escape")
RES = os.path.expanduser("~/immune_escape_project/results")
OUT = sys.argv[1] if len(sys.argv) > 1 else f"{RES}/figures"
os.makedirs(OUT, exist_ok=True)
num = lambda s: pd.to_numeric(s, errors="coerce")

# ======================================================= one: filters
f = pd.read_csv(
    f"{WS}/simulation/step10_filter_sweep/results/filter_sweep_summary.tsv",
    sep="\t")
order = ["strict", "default", "relaxed", "permissive", "verypermissive"]
f["_o"] = f.setting.map({k: i for i, k in enumerate(order)})
f = f.sort_values("_o").reset_index(drop=True)
base = float(f[f.setting == "default"].recovered.iloc[0])
ceiling = float(f.seen.iloc[0])

print("=" * 62)
print(" FILTER STRINGENCY")
print("=" * 62)
print(f"\n  {'setting':<16}{'recovered':>11}{'recall':>9}{'vs default':>12}"
      f"{'false':>8}")
for _, r in f.iterrows():
    print(f"  {r.setting:<16}{int(r.recovered):>11}{r.recall:>8.1f}%"
          f"{int(r.recovered - base):>+12}{int(r.not_in_truth):>8}")

fig = figure(12, 4.6)
ax = panel(fig, 0.065, 0.87, bottom=0.20,
           title="recovered at each filter setting, of 1 497 verified")

xs = np.arange(len(f))
ax.bar(xs, f.recall, width=0.5, zorder=3, edgecolor="none",
       color=[SLATE if s == "default" else
              (EMBER if s == "strict" else GOLD) for s in f.setting])

# the detection ceiling: what the caller emitted, however it was filtered
ax.axhline(ceiling, color=BONE, lw=1.4, ls=(0, (5, 4)), alpha=0.6, zorder=5)
ax.text(-0.55, ceiling + 1.5, f"emitted by the caller  {ceiling:.1f}%",
        family=BODY, fontsize=11.5, color=ASH, ha="left", va="bottom",
        zorder=6)

for x, r in zip(xs, f.itertuples()):
    ax.text(x, r.recall + 2.2, f"{r.recall:.1f}%", family=HEAD,
            fontsize=15, color=BONE, ha="center", va="bottom", zorder=6)
    delta = int(r.recovered - base)
    if delta:
        ax.text(x, 3.5, f"{delta:+d}", family=BODY, fontsize=12,
                color="#2A1D05" if r.setting != "strict" else BONE,
                ha="center", va="bottom", zorder=6)
    ax.text(x, r.recall - 5.5, f"{int(r.not_in_truth)} false",
            family=BODY, fontsize=10.5,
            color="#2A1D05" if r.setting != "strict" else BONE,
            ha="center", va="top", zorder=6)

ax.set_xticks(xs)
ax.set_xticklabels(f.setting, fontsize=11.5)
ax.tick_params(axis="x", colors=ASH, labelsize=11.5, length=0, pad=8)
ax.set_ylim(0, 100)
ax.set_xlim(-0.65, len(f) - 0.35)
ax.yaxis.set_major_formatter(lambda v, p: f"{v:.0f}%")
finish(ax)
print("\n" + save(fig, f"{OUT}/sweep_filters.png"))

# ======================================================== two: depth
c = pd.read_csv(f"{RES}/coverage_sweep_summary.tsv", sep="\t")
c = c.sort_values("depth")

print("\n" + "=" * 62)
print(" SEQUENCING DEPTH")
print("=" * 62)
print(f"\n  {'depth':<9}{'verified':>10}{'found':>8}{'recall':>9}")
for _, r in c.iterrows():
    print(f"  {int(r.depth):<9}{int(r.verified):>10}{int(r.found):>8}"
          f"{r.recall:>8.1f}%")

fig = figure(12, 4.6)

axl = panel(fig, 0.065, 0.375, bottom=0.20,
            title="recall against sequencing depth")

# the full-depth point belongs on the curve: it is the same measurement
xs = list(c.depth) + [34]
ys = list(c.recall) + [74.5]

axl.plot(xs, ys, color=EMBER, lw=2.4, marker="o", ms=8, zorder=4)
axl.scatter([34], [74.5], s=140, facecolor=GOLD, edgecolor=BONE,
            linewidth=1.6, zorder=6)

for x, y in zip(xs, ys):
    axl.text(x, y + 3.2, f"{y:.0f}%", family=HEAD, fontsize=13,
             color=BONE, ha="center", va="bottom", zorder=6)
axl.text(34, 74.5 - 4.5, "as sequenced", family=BODY, fontsize=11,
         color=GOLD, ha="center", va="top", zorder=6)

axl.set_ylim(0, 92)
axl.set_xlim(6, 39)
axl.set_xticks([10, 15, 20, 30, 34])
axl.set_xticklabels(["10x", "15x", "20x", "30x", "34x"])
axl.tick_params(axis="x", colors=ASH, labelsize=12, length=0, pad=8)
axl.yaxis.set_major_formatter(lambda v, p: f"{v:.0f}%")
finish(axl)

# what each halving costs
axr = panel(fig, 0.585, 0.375, bottom=0.20,
            title="mutations lost against full depth")

lost = [1497 - int(r.found) for _, r in c.iterrows()]
xs2 = np.arange(len(c))
axr.bar(xs2, lost, width=0.55, zorder=3, edgecolor="none",
        color=[EMBER if v > 700 else (GOLD if v > 400 else SLATE)
               for v in lost])
top = max(lost)
for x, v in zip(xs2, lost):
    axr.text(x, v + top * 0.04, f"{v}", family=HEAD, fontsize=15,
             color=BONE, ha="center", va="bottom", zorder=6)
axr.set_xticks(xs2)
axr.set_xticklabels([f"{int(v)}x" for v in c.depth])
axr.tick_params(axis="x", colors=ASH, labelsize=12, length=0, pad=8)
axr.set_ylim(0, top * 1.20)
axr.set_xlim(-0.65, len(c) - 0.35)
finish(axr)

print("\n" + save(fig, f"{OUT}/sweep_depth.png"))

# ============================================== three: what does not matter
fig = figure(12, 4.4)
ax = blank(fig, (0.03, 0.05, 0.94, 0.86))
fig.text(0.03, 0.925, "THREE THINGS THAT TURNED OUT NOT TO MATTER",
         family=HEAD, fontsize=14.5, color=BONE, fontweight="bold",
         va="bottom")

d3 = pd.read_csv(f"{RES}/step3_rescored_full.tsv", sep="\t")
g = d3.groupby("sample").agg(n=("pass", "size"), recall=("pass", "mean"),
                             vaf=("observed_vaf", "median")).reset_index()
r_n = float(np.corrcoef(g.n, g.recall)[0, 1])
r_v = float(np.corrcoef(g.vaf, g.recall)[0, 1])

imgt = pd.read_csv(
    f"{WS}/simulation/step16_imgt_sweep/results/imgt_sweep.tsv", sep="\t")
n_cfg = len(imgt)
n_p = len(num(imgt.pval).dropna().unique())

leng = pd.read_csv(f"{RES}/binders_by_length.tsv", sep="\t") \
    if os.path.exists(f"{RES}/binders_by_length.tsv") else None

CARDS = [
    ("MUTATION BURDEN", f"r = {r_n:.2f}",
     f"against r = {r_v:.2f} for allele fraction;\nrecall follows the fractions, not the count"),
    ("IMGT SUBTYPES", f"{n_cfg} configurations",
     f"{n_p} distinct p-values, identical to six decimals;\nLOHHLA uses one sequence per allele name"),
    ("PEPTIDE LENGTH", "9-mers dominate",
     "1 187 strong of 1 830, four times any other length;\nthe rest are included, not discarded"),
]

CW, CH, CY = 0.315, 0.44, 0.28
x = 0.0
for title, big, sub in CARDS:
    ax.add_patch(FancyBboxPatch((x, CY), CW, CH,
                                boxstyle="round,pad=0,rounding_size=0.02",
                                facecolor=INK, edgecolor="#4A3B38", lw=1.4,
                                zorder=3))
    ax.text(x + 0.025, CY + CH - 0.07, title, family=HEAD, fontsize=12.5,
            color=ASH, va="center", fontweight="bold", zorder=6)
    ax.text(x + 0.025, CY + CH - 0.175, big, family=HEAD, fontsize=17,
            color=GOLD, va="center", fontweight="bold", zorder=6)
    ax.text(x + 0.025, CY + 0.11, sub, family=BODY, fontsize=11,
            color=BONE, va="center", linespacing=1.6, zorder=6)
    x += CW + 0.0275

ax.text(0.0, 0.13, "A negative result costs the same to obtain as a "
        "positive one and is worth stating.",
        family=BODY, fontsize=12, color=DUSK, va="center")

print("\n" + "=" * 62)
print(" WHAT DOES NOT MATTER")
print("=" * 62)
print(f"\n  recall vs mutation count   r = {r_n:.3f}")
print(f"  recall vs median VAF       r = {r_v:.3f}")
print(f"  IMGT configurations        {n_cfg}, {n_p} distinct p-values")
if leng is not None:
    print(f"\n  peptide length:")
    print(leng.to_string(index=False))

print("\n" + save(fig, f"{OUT}/sweep_negatives.png"))
