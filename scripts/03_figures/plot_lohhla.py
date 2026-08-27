#!/usr/bin/env python3
"""
Two figures for the LOHHLA slides.

The first is the result: how often the method detects a loss that was
made, and how often it reports one that was not. Sensitivity rises with
the amount removed, and specificity is absolute — not one control locus
reports a loss at any threshold tested.

The second is the reason sensitivity is what it is. LOHHLA does not fail
for want of reads; it fails for want of positions where the two alleles
differ and both are covered. 338 of 376 attempts never reached the test
at all, and the two commonest reasons say why.

Usage:
  python plot_lohhla.py [outdir]
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
S9 = f"{WS}/simulation/step9_lohhla_panel/results"

num = lambda s: pd.to_numeric(s, errors="coerce")

d = pd.read_csv(f"{S9}/step9_all_loci_mc5.tsv", sep="\t")
d["pval"] = num(d.pval)
roc = pd.read_csv(f"{S9}/lohhla_roc_mc5.tsv", sep="\t")

os.makedirs(OUT, exist_ok=True)

# ==================================================== figure one: result
fig = figure(12, 4.6)

# ------------------------------- left: detections by how much was removed
axl = panel(fig, 0.065, 0.375, bottom=0.20,
            title="loci reporting a loss, at p < 0.05")

order = ["A_untouched", "B_loss20", "C_loss35", "D_loss50"]
labels = {"A_untouched": "nothing\nremoved", "B_loss20": "kept 20%",
          "C_loss35": "kept 35%", "D_loss50": "kept 50%"}
rows = []
for c in order:
    g = d[d.collection == c]
    if not len(g):
        continue
    sg = g[g.pval < 0.05]
    rows.append({"c": c, "n": len(g),
                 "got": int(g.pval.notna().sum()),
                 "sig": len(sg),
                 "on_target": int((sg.role == "TARGET").sum())
                 if "role" in g else 0})
t = pd.DataFrame(rows)

print("=" * 62)
print(" DETECTIONS BY AMOUNT REMOVED")
print("=" * 62)
print(f"\n  {'':<14}{'loci':>6}{'tested':>8}{'p<0.05':>8}{'on target':>11}")
for _, r in t.iterrows():
    print(f"  {r.c:<14}{int(r.n):>6}{int(r.got):>8}{int(r.sig):>8}"
          f"{int(r.on_target):>11}")

xs = np.arange(len(t))
axl.bar(xs, t.sig, width=0.55, zorder=3, edgecolor="none",
        color=[SLATE if v == 0 else EMBER for v in t.sig])

top = max(t.sig.max(), 1)
for x, r in zip(xs, t.itertuples()):
    axl.text(x, r.sig + top * 0.06, f"{int(r.sig)}", family=HEAD,
             fontsize=16, color=BONE, ha="center", va="bottom", zorder=6)
    if r.sig:
        axl.text(x, top * 0.05, "all on target", family=BODY, fontsize=10.5,
                 color="#2A1008", ha="center", va="bottom", zorder=6)

axl.set_xticks(xs)
axl.set_xticklabels([labels[c] for c in t.c], fontsize=11)
axl.tick_params(axis="x", colors=ASH, labelsize=11, length=0, pad=8)
axl.set_ylim(0, top * 1.32)
axl.set_xlim(-0.65, len(t) - 0.35)
finish(axl)

# ---------------------------- right: sensitivity against the threshold
axr = panel(fig, 0.585, 0.375, bottom=0.20,
            title="sensitivity, and false positives")

axr.plot(roc.threshold, roc.sensitivity, color=EMBER, lw=2.2,
         marker="o", ms=6, zorder=4)
axr.plot(roc.threshold, roc.fp, color=GREEN, lw=2.2, marker="s", ms=6,
         zorder=4)

axr.set_xscale("log")
axr.set_ylim(-4, 100)
axr.yaxis.set_major_formatter(lambda v, p: f"{v:.0f}%")

# label the two lines where they are furthest apart rather than in a box
axr.text(roc.threshold.iloc[-1] * 0.9, roc.sensitivity.iloc[-1] - 8,
         "detected", family=BODY, fontsize=12, color=EMBER,
         ha="right", va="top", zorder=6)
axr.text(roc.threshold.iloc[-1] * 0.9, 5, "false positives",
         family=BODY, fontsize=12, color=GREEN, ha="right",
         va="bottom", zorder=6)

finish(axr, "p-value threshold")

print(f"\n  false positives at every threshold: {int(roc.fp.sum())}")
print("\n" + save(fig, f"{OUT}/lohhla_result.png"))

# ================================================ figure two: why it fails
fig = figure(12, 4.4)
ax = blank(fig, (0.03, 0.05, 0.94, 0.84))

n_total = len(d)
n_tested = int(d.pval.notna().sum())
n_failed = n_total - n_tested

fig.text(0.03, 0.925, f"{n_total} LOCUS ATTEMPTS", family=HEAD,
         fontsize=15, color=BONE, fontweight="bold", va="bottom")

BY, BH = 0.62, 0.20
frac = n_failed / n_total
ax.add_patch(FancyBboxPatch((0, BY), frac, BH,
                            boxstyle="round,pad=0,rounding_size=0.02",
                            facecolor=RUST, edgecolor="none", zorder=3))
ax.add_patch(FancyBboxPatch((frac + 0.006, BY), 1 - frac - 0.006, BH,
                            boxstyle="round,pad=0,rounding_size=0.02",
                            facecolor=GOLD, edgecolor="none", zorder=3))

ax.text(0.02, BY + BH / 2, "NEVER REACHED THE TEST", family=HEAD,
        fontsize=13, color=BONE, va="center", fontweight="bold", zorder=6)
ax.text(frac - 0.02, BY + BH / 2, f"{n_failed}", family=HEAD, fontsize=17,
        color=BONE, ha="right", va="center", fontweight="bold", zorder=6)
ax.text(0.98, BY + BH / 2, f"{n_tested}", family=HEAD, fontsize=17,
        color="#2A1D05", ha="right", va="center", fontweight="bold",
        zorder=6)

# the reasons, as bars
fig.text(0.03, 0.46, "WHY", family=HEAD, fontsize=13, color=ASH,
         fontweight="bold", va="bottom")

reasons = d[d.reason.notna()].reason.str.slice(0, 42).value_counts().head(4)
widest = reasons.iloc[0]
y = 0.34
for text, k in reasons.items():
    w = 0.42 * k / widest
    ax.add_patch(FancyBboxPatch((0.50, y), max(w, 0.02), 0.06,
                                boxstyle="round,pad=0,rounding_size=0.012",
                                facecolor=RUST if k > 100 else SLATE,
                                edgecolor="none", zorder=3))
    ax.text(0.485, y + 0.03, text, family=BODY, fontsize=12, color=BONE,
            ha="right", va="center", zorder=6)
    ax.text(0.50 + max(w, 0.02) + 0.014, y + 0.03, str(k), family=BODY,
            fontsize=12.5, color=BONE, va="center", zorder=6)
    y -= 0.095

print("\n" + "=" * 62)
print(" WHY LOCI PRODUCE NOTHING")
print("=" * 62)
print(f"\n  {n_failed} of {n_total} never reached the test\n")
for text, k in reasons.items():
    print(f"  {k:>5}   {text}")

r = pd.read_csv(f"{RES}/lohhla_report.tsv", sep="\t")
if "sites" in r.columns:
    s = num(r.sites).dropna()
    s = s[s > 0]
    print(f"\n  where the test did run, distinguishing positions:")
    print(f"    median {s.median():.0f}, range {s.min():.0f}-{s.max():.0f}")

print("\n" + save(fig, f"{OUT}/lohhla_why.png"))
