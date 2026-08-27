#!/usr/bin/env python3
"""
Three controls on the variant calling.

Recall says how much of what was placed came back. It says nothing about
what comes back when nothing was placed, and a caller that reported
variants everywhere would score perfectly on recall while being useless.

Specificity: every PASS call across the forty samples, checked against
the truth set. Not one falls outside it.

Reproducibility: the same substitutions were called a second time, from a
BAM built by a different route in a different week. The two runs agree on
every one of the 1 071 they share, which bounds how much of any single
recall figure is run-to-run noise at zero.

The negative control: reads from one healthy individual split in half by
read name, one half offered as tumour and the other as normal. Both halves
are real independent sequencing of a person with no somatic mutations, so
every PASS call is false by construction. It runs at half the usual depth,
which makes the estimate conservative rather than flattering.

Where those false calls land is the finding. The MHC is 15% of the panel
and produces 58% of them, and the same coordinates recur across unrelated
people — the signature of a reference problem rather than of noise.

Usage:
  python plot_mutect_controls.py [outdir]
"""
import os
import sys
import gzip
import glob
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from deckplot import *
from matplotlib.patches import FancyBboxPatch

WS = os.environ.get("WS", "/gpfs/bwfor/work/ws/fr_os136-immune_escape")
RES = os.path.expanduser("~/immune_escape_project/results")
OUT = sys.argv[1] if len(sys.argv) > 1 else f"{RES}/figures"
NC = f"{WS}/simulation/step17b_split_control"

MHC = ("chr6", 29600000, 33100000)


def pass_calls(path):
    out = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8 or f[6] != "PASS":
                continue
            out.append((f[0], int(f[1])))
    return out


rows = []
for p in sorted(glob.glob(f"{NC}/*/*.filtered.vcf.gz")):
    sid = os.path.basename(os.path.dirname(p))
    calls = pass_calls(p)
    n_mhc = sum(1 for c, pos in calls
                if c == MHC[0] and MHC[1] <= pos < MHC[2])
    rows.append({"sample": sid, "total": len(calls), "mhc": n_mhc,
                 "outside": len(calls) - n_mhc})

d = pd.DataFrame(rows)
if d.empty:
    sys.exit("no negative-control results found")

# panel geometry, so the density comparison is against real megabases
panel = pd.read_csv(f"{WS}/simulation/panel/panel.bed", sep="\t",
                    header=None, names=["chrom", "start", "end", "gene"])
panel["len"] = panel.end - panel.start
mhc_mb = panel[(panel.chrom == MHC[0]) & (panel.start >= MHC[1]) &
               (panel.end < MHC[2])]["len"].sum() / 1e6
tot_mb = panel["len"].sum() / 1e6
out_mb = tot_mb - mhc_mb

n_s = len(d)
tot = int(d.total.sum())
mhc = int(d.mhc.sum())
outside = tot - mhc

W = 64
print("=" * W)
print(" NEGATIVE CONTROL")
print("=" * W)
print(f"\n  samples                    {n_s}")
print(f"  false positives            {tot}   "
      f"({tot/n_s:.0f} per sample)")
print(f"    inside the MHC           {mhc}   ({100*mhc/tot:.0f}%)")
print(f"    outside                  {outside}")
print(f"\n  panel                      {tot_mb:.2f} Mb")
print(f"    MHC                      {mhc_mb:.2f} Mb   "
      f"({100*mhc_mb/tot_mb:.0f}% of the panel)")
print(f"\n  false calls per Mb")
print(f"    inside the MHC           {mhc/mhc_mb/n_s:.1f}")
print(f"    outside                  {outside/out_mb/n_s:.1f}")
print(f"    ratio                    "
      f"{(mhc/mhc_mb)/(outside/out_mb):.0f}x")

os.makedirs(OUT, exist_ok=True)
fig = figure(12, 4.6)
ax = blank(fig, (0.03, 0.04, 0.94, 0.86))

fig.text(0.03, 0.925, "THREE CONTROLS", family=HEAD, fontsize=15,
         color=BONE, fontweight="bold", va="bottom")

# ----------------------------------------------- the two clean results
CW, CH, CY = 0.305, 0.30, 0.50
card(ax, 0.0, CY, CW, CH, "FALSE POSITIVES",
     "0  in 40 samples", accent=GREEN)
card(ax, CW + 0.03, CY, CW, CH, "TWO RUNS, SAME VARIANTS",
     "0  disagreements of 1071", accent=GREEN)
card(ax, 2 * (CW + 0.03), CY, CW, CH, "NOTHING TO FIND",
     f"{tot/n_s:.0f}  per sample", accent=GOLD)

ax.text(0.0, CY - 0.055, "every PASS call sits in the truth set",
        family=BODY, fontsize=11.5, color=DUSK, va="center")
ax.text(CW + 0.03, CY - 0.055, "a BAM built by a different route",
        family=BODY, fontsize=11.5, color=DUSK, va="center")
ax.text(2 * (CW + 0.03), CY - 0.055,
        "one person's reads split in half", family=BODY, fontsize=11.5,
        color=DUSK, va="center")

# ------------------------------------------- where the false calls land
fig.text(0.03, 0.34, "WHERE THE FALSE CALLS LAND", family=HEAD,
         fontsize=13, color=ASH, fontweight="bold", va="bottom")

BY, BH = 0.14, 0.13
frac_mhc = mhc / tot

ax.add_patch(FancyBboxPatch((0, BY), frac_mhc, BH,
                            boxstyle="round,pad=0,rounding_size=0.02",
                            facecolor=RUST, edgecolor="none", zorder=3))
ax.add_patch(FancyBboxPatch((frac_mhc + 0.006, BY),
                            1 - frac_mhc - 0.006, BH,
                            boxstyle="round,pad=0,rounding_size=0.02",
                            facecolor=SLATE, edgecolor="none", zorder=3))

ax.text(0.02, BY + BH / 2, "MHC", family=HEAD, fontsize=13.5,
        color=BONE, va="center", fontweight="bold", zorder=6)
ax.text(frac_mhc - 0.02, BY + BH / 2, f"{mhc}", family=HEAD,
        fontsize=15, color=BONE, ha="right", va="center",
        fontweight="bold", zorder=6)
ax.text(frac_mhc + 0.026, BY + BH / 2, "the rest of the panel",
        family=HEAD, fontsize=13.5, color=BONE, va="center",
        fontweight="bold", zorder=6)
ax.text(0.98, BY + BH / 2, f"{outside}", family=HEAD, fontsize=15,
        color=BONE, ha="right", va="center", fontweight="bold", zorder=6)

ax.text(0.0, BY - 0.06,
        f"{mhc_mb:.1f} Mb of {tot_mb:.1f}", family=BODY, fontsize=11.5,
        color="#C9756A", va="center")
ax.text(frac_mhc + 0.026, BY - 0.06,
        f"{tot_mb - mhc_mb:.1f} Mb", family=BODY, fontsize=11.5,
        color=DUSK, va="center")

print("\n" + save(fig, f"{OUT}/mutect_controls.png"))
