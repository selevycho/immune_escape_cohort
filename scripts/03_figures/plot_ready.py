#!/usr/bin/env python3
"""
What exists on disk before the first mutation is placed.

A transitional figure: three objects, each with the number that was
verified for it earlier in this section. It carries no new measurement —
its job is to close the data-building half of the talk and make the
starting state of the pipeline concrete before step 2 begins.

Values are read from the verification tables rather than typed in, so the
slide cannot drift away from what the checks reported.

Usage:
  python plot_ready.py [outdir]
"""
import os
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.patches import FancyBboxPatch

RES = os.path.expanduser("~/immune_escape_project/results")
WS = os.environ.get("WS", "/gpfs/bwfor/work/ws/fr_os136-immune_escape")
OUT = sys.argv[1] if len(sys.argv) > 1 else f"{RES}/figures"

EMBER, GOLD = "#E8402A", "#F2A623"
BONE, ASH, DUSK = "#F4EFEA", "#A29591", "#7A6C68"

have = {f.name for f in fm.fontManager.ttflist}
HEAD = next((f for f in ["TeX Gyre Heros Cn", "Liberation Sans Narrow",
                         "Liberation Sans", "DejaVu Sans"] if f in have),
            "DejaVu Sans")
BODY = next((f for f in ["Carlito", "Open Sans", "Liberation Sans",
                         "DejaVu Sans"] if f in have), "DejaVu Sans")

v1 = pd.read_csv(f"{RES}/verify_step1.tsv", sep="\t")
audit = pd.read_csv(f"{RES}/panel_audit.tsv", sep="\t")
panel = pd.read_csv(f"{WS}/simulation/panel/panel.bed", sep="\t", header=None,
                    usecols=[1, 2], names=["s", "e"])

n_bam = len(v1)
gb = v1.size_MB.sum() / 1024
depth = v1.panel_mean_depth.median()
mb = (panel.e - panel.s).sum() / 1e6
n_iv = len(panel)
n_mut = int(audit.injected.sum())
n_genes = int(audit.genes_hit.median())

print("=" * 58)
print(" STATE BEFORE STEP 2")
print("=" * 58)
print(f"\n  normal BAMs           {n_bam}, {gb:.1f} GB")
print(f"  median panel depth    {depth:.1f}x")
print(f"  panel                 {mb:.2f} Mb in {n_iv} intervals")
print(f"  mutations ready       {n_mut}")
print(f"  genes hit, median     {n_genes} per patient")

os.makedirs(OUT, exist_ok=True)
fig = plt.figure(figsize=(12, 2.9), dpi=300)
ax = fig.add_axes([0, 0, 1, 1])
W, H = 1000, 1000 * 2.9 / 12
ax.set_xlim(0, W)
ax.set_ylim(0, H)
ax.axis("off")
fig.patch.set_alpha(0)

CARDS = [
    ("Normal BAMs", f"{n_bam}", f"{gb:.1f} GB  \u00b7  median {depth:.1f}\u00d7",
     "#241A18", "#4A3B38", BONE, GOLD),
    ("Panel", f"{mb:.2f} Mb", f"{n_iv:,} intervals".replace(",", " "),
     "#241A18", "#4A3B38", BONE, GOLD),
    ("Mutations waiting", f"{n_mut:,}".replace(",", " "),
     f"median {n_genes} genes per patient",
     "#3A1A13", EMBER, BONE, GOLD),
]

CW, GAP = 300, 20
X0 = (W - (CW * 3 + GAP * 2)) / 2
CH = 140
CY = (H - CH) / 2

for i, (title, big, sub, fill, edge, tc, vc) in enumerate(CARDS):
    x = X0 + i * (CW + GAP)
    lw = 1.6 if i == 2 else 1.2
    if i == 2:
        for pad, a in ((8, 0.06), (4.5, 0.11), (2, 0.18)):
            ax.add_patch(FancyBboxPatch(
                (x - pad, CY - pad), CW + 2 * pad, CH + 2 * pad,
                boxstyle=f"round,pad=0,rounding_size={10 + pad}",
                facecolor="none", edgecolor=EMBER, lw=1.5, alpha=a, zorder=2))
    ax.add_patch(FancyBboxPatch(
        (x, CY), CW, CH, boxstyle="round,pad=0,rounding_size=10",
        facecolor=fill, edgecolor=edge, lw=lw, zorder=3))
    ax.text(x + CW / 2, CY + CH - 30, title.upper(), family=HEAD,
            fontsize=13.5, color=ASH, ha="center", va="center",
            fontweight="bold", zorder=6)
    ax.text(x + CW / 2, CY + CH - 74, big, family=HEAD, fontsize=32,
            color=vc, ha="center", va="center", fontweight="bold", zorder=6)
    ax.text(x + CW / 2, CY + CH - 108, sub, family=BODY, fontsize=12,
            color=DUSK, ha="center", va="center", zorder=6)

p = f"{OUT}/ready_to_inject.png"
fig.savefig(p, dpi=300, transparent=True, bbox_inches="tight", pad_inches=0.14)
plt.close(fig)
print(f"\nwritten to {p}")
