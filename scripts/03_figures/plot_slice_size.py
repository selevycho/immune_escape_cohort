#!/usr/bin/env python3
"""
What remote slicing moves, drawn as areas.

Two earlier attempts failed for the same reason: they compared the panel
to an exome target in megabases, which is 40% and looks like no saving at
all. The saving is not there. It is that the CRAM on the server is
fifteen gigabytes and what crosses the network is a few hundred megabytes
— the panel is a thin fraction of the file, not of the capture design.

So the figure compares file sizes by area. The block is the CRAM; the
stripes inside it are what gets fetched, drawn at their true share; the
block below is the BAM that results, at the same scale. The proportion is
the argument and the two numbers only name it.

CRAM sizes are read from the EBI server by HEAD request, so the large
number is measured rather than assumed. If the network is unavailable the
script falls back to a stated constant and says so.

Usage:
  python plot_slice_size.py [outdir] [--nofetch]
"""
import os
import sys
import subprocess
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

RES = os.path.expanduser("~/immune_escape_project/results")
V1 = f"{RES}/verify_step1.tsv"
WS = os.environ.get("WS", "/gpfs/bwfor/work/ws/fr_os136-immune_escape")
MANIFEST = f"{WS}/simulation/cohort/manifest.tsv"
args = [a for a in sys.argv[1:] if not a.startswith("--")]
OUT = args[0] if args else f"{RES}/figures"
FETCH = "--nofetch" not in sys.argv

EMBER, GOLD = "#E8402A", "#F2A623"
BONE, ASH, DUSK = "#F4EFEA", "#A29591", "#7A6C68"

have = {f.name for f in fm.fontManager.ttflist}
HEAD = next((f for f in ["TeX Gyre Heros Cn", "Liberation Sans Narrow",
                         "Liberation Sans", "DejaVu Sans"] if f in have),
            "DejaVu Sans")
BODY = next((f for f in ["Carlito", "Open Sans", "Liberation Sans",
                         "DejaVu Sans"] if f in have), "DejaVu Sans")

d = pd.read_csv(V1, sep="\t")
bam_mb = d.size_MB.median()

cram_gb, n_probed = 15.0, 0
if FETCH and os.path.exists(MANIFEST):
    m = pd.read_csv(MANIFEST, sep="\t")
    col = next((c for c in m.columns if "cram" in c.lower()), None)
    sizes = []
    for url in m[col].dropna().head(8):
        if not str(url).startswith("http"):
            continue
        r = subprocess.run(
            f"curl -sIL --max-time 20 '{url}' | "
            f"awk 'tolower($1)==\"content-length:\" {{v=$2}} END {{print v+0}}'",
            shell=True, capture_output=True, text=True)
        try:
            v = int(r.stdout.strip())
            if v > 0:
                sizes.append(v / 1024 ** 3)
        except ValueError:
            pass
    if sizes:
        cram_gb = float(np.median(sizes))
        n_probed = len(sizes)

frac = (bam_mb / 1024) / cram_gb

print("=" * 62)
print(" WHAT CROSSES THE NETWORK")
print("=" * 62)
print(f"\n  CRAM on the server     {cram_gb:.1f} GB"
      f"{f'  (median of {n_probed} probed)' if n_probed else '  (assumed)'}")
print(f"  BAM that results       {bam_mb:.0f} MB, "
      f"range {d.size_MB.min():.0f}-{d.size_MB.max():.0f}")
print(f"  fraction transferred   {100 * frac:.1f}%")
print(f"  forty normal BAMs      {d.size_MB.sum() / 1024:.1f} GB")
print(f"  forty whole CRAMs      {cram_gb * 40:.0f} GB")

os.makedirs(OUT, exist_ok=True)
fig = plt.figure(figsize=(12, 5.0), dpi=300)
ax = fig.add_axes([0, 0, 1, 1])
W, H = 1000, 1000 * 5.0 / 12
ax.set_xlim(0, W)
ax.set_ylim(0, H)
ax.axis("off")
fig.patch.set_alpha(0)


def rr(x, y, w, h, fc, ec=None, r=6, z=3, lw=1.2, a=1.0):
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                                boxstyle=f"round,pad=0,rounding_size={r}",
                                facecolor=fc, edgecolor=ec or fc, lw=lw,
                                zorder=z, alpha=a))


# the CRAM: a wide block whose area stands for its size
BX, BW, BH = 60, 880, 150
BY = H - 40 - BH
rr(BX, BY, BW, BH, "#241A18", "#4A3B38", r=10, lw=1.3)
ax.text(BX + 26, BY + BH - 30, "CRAM ON THE EBI SERVER", family=HEAD,
        fontsize=15, color=BONE, va="center", weight="bold", zorder=7)
ax.text(BX + BW - 26, BY + BH - 30, f"{cram_gb:.1f} GB", family=HEAD,
        fontsize=24, color=ASH, ha="right", va="center", weight="bold",
        zorder=7)

# the fetched byte ranges, at their true share of that area
rng = np.random.default_rng(5)
n_seg = 26
seg_w = BW * frac / n_seg
for c in np.sort(rng.uniform(0.03, 0.97, n_seg)):
    rr(BX + BW * c - seg_w / 2, BY + 20, max(seg_w, 1.6), BH - 78,
       EMBER, r=1, z=5)

AY = BY - 20
ax.add_patch(FancyArrowPatch((BX + BW / 2, AY), (BX + BW / 2, AY - 52),
                             arrowstyle="-|>", mutation_scale=18,
                             color=EMBER, lw=2.0, zorder=6,
                             shrinkA=0, shrinkB=0))
ax.text(BX + BW / 2 + 18, AY - 26, "only these byte ranges cross the network",
        family=BODY, fontsize=13, color=GOLD, va="center", zorder=7)

# the BAM: same scale, so its width is the fraction transferred
SH = 90
SY = AY - 52 - 16 - SH
SW = max(BW * frac, 118)
rr(BX, SY, SW, SH, "#3A1A13", EMBER, r=8, lw=1.8, z=4)
ax.text(BX + SW + 26, SY + SH - 32, f"{bam_mb:.0f} MB", family=HEAD,
        fontsize=24, color=BONE, va="center", weight="bold", zorder=7)
ax.text(BX + SW + 26, SY + SH - 62, "BAM ON DISK", family=HEAD,
        fontsize=15, color=ASH, va="center", weight="bold", zorder=7)
ax.text(BX + SW + 190, SY + SH - 32, f"{100 * frac:.1f}%", family=HEAD,
        fontsize=24, color=GOLD, va="center", weight="bold", zorder=7)
ax.text(BX + SW + 190, SY + SH - 62, "OF THE FILE", family=HEAD,
        fontsize=15, color=ASH, va="center", weight="bold", zorder=7)

p = f"{OUT}/slice_size.png"
fig.savefig(p, dpi=300, transparent=True, bbox_inches="tight", pad_inches=0.14)
plt.close(fig)
print(f"\nwritten to {p}")
