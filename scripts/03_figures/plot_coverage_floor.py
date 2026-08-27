#!/usr/bin/env python3
"""
The coverage that came back, and the detection floor it sets.

Two panels, because the second only means something once the first is on
screen. On the left, where the forty samples landed. On the right, what
that depth implies: a mutation present in a given fraction of the cells is
carried by depth times fraction reads, and below about three reads no
caller can tell it from a sequencing error.

The right panel is the reason the sensitivity curve in step 3 looks the
way it does, and drawing it here rather than there makes the point before
anyone can suspect it was fitted afterwards.

Usage:
  python plot_coverage_floor.py [outdir]
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm

RES = os.path.expanduser("~/immune_escape_project/results")
V1 = f"{RES}/verify_step1.tsv"
OUT = sys.argv[1] if len(sys.argv) > 1 else f"{RES}/figures"

EMBER, GOLD = "#E8402A", "#F2A623"
BONE, ASH, DUSK = "#F4EFEA", "#A29591", "#7A6C68"
BRCA_C, OV_C = "#C0341F", "#7E4335"
DEAD = "#6E2A20"

have = {f.name for f in fm.fontManager.ttflist}
HEAD = next((f for f in ["TeX Gyre Heros Cn", "Liberation Sans Narrow",
                         "Liberation Sans", "DejaVu Sans"] if f in have),
            "DejaVu Sans")
BODY = next((f for f in ["Carlito", "Open Sans", "Liberation Sans",
                         "DejaVu Sans"] if f in have), "DejaVu Sans")

if not os.path.exists(V1):
    sys.exit(f"missing {V1}\nrun verify_step1_slicing.py first")

d = pd.read_csv(V1, sep="\t").sort_values("panel_mean_depth").reset_index(drop=True)
med = d.panel_mean_depth.median()

print("=" * 62)
print(" COVERAGE AND WHAT IT ALLOWS")
print("=" * 62)
print(f"\n  samples              {len(d)}")
print(f"  panel depth          median {med:.1f}x, "
      f"range {d.panel_mean_depth.min():.1f}-{d.panel_mean_depth.max():.1f}x")
for c in ["mapped_pct", "dup_pct", "mean_mapq", "below_10x_pct"]:
    if c in d.columns and d[c].notna().any():
        v = d[c].dropna()
        print(f"  {c:<20} median {v.median():.2f}, "
              f"range {v.min():.2f}-{v.max():.2f}")
for c in ["hla_A_depth", "hla_B_depth", "hla_C_depth"]:
    if c in d.columns and d[c].notna().any():
        v = d[c].dropna()
        print(f"  {c:<20} median {v.median():.1f}x, "
              f"range {v.min():.1f}-{v.max():.1f}x")

VAFS = [0.05, 0.10, 0.15, 0.20, 0.30]
print(f"\n  at {med:.1f}x, a mutation is carried by:")
for v in VAFS:
    print(f"    VAF {v:.0%}   {med * v:.1f} reads")

os.makedirs(OUT, exist_ok=True)
fig = plt.figure(figsize=(12, 4.4), dpi=300)
fig.patch.set_alpha(0)

# ------------------------------------------------------------- left
ax = fig.add_axes([0.055, 0.17, 0.44, 0.74])
ax.set_facecolor((0, 0, 0, 0))
for s in ax.spines.values():
    s.set_visible(False)

x = np.arange(len(d))
ax.bar(x, d.panel_mean_depth, width=0.74, zorder=3, edgecolor="none",
       color=[BRCA_C if c == "brca" else OV_C for c in d.cohort])
ax.axhline(med, color=GOLD, lw=1.3, ls=(0, (5, 4)), zorder=5)
ax.text(len(d) * 0.02, med + 1.4, f"median {med:.1f}\u00d7", color=GOLD,
        family=BODY, fontsize=12.5, weight="bold", va="bottom", zorder=6)

ax.set_xticks([])
ax.set_xlim(-1.4, len(d) + 0.4)
ax.set_ylim(0, d.panel_mean_depth.max() * 1.22)
ax.tick_params(axis="y", colors=DUSK, labelsize=11, length=0)
for t in ax.get_yticklabels():
    t.set_fontfamily(BODY)
ax.yaxis.set_major_formatter(lambda v, p: f"{v:.0f}\u00d7")
ax.grid(axis="y", color="#3A2F2D", lw=0.6, zorder=0)
ax.set_axisbelow(True)
ax.text(0, d.panel_mean_depth.max() * 1.30, "DEPTH ACROSS THE PANEL",
        family=HEAD, fontsize=13.5, color=BONE, weight="bold", va="center")

# ------------------------------------------------------------ right
ax2 = fig.add_axes([0.575, 0.17, 0.40, 0.74])
ax2.set_facecolor((0, 0, 0, 0))
for s in ax2.spines.values():
    s.set_visible(False)

reads = [med * v for v in VAFS]
cols = [DEAD if r < 3 else (GOLD if r < 5 else EMBER) for r in reads]
xs = np.arange(len(VAFS))
ax2.bar(xs, reads, width=0.62, color=cols, zorder=3, edgecolor="none")
for xi, r in zip(xs, reads):
    ax2.text(xi, r + 0.35, f"{r:.1f}", family=BODY, fontsize=12,
             color=BONE, ha="center", va="bottom", zorder=6)

ax2.axhline(3, color=GOLD, lw=1.3, ls=(0, (5, 4)), zorder=5)
ax2.text(len(VAFS) - 0.45, 3.3, "three reads", color=GOLD, family=BODY,
         fontsize=12, weight="bold", ha="right", va="bottom", zorder=6)

ax2.set_xticks(xs)
ax2.set_xticklabels([f"{v:.0%}" for v in VAFS])
ax2.tick_params(axis="x", colors=ASH, labelsize=12, length=0, pad=8)
ax2.tick_params(axis="y", colors=DUSK, labelsize=11, length=0)
for t in ax2.get_xticklabels() + ax2.get_yticklabels():
    t.set_fontfamily(BODY)
ax2.set_ylim(0, max(reads) * 1.24)
ax2.set_xlim(-0.7, len(VAFS) - 0.3)
ax2.grid(axis="y", color="#3A2F2D", lw=0.6, zorder=0)
ax2.set_axisbelow(True)
ax2.text(-0.7, max(reads) * 1.32, "READS CARRYING THE MUTATION",
         family=HEAD, fontsize=13.5, color=BONE, weight="bold", va="center")
ax2.text(len(VAFS) / 2 - 0.5, -max(reads) * 0.15,
         "variant allele fraction", family=BODY, fontsize=12, color=ASH,
         ha="center", va="top")

p = f"{OUT}/coverage_floor.png"
fig.savefig(p, dpi=300, transparent=True, bbox_inches="tight", pad_inches=0.14)
plt.close(fig)
print(f"\nwritten to {p}")
