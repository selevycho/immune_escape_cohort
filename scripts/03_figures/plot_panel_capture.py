#!/usr/bin/env python3
"""
What share of a donor's mutations the panel actually captures.

The slide claims a median of 12.2% with a range of 6.2 to 23.9. Those
numbers came out of check_panel.py in an earlier session; this recomputes
them from panel_audit.tsv and draws the figure from the same table, so the
picture and the caption cannot drift apart.

Bars are sorted by capture rather than by sample name. Sorted, the shape of
the distribution is the point: most donors sit in a narrow band and the
extremes are visibly extremes, which a name-ordered chart hides.

Output is a transparent PNG in the deck's palette.

Usage:
  python plot_panel_capture.py [outdir]
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.patches import FancyBboxPatch

AUDIT = os.path.expanduser("~/immune_escape_project/results/panel_audit.tsv")
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser(
    "~/immune_escape_project/results/figures")

EMBER, GOLD = "#E8402A", "#F2A623"
BONE, ASH, DUSK = "#F4EFEA", "#A29591", "#7A6C68"
BRCA_C, OV_C = "#C0341F", "#7E4335"

# the deck uses Bebas Neue over Open Sans; neither is guaranteed on a
# cluster, so fall back through the nearest condensed and humanist faces
have = {f.name for f in fm.fontManager.ttflist}
HEAD = next((f for f in ["TeX Gyre Heros Cn", "Liberation Sans Narrow",
                         "DejaVu Sans Condensed", "Liberation Sans",
                         "DejaVu Sans"] if f in have), "DejaVu Sans")
BODY = next((f for f in ["Carlito", "Open Sans", "Liberation Sans",
                         "DejaVu Sans"] if f in have), "DejaVu Sans")

if not os.path.exists(AUDIT):
    sys.exit(f"missing {AUDIT}\nrun check_panel.py first")

d = pd.read_csv(AUDIT, sep="\t").sort_values("pct_captured").reset_index(drop=True)

print("=" * 68)
print(" PANEL CAPTURE, RECOMPUTED FROM panel_audit.tsv")
print("=" * 68)
print(f"\n  samples                {len(d)}")
print(f"  median capture         {d.pct_captured.median():.1f}%")
print(f"  range                  {d.pct_captured.min():.1f}% "
      f"({d.loc[d.pct_captured.idxmin(), 'sample']})  to  "
      f"{d.pct_captured.max():.1f}% "
      f"({d.loc[d.pct_captured.idxmax(), 'sample']})")
print(f"  mutations in panel     {int(d.in_panel.sum())}")
print(f"  donor mutations usable {int(d.vaf_filtered.sum())}")
print(f"  injected total         {int(d.injected.sum())}")
print(f"  genes hit per patient  median {d.genes_hit.median():.0f}, "
      f"range {d.genes_hit.min()}-{d.genes_hit.max()}")
for coh, g in d.groupby("cohort"):
    print(f"    {coh.upper():<6} median {g.pct_captured.median():.1f}%, n={len(g)}")

lo, hi = d.iloc[0], d.iloc[-1]
med = d.pct_captured.median()
top = d.pct_captured.max()

os.makedirs(OUT, exist_ok=True)
fig = plt.figure(figsize=(12, 4.6), dpi=300)
ax = fig.add_axes([0.055, 0.16, 0.92, 0.78])
for s in ax.spines.values():
    s.set_visible(False)
ax.set_facecolor((0, 0, 0, 0))
fig.patch.set_alpha(0)

x = np.arange(len(d))
ax.bar(x, d.pct_captured, width=0.72, zorder=3, edgecolor="none",
       color=[BRCA_C if c == "brca" else OV_C for c in d.cohort])

ax.axhline(med, color=GOLD, lw=1.3, ls=(0, (5, 4)), zorder=4)
ax.text(-2.0, med + 0.5, f"median {med:.1f}%", color=GOLD, fontsize=12,
        family=BODY, ha="left", va="bottom", weight="bold")

for i, r in [(0, lo), (len(d) - 1, hi)]:
    ax.text(i, r.pct_captured + 0.6, f"{r.pct_captured:.1f}%", color=BONE,
            fontsize=11, family=BODY, ha="center", va="bottom")
    ax.text(i, -top * 0.10, r["sample"], color=DUSK, fontsize=10.5,
            family=BODY, ha="center", va="top")

ax.set_xticks([])
ax.set_xlim(-2.4, len(d) + 1.4)
ax.set_ylim(0, top * 1.28)
ax.set_ylabel("mutations inside the panel", color=ASH, fontsize=12,
              family=BODY, labelpad=10)
ax.tick_params(axis="y", colors=DUSK, labelsize=11, length=0)
for t in ax.get_yticklabels():
    t.set_fontfamily(BODY)
ax.yaxis.set_major_formatter(lambda v, p: f"{v:.0f}%")
ax.grid(axis="y", color="#3A2F2D", lw=0.6, zorder=0)
ax.set_axisbelow(True)

lx = len(d) * 0.02
for name, c in [("BRCA", BRCA_C), ("OV", OV_C)]:
    ax.add_patch(FancyBboxPatch((lx, top * 1.14), 0.8, top * 0.045,
                                boxstyle="round,pad=0,rounding_size=0.2",
                                facecolor=c, edgecolor="none", zorder=5))
    ax.text(lx + 1.3, top * 1.16, name, color=ASH, fontsize=11.5,
            family=BODY, va="center")
    lx += 4.2

p = f"{OUT}/panel_capture.png"
fig.savefig(p, dpi=300, transparent=True, bbox_inches="tight", pad_inches=0.12)
plt.close(fig)
print(f"\nwritten to {p}")
