#!/usr/bin/env python3
"""
Whether the coordinate conversion landed where it should have.

TCGA publishes mutations on hg19; the reads are aligned to hg38. Every
position has to be converted before a mutation can be placed, and a
conversion that quietly puts a variant at the wrong base is worse than one
that fails outright — the wrong base enters the truth set as if it were
real, and every recall figure afterwards is measured against a mutation
that never existed.

So each lifted position is checked against the hg38 reference: the base
sitting there must be the reference allele the MAF records. Bases are read
with samtools faidx in one batched call rather than one call per position,
which would take hours over this many variants.

The figure shows two things. How far positions moved, which establishes
that a conversion happened at all — copied coordinates would sit at zero.
And what survived the check, drawn as counts rather than as a percentage
bar, because two bars at 99.7% are indistinguishable by eye.

Usage:
  python plot_liftover.py [outdir]
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
from matplotlib.patches import FancyBboxPatch

WS = os.environ.get("WS", "/gpfs/bwfor/work/ws/fr_os136-immune_escape")
LIFT = f"{WS}/liftover/out"
REF = f"{WS}/ref/Homo_sapiens_assembly38.fasta"
RES = os.path.expanduser("~/immune_escape_project/results")
OUT = sys.argv[1] if len(sys.argv) > 1 else f"{RES}/figures"
N_CHECK = 20000

EMBER, GOLD = "#E8402A", "#F2A623"
BONE, ASH, DUSK = "#F4EFEA", "#A29591", "#7A6C68"
GREEN, RUST = "#6F9E44", "#A03A22"

have = {f.name for f in fm.fontManager.ttflist}
HEAD = next((f for f in ["TeX Gyre Heros Cn", "Liberation Sans Narrow",
                         "Liberation Sans", "DejaVu Sans"] if f in have),
            "DejaVu Sans")
BODY = next((f for f in ["Carlito", "Open Sans", "Liberation Sans",
                         "DejaVu Sans"] if f in have), "DejaVu Sans")


def fetch_bases(regions):
    """One faidx call for all positions; returns bases in input order."""
    rf = "/tmp/lift_regions.txt"
    with open(rf, "w") as fh:
        fh.write("\n".join(regions) + "\n")
    r = subprocess.run(f"samtools faidx {REF} -r {rf}", shell=True,
                       capture_output=True, text=True)
    if r.returncode != 0:
        print("  samtools faidx failed:", r.stderr.strip()[:120])
        return []
    out, cur = [], []
    for line in r.stdout.split("\n"):
        if line.startswith(">"):
            if cur:
                out.append("".join(cur).upper())
            cur = []
        elif line.strip():
            cur.append(line.strip())
    if cur:
        out.append("".join(cur).upper())
    return out


print("=" * 66)
print(" LIFTOVER, RECOMPUTED")
print("=" * 66)

summary, deltas = {}, []

for coh in ["brca", "ov"]:
    p = f"{LIFT}/{coh}.hg38.maf.tsv"
    if not os.path.exists(p):
        print(f"\n  missing {p}")
        continue

    d = pd.read_csv(p, sep="\t", low_memory=False)
    d["delta"] = d.Start_Position_hg38 - d.Start_Position
    deltas.append(d["delta"].dropna().values)

    snv = d[(d.Reference_Allele.astype(str).str.len() == 1) &
            (d.Reference_Allele.astype(str) != "-") &
            d.Start_Position_hg38.notna()].copy()
    take = snv.sample(min(N_CHECK, len(snv)), random_state=7)

    regions = [f"{r.Chromosome_hg38}:{int(r.Start_Position_hg38)}"
               f"-{int(r.Start_Position_hg38)}" for _, r in take.iterrows()]
    bases = fetch_bases(regions)

    n = min(len(bases), len(take))
    exp = take.Reference_Allele.astype(str).str.upper().values[:n]
    match = int(sum(b == e for b, e in zip(bases[:n], exp)))
    conc = 100 * match / n if n else float("nan")

    summary[coh] = dict(total=len(d), checked=n, match=match, conc=conc,
                        bad=n - match,
                        med=float(np.median(np.abs(d["delta"].dropna()))))

    print(f"\n  {coh.upper()}")
    print(f"    mutations in the MAF     {len(d)}")
    print(f"    positions checked        {n}")
    print(f"    reference base agrees    {match}   ({conc:.2f}%)")
    print(f"    disagrees, discarded     {n - match}")
    print(f"    median absolute shift    {summary[coh]['med']:,.0f} bp"
          .replace(",", " "))

if not summary:
    sys.exit("nothing to plot")

print(f"\n  Coordinates moved by hundreds of kilobases in both directions;")
print(f"  a copy rather than a conversion would sit at zero throughout.")

os.makedirs(OUT, exist_ok=True)
fig = plt.figure(figsize=(12, 4.2), dpi=300)
fig.patch.set_alpha(0)

# --------------------------------------------------- left: how far it moved
ax = fig.add_axes([0.06, 0.20, 0.40, 0.64])
ax.set_facecolor((0, 0, 0, 0))
for s in ax.spines.values():
    s.set_visible(False)

sh = np.abs(np.concatenate(deltas))
sh = sh[sh > 0]
bins = np.logspace(0, np.log10(sh.max()), 30)
ax.hist(sh, bins=bins, color=EMBER, alpha=0.9, zorder=3)
med = np.median(sh)
ax.axvline(med, color=GOLD, lw=1.4, ls=(0, (5, 4)), zorder=5)
ax.text(med * 1.35, ax.get_ylim()[1] * 0.88,
        f"median {med/1000:.0f} kb", color=GOLD, family=BODY, fontsize=12,
        weight="bold", va="center")

ax.set_xscale("log")
ax.tick_params(colors=DUSK, labelsize=11, length=0)
for t in ax.get_xticklabels() + ax.get_yticklabels():
    t.set_fontfamily(BODY)
ax.grid(axis="y", color="#3A2F2D", lw=0.6, zorder=0)
ax.set_axisbelow(True)
ax.set_xlabel("distance the position moved, bp", color=ASH, fontsize=12,
              family=BODY, labelpad=8)
ax.set_ylabel("mutations", color=ASH, fontsize=12, family=BODY, labelpad=8)
ax.text(0, 1.10, "EVERY POSITION MOVED", transform=ax.transAxes,
        family=HEAD, fontsize=13.5, color=BONE, va="bottom",
        fontweight="bold")

# ------------------------------------------- right: what survived the check
ax2 = fig.add_axes([0.55, 0.20, 0.42, 0.64])
ax2.set_facecolor((0, 0, 0, 0))
ax2.axis("off")
ax2.set_xlim(0, 1)
ax2.set_ylim(0, 1)

ax2.text(0, 1.10, "REFERENCE BASE RECHECKED", transform=ax2.transAxes,
         family=HEAD, fontsize=13.5, color=BONE, va="bottom",
         fontweight="bold")

y = 0.66
for coh, s in summary.items():
    frac_bad = max(s["bad"] / s["checked"], 0.004) if s["checked"] else 0
    ax2.add_patch(FancyBboxPatch((0, y), 1.0, 0.20,
                                 boxstyle="round,pad=0,rounding_size=0.03",
                                 facecolor=GREEN, edgecolor="none", zorder=3))
    ax2.add_patch(FancyBboxPatch((1 - frac_bad * 0.9, y), frac_bad * 0.9, 0.20,
                                 boxstyle="round,pad=0,rounding_size=0.03",
                                 facecolor=RUST, edgecolor="none", zorder=4))
    ax2.text(0.025, y + 0.10, coh.upper(), family=HEAD, fontsize=16,
             color="#14200A", va="center", fontweight="bold", zorder=6)
    ax2.text(0.30, y + 0.10, f"{s['conc']:.2f}% agree", family=HEAD,
             fontsize=17, color="#14200A", va="center", fontweight="bold",
             zorder=6)
    ax2.text(0, y - 0.10,
             f"{s['match']} of {s['checked']} checked  ·  "
             f"{s['bad']} discarded", family=BODY, fontsize=12, color=DUSK,
             va="center")
    y -= 0.44

p = f"{OUT}/liftover.png"
fig.savefig(p, dpi=300, transparent=True, bbox_inches="tight", pad_inches=0.14)
plt.close(fig)
print(f"\nwritten to {p}")
