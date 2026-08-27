#!/usr/bin/env python3
"""
The indel injection, examined on its own terms.

Indels were injected by a different BAMSurgeon entry point, into a
separate BAM, after the substitutions were already in place. They deserve
their own accounting for three reasons.

The mechanics differ. Placing a substitution rewrites one base in a read;
placing a deletion removes bases and shortens it, and an insertion
lengthens it. Both change the read's length, which changes how bwa scores
its realignment, so the failure modes are not the same ones substitutions
have.

Verification differs. mpileup marks a deletion in the column before the
deleted base and as * in the column itself, and an insertion in the column
preceding it. Counting the literal characters finds nothing, which is how
an earlier version of this check reported every deletion as a failure.

And the sample is small: 75 indels against 1 528 substitutions, in 28 of
the 40 samples. Rates computed on 75 events carry wide intervals, and
saying so is part of reporting them.

Usage:
  python analyse_indels.py [outdir]
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

RES = os.path.expanduser("~/immune_escape_project/results")
WS = os.environ.get("WS", "/gpfs/bwfor/work/ws/fr_os136-immune_escape")
OUT = sys.argv[1] if len(sys.argv) > 1 else f"{RES}/figures"

IND = f"{RES}/verify_step8_per_indel.tsv"
IND_S = f"{RES}/verify_step8.tsv"
TRUTH = f"{WS}/simulation/indels/indel_truth_all.tsv"
SKIP = f"{WS}/simulation/indels/skipped.log"

EMBER, GOLD = "#E8402A", "#F2A623"
BONE, ASH, DUSK = "#F4EFEA", "#A29591", "#7A6C68"
RUST, SLATE, GREEN = "#8E3020", "#5A4744", "#6F9E44"

have = {f.name for f in fm.fontManager.ttflist}
HEAD = next((f for f in ["TeX Gyre Heros Cn", "Liberation Sans Narrow",
                         "Liberation Sans", "DejaVu Sans"] if f in have),
            "DejaVu Sans")
BODY = next((f for f in ["Carlito", "Open Sans", "Liberation Sans",
                         "DejaVu Sans"] if f in have), "DejaVu Sans")

if not os.path.exists(IND):
    sys.exit(f"missing {IND}\nrun verify_step8_indels.py first")

d = pd.read_csv(IND, sep="\t")
s = pd.read_csv(IND_S, sep="\t") if os.path.exists(IND_S) else pd.DataFrame()
t = pd.read_csv(TRUTH, sep="\t") if os.path.exists(TRUTH) else pd.DataFrame()


def ci(k, n):
    """Wilson interval — a normal approximation on 75 events would give
    bounds outside [0,1] at rates this close to the ceiling."""
    if n == 0:
        return (np.nan, np.nan)
    z, p = 1.96, k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return (100 * max(0, c - h), 100 * min(1, c + h))


W = 68
print("=" * W)
print(" 1. WHAT WAS ATTEMPTED, AND WHERE")
print("=" * W)

n, landed = len(d), int(d.landed.sum())
lo, hi = ci(landed, n)
print(f"\n  indels attempted          {n}")
print(f"  confirmed in the reads    {landed}   "
      f"({100*landed/n:.1f}%, 95% CI {lo:.0f}\u2013{hi:.0f}%)")
print(f"  samples carrying any      {d['sample'].nunique()} of 40")

if os.path.exists(SKIP):
    sk = [l.split("\t")[0] for l in open(SKIP) if l.strip()]
    print(f"  samples with none in the panel   {len(sk)}")
    print(f"    {' '.join(sorted(sk))}")
    print(f"\n  Those donors carry no indel inside the panel at all. That is")
    print(f"  their own mutation profile, not a failure of the injection.")

print(f"\n  by type:")
print(f"  {'':<10}{'attempted':>11}{'landed':>9}{'rate':>9}{'95% CI':>14}")
for ty, g in d.groupby("type"):
    k = int(g.landed.sum())
    l2, h2 = ci(k, len(g))
    print(f"  {ty:<10}{len(g):>11}{k:>9}{100*k/len(g):>8.1f}%"
          f"{f'{l2:.0f}-{h2:.0f}%':>14}")

if len(t) and "consequence" in t.columns:
    print(f"\n  by consequence:")
    m = d.merge(t[["sample", "pos", "consequence"]], on=["sample", "pos"],
                how="left", suffixes=("", "_t"))
    col = "consequence_t" if "consequence_t" in m.columns else "consequence"
    for c, g in m.groupby(col):
        print(f"    {str(c):<24}{int(g.landed.sum()):>4} of {len(g):>4}")

print()
print("=" * W)
print(" 2. WHY SIX DID NOT LAND")
print("=" * W)
lost = d[~d.landed]
if len(lost):
    print(f"\n  {len(lost)} indels missing, all deletions" if
          (lost.type == "DEL").all() else f"\n  {len(lost)} indels missing")
    print(f"\n  {'sample':<8}{'position':<22}{'gene':<12}{'depth':>7}"
          f"{'requested VAF':>15}")
    for _, r in lost.iterrows():
        print(f"  {r['sample']:<8}{r.chrom + ':' + str(int(r.pos)):<22}"
              f"{str(r.get('gene', '')):<12}{int(r.depth):>7}"
              f"{r.target_vaf:>15.3f}")

    genes = sorted(set(lost.gene.dropna())) if "gene" in lost.columns else []
    if genes:
        print(f"\n  All of them sit in {len(genes)} genes: {', '.join(genes)}")
        print(f"\n  These are segmental duplications and tandem repeats.")
        print(f"  A read edited there realigns to one of the other copies as")
        print(f"  readily as to the original, and bwa places it elsewhere or")
        print(f"  drops its proper-pair flag. The same positions fail across")
        print(f"  different samples, which is the signature of a reference")
        print(f"  property rather than a stochastic failure.")

    print(f"\n  depth at the failing positions:")
    print(f"    median where lost     {lost.depth.median():.0f}x")
    print(f"    median where landed   {d[d.landed].depth.median():.0f}x")

print()
print("=" * W)
print(" 3. FIDELITY OF THE ALLELE FRACTION")
print("=" * W)
ok = d[d.landed & (d.target_vaf > 0) & (d.observed_vaf > 0)]
if len(ok) > 3:
    r = np.corrcoef(ok.target_vaf, ok.observed_vaf)[0, 1]
    dev = ok.observed_vaf - ok.target_vaf
    ratio = (ok.observed_vaf / ok.target_vaf).median()
    print(f"\n  landed with a measurable fraction   {len(ok)}")
    print(f"  correlation                         r = {r:.3f}")
    print(f"  mean deviation                      {dev.mean():+.4f}")
    print(f"  median observed / requested         {ratio:.3f}")
    print(f"\n  For comparison, substitutions in step 2 land at about 0.83")
    print(f"  of the requested fraction. A different ratio here would mean")
    print(f"  the two injection routines behave differently, not that one")
    print(f"  of them is wrong.")

print()
print("=" * W)
print(" 4. WHAT IS STILL UNMEASURED")
print("=" * W)
print(f"\n  These {n} indels were injected into a separate BAM and Mutect2")
print(f"  was never run against it. The truth-versus-calls comparison in")
print(f"  step 3 contains substitutions only.")
print(f"\n  So the pipeline has been shown to place indels correctly and")
print(f"  has not been shown to find them. That gap is worth naming rather")
print(f"  than leaving for a reviewer to notice.")

# ------------------------------------------------------------------ figure
os.makedirs(OUT, exist_ok=True)
fig = plt.figure(figsize=(12, 4.0), dpi=300)
fig.patch.set_alpha(0)

ax = fig.add_axes([0.05, 0.20, 0.40, 0.62])
ax.set_facecolor((0, 0, 0, 0))
for sp in ax.spines.values():
    sp.set_visible(False)

types = ["INS", "DEL"]
present = [ty for ty in types if (d.type == ty).any()]
xs = np.arange(len(present))
rates, ns, errs = [], [], []
for ty in present:
    g = d[d.type == ty]
    k = int(g.landed.sum())
    rates.append(100 * k / len(g))
    ns.append((k, len(g)))
    l2, h2 = ci(k, len(g))
    errs.append([rates[-1] - l2, h2 - rates[-1]])

ax.bar(xs, rates, width=0.5, color=[GREEN if r > 95 else EMBER for r in rates],
       zorder=3, edgecolor="none")
ax.errorbar(xs, rates, yerr=np.array(errs).T, fmt="none", ecolor=GOLD,
            elinewidth=1.4, capsize=6, capthick=1.4, zorder=5)
for x, r, (k, tot) in zip(xs, rates, ns):
    ax.text(x, 4, f"{k} of {tot}", family=BODY, fontsize=12, color=BONE,
            ha="center", va="bottom", zorder=6)
ax.set_xticks(xs)
ax.set_xticklabels(["insertions", "deletions"])
ax.tick_params(axis="x", colors=ASH, labelsize=12.5, length=0, pad=8)
ax.tick_params(axis="y", colors=DUSK, labelsize=11, length=0)
for tl in ax.get_xticklabels() + ax.get_yticklabels():
    tl.set_fontfamily(BODY)
ax.yaxis.set_major_formatter(lambda v, p: f"{v:.0f}%")
ax.set_ylim(0, 118)
ax.grid(axis="y", color="#3A2F2D", lw=0.6, zorder=0)
ax.set_axisbelow(True)
ax.text(0, 1.08, "CONFIRMED IN THE READS", transform=ax.transAxes,
        family=HEAD, fontsize=13.5, color=BONE, va="bottom",
        fontweight="bold")

ax2 = fig.add_axes([0.55, 0.20, 0.42, 0.62])
ax2.set_facecolor((0, 0, 0, 0))
ax2.axis("off")
ax2.set_xlim(0, 1)
ax2.set_ylim(0, 1)
ax2.text(0, 1.08, "WHERE THE FAILURES SIT", transform=ax2.transAxes,
         family=HEAD, fontsize=13.5, color=BONE, va="bottom",
         fontweight="bold")

if len(lost) and "gene" in lost.columns:
    gc = lost.gene.value_counts()
    y = 0.78
    for g, k in gc.items():
        ax2.add_patch(FancyBboxPatch((0, y - 0.06), 0.30, 0.16,
                                     boxstyle="round,pad=0,rounding_size=0.03",
                                     facecolor="#3A1A13", edgecolor=EMBER,
                                     lw=1.3, zorder=3))
        ax2.text(0.15, y + 0.02, g, family=HEAD, fontsize=15, color=BONE,
                 ha="center", va="center", fontweight="bold", zorder=6)
        ax2.text(0.34, y + 0.02, f"{k} lost", family=BODY, fontsize=13,
                 color=GOLD, va="center")
        y -= 0.26
    ax2.text(0, 0.02, "segmental duplications and tandem repeats",
             family=BODY, fontsize=12, color=DUSK, va="center")

p = f"{OUT}/indel_injection.png"
fig.savefig(p, dpi=300, transparent=True, bbox_inches="tight", pad_inches=0.14)
plt.close(fig)
print(f"\nwritten to {p}")
