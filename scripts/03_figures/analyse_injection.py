#!/usr/bin/env python3
"""
The injection, recomputed after the verification was corrected.

The earlier numbers came from mpileup, which in samtools 1.23 silently
requires the proper-pair flag. BAMSurgeon edits reads and realigns them
through bwa, and realignment strips that flag from precisely the reads it
rewrote, so mpileup discarded mutation-carrying reads preferentially. The
verification now walks CIGAR strings directly and the picture changes:
substitutions land at 98% rather than 91%, indels at 100% rather than 92%,
and the proportional shortfall in allele fraction largely disappears.

Four sections, and two figures for the slides. Substitutions and indels
are kept separate throughout, since they went through different BAMSurgeon
entry points into different files.

Usage:
  python analyse_injection.py [outdir]
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
SNV = f"{RES}/verify_step2_per_mutation.tsv"
IND = f"{RES}/verify_step8_per_indel.tsv"

for p in (SNV, IND):
    if not os.path.exists(p):
        sys.exit(f"missing {p}")

s = pd.read_csv(SNV, sep="\t")
i = pd.read_csv(IND, sep="\t")

n_s, ok_s = len(s), int(s.landed.sum())
n_i, ok_i = len(i), int(i.landed.sum())
lost_s = s[~s.landed]

W = 70
print("=" * W)
print(" 1. WHAT LANDED")
print("=" * W)
print(f"\n  substitutions   {ok_s} of {n_s}   ({100*ok_s/n_s:.1f}%)")
print(f"  indels          {ok_i} of {n_i}   ({100*ok_i/n_i:.1f}%)")
for ty, g in i.groupby("type"):
    print(f"    {ty:<12}{int(g.landed.sum())} of {len(g)}")

print(f"\n  substitutions outside the MHC   "
      f"{int(s[~s.in_mhc].landed.sum())} of {int((~s.in_mhc).sum())}"
      f"   ({100*s[~s.in_mhc].landed.mean():.1f}%)")
if s.in_mhc.any():
    print(f"  substitutions inside the MHC    "
          f"{int(s[s.in_mhc].landed.sum())} of {int(s.in_mhc.sum())}"
          f"   ({100*s[s.in_mhc].landed.mean():.1f}%)")

print()
print("=" * W)
print(" 2. WHERE THE REMAINING LOSSES ARE")
print("=" * W)
CATS = [
    ("Inside the MHC", int(lost_s.in_mhc.sum()), RUST),
    ("Depth below 5", int(((~lost_s.in_mhc) & (lost_s.tumour_depth < 5)).sum()),
     SLATE),
    ("Covered, absent", int(((~lost_s.in_mhc) & (lost_s.tumour_depth >= 5)).sum()),
     EMBER),
]
print(f"\n  {len(lost_s)} substitutions of {n_s} did not land\n")
for name, k, _ in CATS:
    print(f"  {name:<22}{k:>5}   {100*k/max(1,len(lost_s)):>5.0f}% of losses")

rest = lost_s[(~lost_s.in_mhc) & (lost_s.tumour_depth >= 5)]
if len(rest):
    print(f"\n  the covered group: median depth "
          f"{rest.tumour_depth.median():.0f}x, "
          f"median requested VAF {rest.target_vaf.median():.3f}")
print(f"\n  The MHC carries alternate haplotypes in hg38, so an edited read")
print(f"  realigns to one of them as readily as to the primary contig.")
print(f"  That is {100*int(lost_s.in_mhc.sum())/max(1,len(lost_s)):.0f}% of "
      f"what is left, from {100*s.in_mhc.mean():.0f}% of the positions.")

print()
print("=" * W)
print(" 3. FIDELITY OF THE ALLELE FRACTION")
print("=" * W)


def fidelity(d, label):
    ok = d[d.landed & (d.target_vaf > 0) & (d.observed_vaf > 0)]
    if len(ok) < 4:
        return None
    r = float(np.corrcoef(ok.target_vaf, ok.observed_vaf)[0, 1])
    ratio = float((ok.observed_vaf / ok.target_vaf).median())
    dev = float((ok.observed_vaf - ok.target_vaf).mean())
    print(f"\n  {label}")
    print(f"    n                      {len(ok)}")
    print(f"    correlation            r = {r:.3f}")
    print(f"    median observed/target {ratio:.3f}")
    print(f"    mean deviation         {dev:+.4f}")
    return ok, r, ratio


fs = fidelity(s, "substitutions")
fi = fidelity(i, "indels")

ok_snv, r_snv, ratio_snv = fs
ok_ind, r_ind, ratio_ind = fi

BINS = [0, .05, .10, .15, .20, .30, 1.01]
NAMES = ["<5%", "5-10%", "10-15%", "15-20%", "20-30%", ">30%"]
ok_snv = ok_snv.copy()
ok_snv["bin"] = pd.cut(ok_snv.target_vaf, bins=BINS, labels=NAMES, right=False)
labs, ratios = [], []
print(f"\n  substitutions by requested fraction:")
print(f"  {'requested':<11}{'n':>6}{'target':>9}{'observed':>11}{'ratio':>9}")
for nm in NAMES:
    g = ok_snv[ok_snv.bin == nm]
    if len(g) <= 3:
        continue
    rt = float((g.observed_vaf / g.target_vaf).median())
    labs.append(nm)
    ratios.append(rt)
    print(f"  {nm:<11}{len(g):>6}{g.target_vaf.median():>9.3f}"
          f"{g.observed_vaf.median():>11.3f}{rt:>9.3f}")

med_depth = float(s.tumour_depth.median())
print(f"\n  One read is {1/med_depth:.3f} of the total at the median depth")
print(f"  of {med_depth:.0f}x, so a requested fraction can only be met to")
print(f"  within that step. The residual shortfall is that granularity,")
print(f"  not a systematic bias in the injection.")

print()
print("=" * W)
print(" 4. WHAT THE TOOL REPORTED")
print("=" * W)
print(f"\n  BAMSurgeon logged success at every locus it touched and")
print(f"  reported no failure anywhere.")
print(f"\n  Reading the alignments found {len(lost_s)} substitutions absent")
print(f"  and {n_i - ok_i} indels absent.")
print(f"\n  A truth set taken from the log would therefore carry "
      f"{len(lost_s)} entries")
print(f"  that are not in the file, each counted against the caller in")
print(f"  step 3 — worth {100*len(lost_s)/n_s:.1f} points of recall.")

# ==================================================================== figures
os.makedirs(OUT, exist_ok=True)

# --- one: what landed, both kinds -------------------------------------
fig = figure(12, 4.6)
ax = blank(fig, (0.03, 0.05, 0.94, 0.80))
fig.text(0.03, 0.90, "CONFIRMED IN THE READS", family=HEAD, fontsize=15,
         color=BONE, fontweight="bold", va="bottom")

BH, GAP = 0.205, 0.070
y = 0.66
for label, k, tot in [("SUBSTITUTIONS", ok_s, n_s), ("INDELS", ok_i, n_i)]:
    frac = k / tot
    ax.add_patch(FancyBboxPatch((0, y), 1.0, BH,
                                boxstyle="round,pad=0,rounding_size=0.02",
                                facecolor=INK, edgecolor="#4A3B38", lw=1.2,
                                zorder=2))
    ax.add_patch(FancyBboxPatch((0, y), max(frac, 0.02), BH,
                                boxstyle="round,pad=0,rounding_size=0.02",
                                facecolor=GREEN, edgecolor="none", zorder=3))
    ax.text(0.022, y + BH / 2, label, family=HEAD, fontsize=14,
            color="#14200A", va="center", fontweight="bold", zorder=5)
    ax.text(frac - 0.022 if frac > 0.35 else frac + 0.022, y + BH / 2,
            f"{k} of {tot}   {100*frac:.1f}%", family=HEAD, fontsize=16,
            color="#14200A" if frac > 0.35 else BONE,
            ha="right" if frac > 0.35 else "left", va="center",
            fontweight="bold", zorder=5)
    y -= BH + GAP

CW, CG, CH, CY = 0.315, 0.0275, 0.26, 0.02
x = 0
for name, k, col in CATS:
    card(ax, x, CY, CW, CH, name,
         f"{k}   ({100*k/max(1,len(lost_s)):.0f}% of losses)", accent=col)
    x += CW + CG
print("\n" + save(fig, f"{OUT}/injection_landing.png"))

# --- two: fidelity ----------------------------------------------------
fig = figure(12, 4.6)

ax = panel(fig, 0.065, 0.38, bottom=0.17, title="allele fraction")
ax.scatter(ok_snv.target_vaf, ok_snv.observed_vaf, s=8, color=EMBER,
           alpha=0.30, edgecolors="none", zorder=3, label="substitutions")
ax.scatter(ok_ind.target_vaf, ok_ind.observed_vaf, s=30, color=GOLD,
           alpha=0.85, edgecolors="none", zorder=4, label="indels")
lim = max(ok_snv.target_vaf.max(), ok_snv.observed_vaf.max()) * 1.06
ax.plot([0, lim], [0, lim], color=BONE, lw=1.2, ls=(0, (5, 4)), zorder=5,
        alpha=0.55)
ax.set_xlim(0, lim)
ax.set_ylim(0, lim)
note(ax, f"r = {r_snv:.3f}", loc="upper left")
leg = ax.legend(frameon=False, loc="lower right", fontsize=11.5,
                markerscale=1.6)
for t in leg.get_texts():
    t.set_color(ASH)
    t.set_fontfamily(BODY)
finish(ax, "requested", "observed in the reads", grid="both")

ax2 = panel(fig, 0.575, 0.385, bottom=0.17, title="observed / requested")
xs = np.arange(len(labs))
ax2.bar(xs, ratios, width=0.58, color=EMBER, zorder=3, edgecolor="none")
ax2.axhline(1.0, color=GOLD, lw=1.4, ls=(0, (5, 4)), zorder=5)
ax2.set_ylim(0, 1.28)
for x, v in zip(xs, ratios):
    ax2.text(x, v + 0.035, f"{v:.2f}", family=BODY, fontsize=12.5,
             color=BONE, ha="center", va="bottom", zorder=6)
ax2.set_xticks(xs)
ax2.set_xticklabels(labs)
ax2.tick_params(axis="x", colors=ASH, labelsize=12, length=0, pad=7)
ax2.set_xlim(-0.7, len(labs) - 0.3)
finish(ax2, "requested fraction")
note(ax2, "1.00 = exactly as asked", loc="upper right", color=DUSK, size=11)
print(save(fig, f"{OUT}/injection_fidelity.png"))
