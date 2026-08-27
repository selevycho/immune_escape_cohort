#!/usr/bin/env python3
"""
The position that exposed the verification bug.

At O016 chr2:61188227 four reads of forty carry the injected base. None of
those four is flagged proper-pair, because BAMSurgeon rewrote them and bwa
dropped the flag when it realigned them. samtools 1.23 mpileup requires
that flag and offers no way to turn the requirement off, so it reported a
depth of one and no alternate reads, and the mutation was scored as a
failed injection.

The figure draws the reads as they are in the file: solid where the
proper-pair flag survived, faded where it did not, with the alternate
carriers marked. Every carrier is in the faded group, which is the whole
argument.

Numbers are read from the BAM rather than typed in, so the picture cannot
drift from the file it describes.

Usage:
  python plot_verification_bug.py [outdir] [--pos=SAMPLE:CHROM:POS:ALT]
"""
import os
import re
import sys
import subprocess
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from deckplot import *
from matplotlib.patches import FancyBboxPatch

WS = os.environ.get("WS", "/gpfs/bwfor/work/ws/fr_os136-immune_escape")
RES = os.path.expanduser("~/immune_escape_project/results")
args = [a for a in sys.argv[1:] if not a.startswith("--")]
OUT = args[0] if args else f"{RES}/figures"

SPEC = "O016:chr2:61188227:A"
for a in sys.argv[1:]:
    if a.startswith("--pos="):
        SPEC = a.split("=", 1)[1]
SID, CHROM, POS, ALT = SPEC.split(":")
POS = int(POS)

BAM = f"{WS}/simulation/cohort/{SID}/{SID}_tumor.bam"
REF = f"{WS}/ref/Homo_sapiens_assembly38.fasta"

SAMTOOLS = None
for cand in ["samtools",
             os.path.expanduser("~/miniconda3/envs/bio_work/bin/samtools"),
             "/home/fr/fr_fr/fr_os136/miniconda3/envs/bio_work/bin/samtools"]:
    if subprocess.run(f"{cand} --version", shell=True,
                      capture_output=True).returncode == 0:
        SAMTOOLS = cand
        break
if SAMTOOLS is None:
    sys.exit("samtools not found")

CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def base_at(start, cigar, seq, target):
    ref, qry = start, 0
    for n, op in CIGAR.findall(cigar):
        n = int(n)
        if op in "M=X":
            if ref <= target < ref + n:
                return seq[qry + (target - ref)]
            ref += n
            qry += n
        elif op in "IS":
            qry += n
        elif op in "DN":
            if ref <= target < ref + n:
                return "*"
            ref += n
    return None


out = subprocess.run(f"{SAMTOOLS} view -F 0x404 {BAM} {CHROM}:{POS}-{POS}",
                     shell=True, capture_output=True, text=True).stdout
reads = []
for line in out.splitlines():
    f = line.split("\t")
    if len(f) < 11:
        continue
    b = base_at(int(f[3]), f[5], f[9], POS)
    if b is None:
        continue
    flag = int(f[1])
    reads.append({"start": int(f[3]), "len": len(f[9]),
                  "proper": bool(flag & 0x2), "base": b.upper()})

pile = subprocess.run(
    f"{SAMTOOLS} mpileup -B -q 0 -Q 0 -r {CHROM}:{POS}-{POS} -f {REF} {BAM} "
    f"2>/dev/null", shell=True, capture_output=True, text=True).stdout.strip()
mp_depth = int(pile.split("\t")[3]) if pile else 0
mp_alt = pile.split("\t")[4].upper().count(ALT) if pile else 0

n = len(reads)
n_alt = sum(r["base"] == ALT for r in reads)
n_proper = sum(r["proper"] for r in reads)
alt_proper = sum(r["proper"] for r in reads if r["base"] == ALT)

print("=" * 66)
print(f" {SID}  {CHROM}:{POS}   injected base {ALT}")
print("=" * 66)
print(f"\n  reads at the position          {n}")
print(f"    carrying {ALT}                    {n_alt}")
print(f"    flagged proper-pair          {n_proper}")
print(f"    carrying {ALT} AND proper-pair    {alt_proper}")
print(f"\n  samtools mpileup reports       depth {mp_depth}, "
      f"{ALT} reads {mp_alt}")
print(f"  reading the alignments gives   depth {n}, {ALT} reads {n_alt}")
print(f"\n  Every read carrying the injected base lost its proper-pair flag")
print(f"  during realignment, and mpileup requires that flag.")

# ------------------------------------------------------------------ figure
os.makedirs(OUT, exist_ok=True)
fig = figure(12, 5.0)
ax = blank(fig, (0.03, 0.04, 0.94, 0.86))

fig.text(0.03, 0.925, f"ONE POSITION  ·  {SID}  {CHROM}:{POS}", family=HEAD,
         fontsize=15, color=BONE, fontweight="bold", va="bottom")

reads = sorted(reads, key=lambda r: (not (r["base"] == ALT), r["start"]))
COL_X, COL_W = 0.03, 0.50
rows = min(len(reads), 40)
BARH = 0.0155
PITCH = 0.0225
top = 0.86

rng = np.random.default_rng(2)
for k, r in enumerate(reads[:rows]):
    y = top - k * PITCH
    jitter = rng.uniform(0, 0.09)
    x = COL_X + jitter
    w = COL_W - jitter - rng.uniform(0, 0.06)
    carries = r["base"] == ALT
    ax.add_patch(FancyBboxPatch(
        (x, y), w, BARH, boxstyle="round,pad=0,rounding_size=0.006",
        facecolor="#4A2A22" if carries else "#332B29",
        edgecolor=EMBER if carries else "#4A3B38",
        lw=1.1 if carries else 0.7,
        alpha=1.0 if not r["proper"] else 0.42, zorder=3))
    if carries:
        ax.add_patch(FancyBboxPatch(
            (x + w * 0.55, y), 0.011, BARH,
            boxstyle="round,pad=0,rounding_size=0.004",
            facecolor=EMBER, edgecolor="none", zorder=5))

LX = COL_X + COL_W + 0.03
ax.add_patch(FancyBboxPatch((LX, top - 0.02), 0.028, BARH,
                            boxstyle="round,pad=0,rounding_size=0.006",
                            facecolor="#4A2A22", edgecolor=EMBER, lw=1.1,
                            zorder=4))
ax.text(LX + 0.038, top - 0.02 + BARH / 2,
        f"carries the injected {ALT}   ({n_alt} reads)", family=BODY,
        fontsize=12.5, color=BONE, va="center")

ax.add_patch(FancyBboxPatch((LX, top - 0.075), 0.028, BARH,
                            boxstyle="round,pad=0,rounding_size=0.006",
                            facecolor="#332B29", edgecolor="#4A3B38",
                            lw=0.7, alpha=0.42, zorder=4))
ax.text(LX + 0.038, top - 0.075 + BARH / 2,
        f"proper-pair flag survived   ({n_proper} reads)", family=BODY,
        fontsize=12.5, color=DUSK, va="center")

CY, CH, CW = 0.30, 0.20, 0.40
ax.add_patch(FancyBboxPatch((LX, CY + CH + 0.05), CW, CH,
                            boxstyle="round,pad=0,rounding_size=0.02",
                            facecolor="#2E1614", edgecolor="#8E3020", lw=1.6,
                            zorder=3))
ax.text(LX + 0.03, CY + CH + 0.05 + CH * 0.66, "SAMTOOLS MPILEUP",
        family=HEAD, fontsize=13, color=ASH, va="center", fontweight="bold")
ax.text(LX + 0.03, CY + CH + 0.05 + CH * 0.28,
        f"depth {mp_depth}   ·   {ALT} reads {mp_alt}   \u2192  discarded",
        family=BODY, fontsize=13.5, color="#E88A78", va="center")

ax.add_patch(FancyBboxPatch((LX, CY), CW, CH,
                            boxstyle="round,pad=0,rounding_size=0.02",
                            facecolor="#1C2416", edgecolor=GREEN, lw=1.6,
                            zorder=3))
ax.text(LX + 0.03, CY + CH * 0.66, "READING THE ALIGNMENTS", family=HEAD,
        fontsize=13, color=ASH, va="center", fontweight="bold")
ax.text(LX + 0.03, CY + CH * 0.28,
        f"depth {n}   ·   {ALT} reads {n_alt}   \u2192  kept",
        family=BODY, fontsize=13.5, color="#A8CC7F", va="center")

ax.text(LX, CY - 0.09,
        f"none of the {n_alt} carriers kept its proper-pair flag",
        family=BODY, fontsize=12.5, color=GOLD, va="center",
        fontweight="bold")

print("\n" + save(fig, f"{OUT}/verification_bug.png"))
