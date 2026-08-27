#!/usr/bin/env python3
"""
Is the Mutect2 block ready to be written up?

Step 3 was verified before the injection check was corrected, so its
recall was computed against a truth set of 1 528 substitutions of which
136 were thought to be absent from the files. That number is now 31, and
recall computed against the wrong denominator is wrong by the difference.

This does not rerun anything. It asks what exists, what it was measured
against, and what a reviewer could still ask that has no answer yet.

Five questions.

  Are the per-sample comparisons present for all forty?
  Which truth set were they scored against, and what does recall become
  against the corrected one?
  Is precision actually measured, or only assumed? A caller with no false
  positives is a claim that needs a denominator of its own.
  Does the matched normal behave as a control - does it stay clean where
  the tumour carries a mutation?
  What is not covered at all?

Usage:
  python audit_step3_readiness.py [workspace]
"""
import os
import sys
import glob
import numpy as np
import pandas as pd

WS = sys.argv[1] if len(sys.argv) > 1 else os.environ.get(
    "WS", "/gpfs/bwfor/work/ws/fr_os136-immune_escape")
COHORT = f"{WS}/simulation/cohort"
INDELS = f"{WS}/simulation/indels"
RES = os.path.expanduser("~/immune_escape_project/results")

W = 74
print("=" * W)
print(" 1. WHAT EXISTS")
print("=" * W)

man = pd.read_csv(f"{COHORT}/manifest.tsv", sep="\t")
files = {
    "tumour BAM": "{c}/{s}/{s}_tumor.bam",
    "normal BAM": "{c}/{s}/{s}_normal.bam",
    "raw VCF": "{c}/{s}/mutect2/{s}.raw.vcf.gz",
    "filtered VCF": "{c}/{s}/mutect2/{s}.filtered.vcf.gz",
    "comparison": "{c}/{s}/comparison/truth_vs_calls.tsv",
}
print()
for label, pat in files.items():
    n = sum(os.path.exists(pat.format(c=COHORT, s=s)) for s in man.sample_id)
    mark = "\u2713" if n == len(man) else " "
    print(f"  {mark} {label:<18}{n} of {len(man)}")

n_ind = len(glob.glob(f"{INDELS}/*/calls/*.filtered.vcf.gz"))
print(f"\n  indel calling (step 8b)   {n_ind} of 28 samples done")

print()
print("=" * W)
print(" 2. WHAT RECALL WAS MEASURED AGAINST")
print("=" * W)

comp = []
for s in man.sample_id:
    p = f"{COHORT}/{s}/comparison/truth_vs_calls.tsv"
    if os.path.exists(p):
        t = pd.read_csv(p, sep="\t")
        t["sample"] = s
        comp.append(t)
if not comp:
    sys.exit("\n  no comparison files; step 3 has not been scored")
T = pd.concat(comp, ignore_index=True)

ver_path = f"{RES}/verify_step2_per_mutation.tsv"
if not os.path.exists(ver_path):
    sys.exit(f"\n  missing {ver_path}; rerun verify_step2_injection.py")
v = pd.read_csv(ver_path, sep="\t")
v["alt"] = v.alt.astype(str)
v = v[(v.alt.str.len() == 1) & (v.alt != "-")]
v["key"] = v["sample"] + ":" + v.chrom + ":" + v.pos.astype(str)
present = set(v[v.landed].key)
absent = set(v[~v.landed].key)

pos_col = "Start_Position_hg38" if "Start_Position_hg38" in T.columns else "pos"
chr_col = "Chromosome_hg38" if "Chromosome_hg38" in T.columns else "chrom"
T["key"] = T["sample"] + ":" + T[chr_col] + ":" + T[pos_col].astype(str)

found = int(T.detected_pass.sum())
print(f"\n  rows in the comparison            {len(T)}")
print(f"    of which present in the BAM     {int(T.key.isin(present).sum())}")
print(f"    of which absent from the BAM    {int(T.key.isin(absent).sum())}")
print(f"\n  PASS calls matching the truth set {found}")

r_all = 100 * found / len(T)
Tv = T[T.key.isin(present)]
r_ver = 100 * Tv.detected_pass.mean()
print(f"\n  recall against all {len(T)} rows        {r_all:.1f}%   <- what is on record")
print(f"  recall against the {len(Tv)} verified   {r_ver:.1f}%   <- correct")
print(f"  difference                          {r_ver - r_all:+.1f} points")

ghost = int(T[T.detected_pass].key.isin(absent).sum())
if ghost:
    print(f"\n  PASS calls at positions the check calls absent: {ghost}")
    print(f"  These were the disagreement that exposed the mpileup fault.")
    print(f"  If the check is now right, this should be near zero.")

if "detected_any" in T.columns:
    seen = int(T.detected_any.sum())
    print(f"\n  seen by Mutect2 at all            {seen}"
          f"   ({100*seen/len(Tv):.1f}% of verified)")
    print(f"  seen but filtered out             {seen - found}")

print()
print("=" * W)
print(" 3. IS PRECISION MEASURED")
print("=" * W)
print(f"\n  Recall needs a truth set; precision needs every PASS call,")
print(f"  including ones the truth set says nothing about. A comparison")
print(f"  table built by walking the truth set cannot see those.\n")

extra_col = next((c for c in T.columns
                  if "extra" in c.lower() or "outside" in c.lower()
                  or "false" in c.lower()), None)
if extra_col:
    print(f"  the table carries a column for it: {extra_col}"
          f"   total {int(T[extra_col].sum())}")
else:
    print(f"  no column in the comparison counts calls outside the truth set")

# count PASS records directly out of the VCFs and see how many are unmatched
import subprocess
truth_keys = set(T.key)
checked = 0
extra_total = 0
for s in man.sample_id[:5]:
    vcf = f"{COHORT}/{s}/mutect2/{s}.filtered.vcf.gz"
    if not os.path.exists(vcf):
        continue
    out = subprocess.run(
        f"zcat {vcf} | grep -v '^#' | awk -F'\\t' '$7==\"PASS\"' "
        f"| awk -F'\\t' '{{print $1\":\"$2}}'",
        shell=True, capture_output=True, text=True).stdout.split()
    n_extra = sum(1 for k in out if f"{s}:{k}" not in truth_keys)
    extra_total += n_extra
    checked += 1
    print(f"    {s}   {len(out)} PASS,   {n_extra} outside the truth set")
if checked:
    print(f"\n  across {checked} samples: {extra_total} PASS calls the truth")
    print(f"  set does not contain. Those are the precision denominator.")

print()
print("=" * W)
print(" 4. DOES THE NORMAL BEHAVE AS A CONTROL")
print("=" * W)
if "normal_alt" in v.columns:
    n_contam = int((v.normal_alt > 0).sum())
    print(f"\n  injected positions where the normal also carries the alt")
    print(f"    {n_contam} of {len(v)}   ({100*n_contam/len(v):.1f}%)")
    print(f"\n  A mutation present in the normal is not somatic, and Mutect2")
    print(f"  is right to reject it. These belong in the analysis as")
    print(f"  expected rejections rather than as misses.")
    if n_contam:
        sub = v[v.normal_alt > 0]
        print(f"    median alt reads in the normal   {sub.normal_alt.median():.0f}")
else:
    print(f"\n  the verification does not record the normal; rerun")
    print(f"  verify_step2_injection.py to add it")

print()
print("=" * W)
print(" 5. WHAT IS NOT COVERED")
print("=" * W)
gaps = []
if not extra_col and not checked:
    gaps.append("precision has no denominator anywhere in the results")
if n_ind < 28:
    gaps.append(f"indel calling incomplete ({n_ind} of 28)")
if abs(r_ver - r_all) > 0.5:
    gaps.append("recall on record uses the pre-correction truth set")
if "detected_any" not in T.columns:
    gaps.append("no record of variants seen but filtered, so the filter "
                "sweep has nothing to sit on")

if gaps:
    for g in gaps:
        print(f"\n  - {g}")
else:
    print(f"\n  nothing outstanding")

print(f"\n  Not gaps, but worth stating on the slide:")
print(f"    the sensitivity curve is binned by requested fraction, and")
print(f"    the injection lands at about 0.94 of it")
print(f"    31 substitutions are absent from the files and 18 of those")
print(f"    are in the MHC")
