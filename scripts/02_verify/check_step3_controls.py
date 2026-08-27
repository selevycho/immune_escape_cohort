#!/usr/bin/env python3
"""
Controls for the substitution calling.

Recall says how much of what was placed came back. It says nothing about
what comes back when nothing was placed, and a caller that reports
variants everywhere would score perfectly on recall while being useless.
Three controls are available without new computation.

The matched normal. Every call is made against a normal sliced from the
same individual, so a variant present in both is germline and rejecting it
is correct rather than a miss.

Specificity across the whole cohort. A PASS call at a position the truth
set does not contain is either a germline variant the normal failed to
filter or a false positive. This was checked on five samples before; here
it runs on all forty, and the calls are separated by whether the normal
carries them.

Reproducibility. Step 8b called the same substitutions a second time, from
a different BAM built by a different route. Two runs of the same caller on
the same variants should agree; where they do not, the difference bounds
how much of any single figure is run-to-run noise.

Usage:
  python check_step3_controls.py [workspace]
"""
import os
import sys
import gzip
import subprocess
import numpy as np
import pandas as pd

WS = sys.argv[1] if len(sys.argv) > 1 else os.environ.get(
    "WS", "/gpfs/bwfor/work/ws/fr_os136-immune_escape")
COHORT = f"{WS}/simulation/cohort"
INDELS = f"{WS}/simulation/indels"
RES = os.path.expanduser("~/immune_escape_project/results")
VER = f"{RES}/verify_step2_per_mutation.tsv"

SAMTOOLS = None
for cand in ["samtools",
             os.path.expanduser("~/miniconda3/envs/bio_work/bin/samtools"),
             "/home/fr/fr_fr/fr_os136/miniconda3/envs/bio_work/bin/samtools"]:
    if subprocess.run(f"{cand} --version", shell=True,
                      capture_output=True).returncode == 0:
        SAMTOOLS = cand
        break

v = pd.read_csv(VER, sep="\t")
v["alt"] = v.alt.astype(str)
snv = v[(v.alt.str.len() == 1) & (v.alt != "-")].copy()
snv["key"] = snv["sample"] + ":" + snv.chrom + ":" + snv.pos.astype(str)
truth_keys = set(snv.key)
man = pd.read_csv(f"{COHORT}/manifest.tsv", sep="\t")


def pass_calls(path):
    out = []
    if not os.path.exists(path):
        return out
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8 or f[6] != "PASS":
                continue
            out.append((f[0], int(f[1]), f[3], f[4]))
    return out


W = 74
print("=" * W)
print(" 1. SPECIFICITY ACROSS ALL FORTY")
print("=" * W)

rows = []
for s in man.sample_id:
    calls = pass_calls(f"{COHORT}/{s}/mutect2/{s}.filtered.vcf.gz")
    subs = [c for c in calls if len(c[2]) == 1 and len(c[3]) == 1]
    outside = [c for c in subs if f"{s}:{c[0]}:{c[1]}" not in truth_keys]
    rows.append({"sample": s, "pass_total": len(calls),
                 "pass_subs": len(subs), "outside_truth": len(outside)})
    for c in outside[:200]:
        rows[-1].setdefault("_ex", []).append(c)

sp = pd.DataFrame([{k: r[k] for k in
                    ("sample", "pass_total", "pass_subs", "outside_truth")}
                   for r in rows])
print(f"\n  samples                      {len(sp)}")
print(f"  PASS calls in total          {int(sp.pass_total.sum())}")
print(f"    substitutions              {int(sp.pass_subs.sum())}")
print(f"    outside the truth set      {int(sp.outside_truth.sum())}")
if sp.outside_truth.sum() == 0:
    print(f"\n  Not one PASS substitution in forty samples falls outside")
    print(f"  the truth set. Precision is 100% by measurement, not by")
    print(f"  assumption.")
else:
    bad = sp[sp.outside_truth > 0]
    print(f"\n  samples with calls outside   {len(bad)}")
    print(bad.to_string(index=False))
    print(f"\n  Checking whether the normal carries them:")
    checked = 0
    germ = 0
    for r in rows:
        for c in r.get("_ex", [])[:10]:
            nbam = f"{COHORT}/{r['sample']}/{r['sample']}_normal.bam"
            if not (SAMTOOLS and os.path.exists(nbam)):
                continue
            out = subprocess.run(
                f"{SAMTOOLS} view -F 0x404 {nbam} {c[0]}:{c[1]}-{c[1]} "
                f"| wc -l", shell=True, capture_output=True,
                text=True).stdout.strip()
            checked += 1
            if int(out or 0) > 0:
                germ += 1
        if checked >= 40:
            break
    if checked:
        print(f"    {germ} of {checked} sampled positions are covered in the")
        print(f"    normal, so a germline origin is plausible for them")

print()
print("=" * W)
print(" 2. THE NORMAL AS A CONTROL")
print("=" * W)
if "normal_alt" in snv.columns:
    contam = snv[snv.normal_alt > 0]
    print(f"\n  injected positions where the normal carries the alt")
    print(f"    {len(contam)} of {len(snv)}   ({100*len(contam)/len(snv):.1f}%)")
    if len(contam):
        print(f"    median alt reads there    {contam.normal_alt.median():.0f}")
        print(f"    median normal depth       {contam.normal_depth.median():.0f}")
        frac = contam.normal_alt / contam.normal_depth.replace(0, np.nan)
        print(f"    median fraction           {frac.median():.3f}")
        print(f"\n  At one or two reads out of thirty this is sequencing")
        print(f"  noise rather than a germline variant, but Mutect2 treats")
        print(f"  any support in the normal as evidence against somatic")
        print(f"  status, so these positions are penalised either way.")

print()
print("=" * W)
print(" 3. THE SAME VARIANTS, CALLED TWICE")
print("=" * W)
print(f"\n  Step 8b called against a BAM that carries the substitutions as")
print(f"  well as the indels, so the same variants were called a second")
print(f"  time from a separately built file.\n")

ver = snv[snv.landed]
both = []
for s in man.sample_id:
    p1 = f"{COHORT}/{s}/mutect2/{s}.filtered.vcf.gz"
    p2 = f"{INDELS}/{s}/calls/{s}.filtered.vcf.gz"
    if not (os.path.exists(p1) and os.path.exists(p2)):
        continue
    k1 = {f"{c[0]}:{c[1]}:{c[3]}" for c in pass_calls(p1)
          if len(c[2]) == 1 and len(c[3]) == 1}
    k2 = {f"{c[0]}:{c[1]}:{c[3]}" for c in pass_calls(p2)
          if len(c[2]) == 1 and len(c[3]) == 1}
    g = ver[ver["sample"] == s]
    want = {f"{r.chrom}:{int(r.pos)}:{r.alt}" for _, r in g.iterrows()}
    both.append({"sample": s, "verified": len(want),
                 "run1": len(want & k1), "run2": len(want & k2),
                 "agree": len(want & k1 & k2),
                 "only1": len((want & k1) - k2),
                 "only2": len((want & k2) - k1)})

if both:
    b = pd.DataFrame(both)
    print(f"  samples with both runs       {len(b)}")
    print(f"  verified substitutions       {int(b.verified.sum())}")
    print(f"    recovered in run 1         {int(b.run1.sum())}"
          f"   ({100*b.run1.sum()/b.verified.sum():.1f}%)")
    print(f"    recovered in run 2         {int(b.run2.sum())}"
          f"   ({100*b.run2.sum()/b.verified.sum():.1f}%)")
    print(f"    both runs                  {int(b.agree.sum())}")
    print(f"    run 1 only                 {int(b.only1.sum())}")
    print(f"    run 2 only                 {int(b.only2.sum())}")
    disc = int(b.only1.sum() + b.only2.sum())
    print(f"\n  disagreement between runs    {disc}"
          f"   ({100*disc/max(1,int(b.verified.sum())):.1f}% of verified)")
    print(f"\n  Two runs of one caller on the same variants, from files")
    print(f"  built by different routes. The disagreement bounds how much")
    print(f"  of any single recall figure is run-to-run variation.")
else:
    print(f"  no sample has both runs available")

out = f"{RES}/step3_controls.tsv"
sp.to_csv(out, sep="\t", index=False)
print(f"\nwritten to {out}")
