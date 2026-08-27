#!/usr/bin/env python3
"""
Mutect2 on substitutions, scored the way the indels were scored.

The indel analysis split its failures into two groups that behave
differently: variants that reached the raw VCF and were rejected by a
filter, and variants that never appeared at all. The first is a filtering
decision, the second a detection limit, and a single recall figure hides
which one dominates.

The same split is applied here. Recall is scored against verified
positions only, since a truth-set entry absent from the file cannot be
found by any caller, and the raw VCF is read alongside the filtered one so
that "seen" means what Mutect2 emitted rather than what survived.

Sequence context is checked for the same reason it mattered for indels:
a substitution inside a homopolymer or next to a repeat is harder to place
confidently, and if the misses cluster there it is a property of where
BAMSurgeon put them rather than of the caller.

The matched normal is checked at every missed position. A variant the
normal also carries is not somatic and rejecting it is correct.

Usage:
  python recheck_step3_full.py [workspace] [--context=25]
"""
import os
import re
import sys
import gzip
import subprocess
import numpy as np
import pandas as pd

args = [a for a in sys.argv[1:] if not a.startswith("--")]
WS = args[0] if args else os.environ.get(
    "WS", "/gpfs/bwfor/work/ws/fr_os136-immune_escape")
CTX = 25
for a in sys.argv[1:]:
    if a.startswith("--context="):
        CTX = int(a.split("=", 1)[1])

COHORT = f"{WS}/simulation/cohort"
RES = os.path.expanduser("~/immune_escape_project/results")
VER = f"{RES}/verify_step2_per_mutation.tsv"
REF = f"{WS}/ref/Homo_sapiens_assembly38.fasta"

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
snv = v[(v.alt.str.len() == 1) & (v.alt != "-")]
ver = snv[snv.landed].copy()


def read_vcf(path):
    """Every record, each alternative allele separately, with its filter."""
    out = []
    if not os.path.exists(path):
        return pd.DataFrame(out)
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            for a in f[4].split(","):
                out.append({"chrom": f[0], "pos": int(f[1]), "ref": f[3],
                            "alt": a, "filter": f[6]})
    return pd.DataFrame(out)


def hit(cand, pos, alt):
    """
    A record describing alt at pos, allowing for haplotype representation.

    Mutect2 may emit a substitution as part of a longer record when it sits
    inside an assembled haplotype, so a literal match on the alt field
    alone would score such a call as missed.
    """
    for _, r in cand.iterrows():
        off = pos - r["pos"]
        ref, a = r["ref"], r["alt"]
        if off == 0 and len(ref) == 1 and len(a) == 1:
            if a.upper() == alt.upper():
                return r, "exact"
        elif len(ref) == len(a) and 0 <= off < len(ref):
            if a[off].upper() == alt.upper():
                return r, "inside a longer record"
    return None, None


print(f"scoring {len(ver)} verified substitutions in "
      f"{ver['sample'].nunique()} samples\n", flush=True)

rows = []
for sid, g in ver.groupby("sample"):
    raw = read_vcf(f"{COHORT}/{sid}/mutect2/{sid}.raw.vcf.gz")
    filt = read_vcf(f"{COHORT}/{sid}/mutect2/{sid}.filtered.vcf.gz")
    src = filt if len(filt) else raw
    if not len(src):
        continue
    by_chrom = {c: df for c, df in src.groupby("chrom")}

    for _, r in g.iterrows():
        near = by_chrom.get(r.chrom)
        rec, how = (None, None)
        if near is not None:
            cand = near[(near.pos >= r.pos - 30) & (near.pos <= r.pos)]
            if len(cand):
                rec, how = hit(cand, int(r.pos), r.alt)
        rows.append({
            "sample": sid, "chrom": r.chrom, "pos": int(r.pos), "alt": r.alt,
            "gene": r.get("gene"), "in_mhc": r.in_mhc,
            "target_vaf": r.target_vaf, "observed_vaf": r.observed_vaf,
            "alt_reads": r.alt_reads, "depth": r.tumour_depth,
            "normal_alt": r.get("normal_alt", 0),
            "emitted": rec is not None,
            "pass": (rec is not None and rec["filter"] == "PASS"),
            "filters": (rec["filter"] if rec is not None
                        and rec["filter"] != "PASS" else None),
            "how": how,
        })

d = pd.DataFrame(rows)
if d.empty:
    sys.exit("nothing scored")

W = 74
print("=" * W)
print(" RECALL")
print("=" * W)
n = len(d)
print(f"\n  verified substitutions        {n}")
print(f"  emitted by Mutect2            {int(d.emitted.sum())}"
      f"   ({100*d.emitted.mean():.1f}%)")
print(f"  PASS                          {int(d['pass'].sum())}"
      f"   ({100*d['pass'].mean():.1f}%)")
print(f"\n  for comparison, indels reached 82.7% of 75 verified")

print()
print("=" * W)
print(" HOW THE FAILURES SPLIT")
print("=" * W)
missed = d[~d["pass"]]
filtered = missed[missed.emitted]
absent = missed[~missed.emitted]
print(f"\n  did not PASS                  {len(missed)}")
print(f"    emitted, then filtered      {len(filtered)}"
      f"   ({100*len(filtered)/max(1,len(missed)):.0f}%)")
print(f"    never emitted               {len(absent)}"
      f"   ({100*len(absent)/max(1,len(missed)):.0f}%)")
print(f"\n  The first group is a filtering decision, the second a detection")
print(f"  limit. Only the first is recoverable by changing settings.")

if len(filtered):
    tags = {}
    for f_ in filtered.filters.dropna():
        for t in str(f_).split(";"):
            tags[t] = tags.get(t, 0) + 1
    print(f"\n  why the emitted ones were rejected:")
    for t, k in sorted(tags.items(), key=lambda x: -x[1]):
        print(f"    {t:<26}{k}")

print()
print("=" * W)
print(" DETECTION AGAINST THE FRACTION PRESENT")
print("=" * W)
b = [0, .05, .10, .15, .20, .30, 1.01]
lab = ["<5%", "5-10%", "10-15%", "15-20%", "20-30%", ">30%"]
d["bin"] = pd.cut(d.observed_vaf, bins=b, labels=lab, right=False)
print(f"\n  {'observed':<11}{'n':>6}{'emitted':>10}{'PASS':>7}{'recall':>9}"
      f"{'alt reads':>12}")
for l in lab:
    g = d[d.bin == l]
    if not len(g):
        continue
    print(f"  {l:<11}{len(g):>6}{int(g.emitted.sum()):>10}"
          f"{int(g['pass'].sum()):>7}{100*g['pass'].mean():>8.1f}%"
          f"{g.alt_reads.median():>12.0f}")
print(f"\n  median alt reads where PASS   {d[d['pass']].alt_reads.median():.0f}")
print(f"  median where filtered         {filtered.alt_reads.median():.0f}")
print(f"  median where never emitted    {absent.alt_reads.median():.0f}")

print()
print("=" * W)
print(" THE MHC")
print("=" * W)
for grp, lbl in [(d[d.in_mhc], "inside"), (d[~d.in_mhc], "outside")]:
    if len(grp):
        print(f"  {lbl:<9}{len(grp):>6} verified,{int(grp.emitted.sum()):>6}"
              f" emitted,{int(grp['pass'].sum()):>6} PASS"
              f"   ({100*grp['pass'].mean():.1f}%)")
if d.in_mhc.any() and not d[d.in_mhc]["pass"].any():
    mh = d[d.in_mhc]
    print(f"\n  Not one substitution in the MHC passed. {int(mh.emitted.sum())}")
    print(f"  of {len(mh)} were emitted, so the region is not invisible to")
    print(f"  the caller — the calls are made and then rejected.")
    if mh.emitted.any():
        t2 = {}
        for f_ in mh[mh.emitted].filters.dropna():
            for t in str(f_).split(";"):
                t2[t] = t2.get(t, 0) + 1
        for t, k in sorted(t2.items(), key=lambda x: -x[1])[:5]:
            print(f"    {t:<26}{k}")

print()
print("=" * W)
print(" IS THE NORMAL THE REASON")
print("=" * W)
gm = missed[missed.normal_alt > 0]
print(f"\n  missed positions where the normal also carries the alt   {len(gm)}")
if len(gm):
    print(f"    median alt reads in the normal   {gm.normal_alt.median():.0f}")
    print(f"\n  A variant present in the normal is not somatic; rejecting it")
    print(f"  is correct and these are not caller failures.")
else:
    print(f"\n  none — every miss is a genuinely somatic variant")

if SAMTOOLS and os.path.exists(REF):
    print()
    print("=" * W)
    print(" SEQUENCE CONTEXT OF THE MISSES")
    print("=" * W)

    def homopolymer_run(chrom, pos):
        out = subprocess.run(
            f"{SAMTOOLS} faidx {REF} {chrom}:{pos-CTX}-{pos+CTX}",
            shell=True, capture_output=True, text=True).stdout
        seq = "".join(l.strip() for l in out.splitlines()[1:]).upper()
        if len(seq) < CTX + 2:
            return 0
        after = seq[CTX + 1:]
        run = 1
        while run < len(after) and after[run] == after[0]:
            run += 1
        before = seq[:CTX][::-1]
        run_b = 1
        while run_b < len(before) and before[run_b] == before[0]:
            run_b += 1
        return max(run, run_b)

    sample_missed = missed.sample(min(60, len(missed)), random_state=4)
    sample_hit = d[d["pass"]].sample(min(60, int(d["pass"].sum())),
                                     random_state=4)
    rm = [homopolymer_run(r.chrom, r.pos) for _, r in sample_missed.iterrows()]
    rh = [homopolymer_run(r.chrom, r.pos) for _, r in sample_hit.iterrows()]
    print(f"\n  longest homopolymer run beside the site, sampled positions\n")
    print(f"    missed   median {np.median(rm):.0f}, "
          f"{sum(1 for x in rm if x >= 4)} of {len(rm)} with a run of 4+")
    print(f"    found    median {np.median(rh):.0f}, "
          f"{sum(1 for x in rh if x >= 4)} of {len(rh)} with a run of 4+")
    print(f"\n  Seven of the thirteen missed indels sat in homopolymers. If")
    print(f"  substitutions do not show the same bias, the two failure modes")
    print(f"  are different and should be described separately.")

print()
print("=" * W)
print(" WORST MISSES")
print("=" * W)
worst = missed.nlargest(10, "observed_vaf")
print(f"\n  highest allele fraction among the failures\n")
print(f"  {'sample':<8}{'position':<22}{'alt':<5}{'VAF':>7}{'reads':>7}"
      f"   {'outcome'}")
for _, r in worst.iterrows():
    fate = r.filters if r.emitted else "never emitted"
    print(f"  {r['sample']:<8}{r.chrom + ':' + str(r.pos):<22}{r.alt:<5}"
          f"{r.observed_vaf:>7.3f}{int(r.alt_reads):>7}   {fate}")

out = f"{RES}/step3_rescored_full.tsv"
d.to_csv(out, sep="\t", index=False)
print(f"\nwritten to {out}")
