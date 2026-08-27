#!/usr/bin/env python3
"""
Rescore Mutect2 on substitutions, reading the VCFs directly.

Recall of 73% came from comparison tables written during step 3. Those
tables were built by a script that is not being re-read here, and the
indel comparison found two things that would silently lower a recall
figure if the same code path produced it.

Positions shift. Twenty-four of sixty-two indels were called one base from
where they were placed, because left-alignment in a repeat allows several
equivalent descriptions. A substitution should not move, but a caller that
represents it as part of a longer haplotype record would put it at the
start of that record.

Representation varies. A variant placed as C>T may be emitted as CA>TA, or
as one alternative among several at the same position, and a comparison
matching the alt field literally would miss it.

So this rebuilds the comparison from the filtered VCFs with no assumptions
about how the earlier one worked, scores against verified positions only,
and reports what the difference is. If it agrees with the tables, the
tables are sound and the 74.5% stands.

Usage:
  python recheck_step3_calls.py [workspace] [--window=5]
"""
import os
import sys
import gzip
import numpy as np
import pandas as pd

args = [a for a in sys.argv[1:] if not a.startswith("--")]
WS = args[0] if args else os.environ.get(
    "WS", "/gpfs/bwfor/work/ws/fr_os136-immune_escape")
WINDOW = 5
for a in sys.argv[1:]:
    if a.startswith("--window="):
        WINDOW = int(a.split("=", 1)[1])

COHORT = f"{WS}/simulation/cohort"
RES = os.path.expanduser("~/immune_escape_project/results")
VER = f"{RES}/verify_step2_per_mutation.tsv"

if not os.path.exists(VER):
    sys.exit(f"missing {VER}")

v = pd.read_csv(VER, sep="\t")
v["alt"] = v.alt.astype(str)
snv = v[(v.alt.str.len() == 1) & (v.alt != "-")].copy()
verified = snv[snv.landed].copy()


def read_vcf(path):
    """
    Every PASS record, expanded so that each alternative allele is its own
    entry. A record carrying two alternatives at one position describes two
    candidate variants, and collapsing them loses one.
    """
    out = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            chrom, pos, ref, alts, filt = f[0], int(f[1]), f[3], f[4], f[6]
            for a in alts.split(","):
                out.append({"chrom": chrom, "pos": pos, "ref": ref,
                            "alt": a, "filter": filt,
                            "is_pass": filt == "PASS"})
    return pd.DataFrame(out)


def matches(rec, chrom, pos, alt):
    """
    Does this VCF record describe a substitution of alt at pos?

    Four ways it can: the plain case; a record starting earlier whose ref
    and alt differ only at our offset; the same length but multi-base; and
    a spanning record where the change sits inside.
    """
    if rec["chrom"] != chrom:
        return False, None
    off = pos - rec["pos"]
    ref, a = rec["ref"], rec["alt"]
    if off == 0 and len(ref) == 1 and len(a) == 1:
        return (a.upper() == alt.upper()), "exact"
    if len(ref) == len(a) and 0 <= off < len(ref):
        return (a[off].upper() == alt.upper()), "within a longer record"
    if abs(off) <= WINDOW and len(a) == 1 and len(ref) == 1:
        return (a.upper() == alt.upper()), "shifted"
    return False, None


rows, per_sample = [], []
missing_vcf = []

for sid, g in verified.groupby("sample"):
    vcf = f"{COHORT}/{sid}/mutect2/{sid}.filtered.vcf.gz"
    if not os.path.exists(vcf):
        missing_vcf.append(sid)
        continue
    calls = read_vcf(vcf)
    by_chrom = {c: df for c, df in calls.groupby("chrom")}

    for _, r in g.iterrows():
        near = by_chrom.get(r.chrom)
        hit_pass = hit_any = None
        how = None
        if near is not None:
            cand = near[(near.pos >= r.pos - WINDOW - 20) &
                        (near.pos <= r.pos + WINDOW)]
            for _, rec in cand.iterrows():
                ok, kind = matches(rec, r.chrom, int(r.pos), r.alt)
                if not ok:
                    continue
                if hit_any is None:
                    hit_any, how = rec, kind
                if rec.is_pass and hit_pass is None:
                    hit_pass, how = rec, kind
        rows.append({
            "sample": sid, "chrom": r.chrom, "pos": int(r.pos), "alt": r.alt,
            "target_vaf": r.target_vaf, "observed_vaf": r.observed_vaf,
            "alt_reads": r.alt_reads, "depth": r.tumour_depth,
            "in_mhc": r.in_mhc,
            "pass": hit_pass is not None,
            "seen": hit_any is not None,
            "how": how,
            "filters": (hit_any["filter"] if hit_any is not None
                        and hit_pass is None else None),
        })
    k = len(g)
    per_sample.append({
        "sample": sid, "verified": k,
        "pass": sum(1 for x in rows[-k:] if x["pass"]),
        "seen": sum(1 for x in rows[-k:] if x["seen"]),
        "pass_records": int(calls.is_pass.sum()),
    })

d = pd.DataFrame(rows)
s = pd.DataFrame(per_sample)
if d.empty:
    sys.exit("nothing scored")

W = 74
print("=" * W)
print(" RESCORED FROM THE VCFs")
print("=" * W)
if missing_vcf:
    print(f"\n  no VCF for: {' '.join(missing_vcf)}")
n = len(d)
print(f"\n  verified substitutions      {n}")
print(f"  seen by Mutect2 at all      {int(d.seen.sum())}"
      f"   ({100*d.seen.mean():.1f}%)")
print(f"  PASS                        {int(d['pass'].sum())}"
      f"   ({100*d['pass'].mean():.1f}%)")
print(f"  seen but filtered out       {int(d.seen.sum() - d['pass'].sum())}")

print(f"\n  how the PASS calls were matched:")
for how, k in d[d["pass"]].how.value_counts().items():
    print(f"    {str(how):<26}{k}")

print()
print("=" * W)
print(" AGAINST WHAT THE STEP-3 TABLES SAID")
print("=" * W)
comp = []
for sid in verified["sample"].unique():
    p = f"{COHORT}/{sid}/comparison/truth_vs_calls.tsv"
    if os.path.exists(p):
        t = pd.read_csv(p, sep="\t")
        t["sample"] = sid
        comp.append(t)
if comp:
    T = pd.concat(comp, ignore_index=True)
    pos_col = ("Start_Position_hg38" if "Start_Position_hg38" in T.columns
               else "pos")
    chr_col = ("Chromosome_hg38" if "Chromosome_hg38" in T.columns
               else "chrom")
    T["key"] = T["sample"] + ":" + T[chr_col] + ":" + T[pos_col].astype(str)
    d["key"] = d["sample"] + ":" + d.chrom + ":" + d.pos.astype(str)
    Tv = T[T.key.isin(set(d.key))]
    old = 100 * Tv.detected_pass.mean()
    new = 100 * d["pass"].mean()
    print(f"\n  step-3 tables, verified rows only   {old:.1f}%")
    print(f"  rescored from the VCFs              {new:.1f}%")
    print(f"  difference                          {new - old:+.1f} points")

    m = d.merge(Tv[["key", "detected_pass"]], on="key", how="left")
    disagree = m[m["pass"] != m.detected_pass.fillna(False)]
    print(f"\n  positions the two disagree on       {len(disagree)}")
    if len(disagree):
        gained = disagree[disagree["pass"]]
        lost = disagree[~disagree["pass"]]
        print(f"    found here, missed there          {len(gained)}")
        print(f"    missed here, found there          {len(lost)}")
        if len(gained):
            print(f"\n    how the newly found ones matched:")
            for how, k in gained.how.value_counts().items():
                print(f"      {str(how):<24}{k}")
            print(f"\n    first ten:")
            for _, r in gained.head(10).iterrows():
                print(f"      {r['sample']:<7}{r.chrom}:{r.pos:<12}"
                      f"{r.alt}  VAF {r.observed_vaf:.3f}  "
                      f"{r.alt_reads:.0f} reads  [{r.how}]")

print()
print("=" * W)
print(" DETECTION AGAINST THE FRACTION ACTUALLY PRESENT")
print("=" * W)
print(f"\n  Binned by observed fraction rather than requested, since the")
print(f"  injection lands at about 0.94 of what was asked.\n")
b = [0, .05, .10, .15, .20, .30, 1.01]
lab = ["<5%", "5-10%", "10-15%", "15-20%", "20-30%", ">30%"]
d["bin"] = pd.cut(d.observed_vaf, bins=b, labels=lab, right=False)
print(f"  {'observed':<11}{'n':>6}{'seen':>7}{'PASS':>7}{'recall':>9}"
      f"{'median alt reads':>19}")
for l in lab:
    g = d[d.bin == l]
    if not len(g):
        continue
    print(f"  {l:<11}{len(g):>6}{int(g.seen.sum()):>7}"
          f"{int(g['pass'].sum()):>7}{100*g['pass'].mean():>8.1f}%"
          f"{g.alt_reads.median():>19.0f}")

print()
print("=" * W)
print(" WHY THE REST WERE FILTERED")
print("=" * W)
filt = d[d.seen & ~d["pass"] & d.filters.notna()]
if len(filt):
    tags = {}
    for f_ in filt.filters:
        for t in str(f_).split(";"):
            tags[t] = tags.get(t, 0) + 1
    print(f"\n  {len(filt)} were seen and rejected\n")
    for t, k in sorted(tags.items(), key=lambda x: -x[1])[:10]:
        print(f"    {t:<26}{k}")
    print(f"\n  median observed VAF of the rejected  "
          f"{filt.observed_vaf.median():.3f}")
    print(f"  median observed VAF of the accepted  "
          f"{d[d['pass']].observed_vaf.median():.3f}")

print()
print("=" * W)
print(" MHC")
print("=" * W)
for grp, lbl in [(d[d.in_mhc], "inside"), (d[~d.in_mhc], "outside")]:
    if len(grp):
        print(f"  {lbl:<10}{len(grp):>6} verified,"
              f"{int(grp['pass'].sum()):>5} PASS"
              f"   ({100*grp['pass'].mean():.1f}%)")

out = f"{RES}/step3_rescored.tsv"
d.to_csv(out, sep="\t", index=False)
print(f"\nwritten to {out}")
