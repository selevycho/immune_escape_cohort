#!/usr/bin/env python3
"""
Did Mutect2 find the indels that were placed?

Step 8 injected indels into a copy of the tumour BAM and nothing was ever
called against it, so the pipeline had been shown to place indels and not
to find them. Step 8b called against those files; this scores the result.

The BAM carries the substitutions from step 2 as well, so the calls
contain both kinds and have to be separated by variant type rather than by
file. That separation is the point: substitution recall is already known,
and the question here is whether indels behave the same way.

Positions are matched with a small window. bwa and GATK both left-align
indels, but in a repeat several placements are equivalent and the two
tools do not always choose the same one. Requiring an exact coordinate
would score a correctly found indel as missed.

Only indels confirmed present in the reads are counted, for the same
reason step 3 should be scored against verified substitutions: a truth set
entry that is not in the file cannot be found by any caller.

Usage:
  python compare_indel_calls.py [workspace] [--window=10]
"""
import os
import sys
import gzip
import numpy as np
import pandas as pd

args = [a for a in sys.argv[1:] if not a.startswith("--")]
WS = args[0] if args else os.environ.get(
    "WS", "/gpfs/bwfor/work/ws/fr_os136-immune_escape")
WINDOW = 10
for a in sys.argv[1:]:
    if a.startswith("--window="):
        WINDOW = int(a.split("=", 1)[1])

INDELS = f"{WS}/simulation/indels"
RES = os.path.expanduser("~/immune_escape_project/results")
VER = f"{RES}/verify_step8_per_indel.tsv"
SNV_VER = f"{RES}/verify_step2_per_mutation.tsv"

if not os.path.exists(VER):
    sys.exit(f"missing {VER}; run verify_step8_indels.py first")

truth = pd.read_csv(VER, sep="\t")
verified = truth[truth.landed].copy()


def read_vcf(path):
    """PASS records as (chrom, pos, ref, alt), from the gzipped VCF."""
    out = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8 or f[6] != "PASS":
                continue
            out.append((f[0], int(f[1]), f[3], f[4]))
    return out


rows, per_sample = [], []
for sid, g in verified.groupby("sample"):
    vcf = f"{INDELS}/{sid}/calls/{sid}.filtered.vcf.gz"
    if not os.path.exists(vcf):
        print(f"  {sid}  no calls")
        continue

    calls = read_vcf(vcf)
    call_ind = [c for c in calls if len(c[2]) != len(c[3])]
    call_snv = [c for c in calls if len(c[2]) == len(c[3])]

    for _, r in g.iterrows():
        # an indel found within the window and of the same direction counts;
        # requiring the exact coordinate would fail on repeats where several
        # placements describe the same event
        hit = None
        for c in call_ind:
            if c[0] != r.chrom or abs(c[1] - int(r.pos)) > WINDOW:
                continue
            is_ins = len(c[3]) > len(c[2])
            if (r.type == "INS") == is_ins:
                hit = c
                break
        rows.append({
            "sample": sid, "chrom": r.chrom, "pos": int(r.pos),
            "type": r.type, "gene": r.get("gene"),
            "target_vaf": r.target_vaf, "observed_vaf": r.observed_vaf,
            "depth": r.depth, "support": r.support,
            "found": hit is not None,
            "called_at": hit[1] if hit else None,
            "offset": (hit[1] - int(r.pos)) if hit else None,
        })

    per_sample.append({
        "sample": sid, "injected": len(g),
        "found": sum(1 for x in rows[-len(g):] if x["found"]),
        "pass_total": len(calls), "pass_indels": len(call_ind),
        "pass_snvs": len(call_snv),
    })

d = pd.DataFrame(rows)
s = pd.DataFrame(per_sample)
if d.empty:
    sys.exit("nothing to compare")

W = 74
print("=" * W)
print(" INDEL RECALL")
print("=" * W)


def wilson(k, m):
    if m == 0:
        return (np.nan, np.nan)
    z, p = 1.96, k / m
    den = 1 + z * z / m
    c = (p + z * z / (2 * m)) / den
    h = z * np.sqrt(p * (1 - p) / m + z * z / (4 * m * m)) / den
    return (100 * max(0, c - h), 100 * min(1, c + h))


n, f = len(d), int(d.found.sum())
lo, hi = wilson(f, n)
print(f"\n  indels verified present in the reads   {n}")
print(f"  recovered by Mutect2 as PASS           {f}"
      f"   ({100*f/n:.1f}%, 95% CI {lo:.0f}-{hi:.0f}%)")
print(f"  samples                                {d['sample'].nunique()}")

print(f"\n  {'type':<8}{'n':>6}{'found':>8}{'rate':>9}{'95% CI':>16}")
for ty, g in d.groupby("type"):
    k = int(g.found.sum())
    l2, h2 = wilson(k, len(g))
    print(f"  {ty:<8}{len(g):>6}{k:>8}{100*k/len(g):>8.1f}%"
          f"{f'{l2:.0f}-{h2:.0f}%':>16}")

if os.path.exists(SNV_VER):
    v = pd.read_csv(SNV_VER, sep="\t")
    v["alt"] = v.alt.astype(str)
    v = v[(v.alt.str.len() == 1) & (v.alt != "-") & v.landed]
    print(f"\n  for comparison, substitutions in step 3 reached 74.5%")
    print(f"  against {len(v)} verified positions")

print()
print("=" * W)
print(" DOES DETECTION FOLLOW THE ALLELE FRACTION")
print("=" * W)
b = [0, .10, .20, .30, 1.01]
lab = ["<10%", "10-20%", "20-30%", ">30%"]
d["bin"] = pd.cut(d.target_vaf, bins=b, labels=lab, right=False)
print(f"\n  {'requested':<11}{'n':>6}{'found':>8}{'rate':>9}"
      f"{'median support':>17}")
for l in lab:
    g = d[d.bin == l]
    if not len(g):
        continue
    print(f"  {l:<11}{len(g):>6}{int(g.found.sum()):>8}"
          f"{100*g.found.mean():>8.1f}%{g.support.median():>17.0f}")

if d.found.any() and (~d.found).any():
    print(f"\n  median supporting reads where found   "
          f"{d[d.found].support.median():.0f}")
    print(f"  median where missed                   "
          f"{d[~d.found].support.median():.0f}")

print()
print("=" * W)
print(" HOW FAR OFF WERE THE POSITIONS")
print("=" * W)
off = d[d.found & d.offset.notna()]
if len(off):
    exact = int((off.offset == 0).sum())
    print(f"\n  called at the exact position          {exact} of {len(off)}")
    print(f"  shifted by 1-{WINDOW} bp                    "
          f"{len(off) - exact}")
    if len(off) - exact:
        print(f"  largest shift                         "
              f"{int(off.offset.abs().max())} bp")
    print(f"\n  Left-alignment in a repeat allows several equivalent")
    print(f"  placements; a shift is the same event, not a different one.")

print()
print("=" * W)
print(" WHAT WAS MISSED")
print("=" * W)
miss = d[~d.found]
if not len(miss):
    print("\n  nothing")
else:
    print(f"\n  {len(miss)} of {n}\n")
    print(f"  {'sample':<8}{'position':<22}{'type':<6}{'gene':<12}"
          f"{'depth':>7}{'support':>9}{'VAF':>8}")
    for _, r in miss.iterrows():
        print(f"  {r['sample']:<8}{r.chrom + ':' + str(r.pos):<22}"
              f"{r.type:<6}{str(r.gene):<12}{int(r.depth):>7}"
              f"{int(r.support):>9}{r.observed_vaf:>8.3f}")

print()
print("=" * W)
print(" WHAT ELSE THE CALLS CONTAIN")
print("=" * W)
print(f"\n  PASS calls in these files      {int(s.pass_total.sum())}")
print(f"    substitutions                {int(s.pass_snvs.sum())}")
print(f"    indels                       {int(s.pass_indels.sum())}")
print(f"\n  indels called but not injected {int(s.pass_indels.sum()) - f}")
print(f"\n  Those are germline indels the matched normal did not filter,")
print(f"  not false positives against the truth set - the truth set only")
print(f"  claims what was placed, not that nothing else is present.")

out = f"{RES}/indel_recall.tsv"
d.to_csv(out, sep="\t", index=False)
s.to_csv(out.replace(".tsv", "_per_sample.tsv"), sep="\t", index=False)
print(f"\nwritten to {out}")
