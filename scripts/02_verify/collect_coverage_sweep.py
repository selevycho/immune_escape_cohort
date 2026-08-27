#!/usr/bin/env python3
"""
The coverage sweep, scored the way step 3 was.

Both BAMs were downsampled to five depths and called again. The question
is the one an experimenter actually faces: how much sequencing is enough.

Recall is measured against mutations verified present in the reads at
full depth, not against the truth set as written. A mutation that was
never in the file cannot be found at any depth, and counting it as a miss
would make every level look worse by the same amount.

Substitutions are matched allowing for a caller representing one as part
of a longer haplotype record; indels are matched within two bases, since
left-alignment in a repeat can shift them.

Usage:
  python collect_coverage_sweep.py [workspace]
"""
import os, sys, gzip, glob
import pandas as pd, numpy as np

WS = sys.argv[1] if len(sys.argv) > 1 else os.environ["WS"]
SW = f"{WS}/simulation/step14_coverage_sweep"
RES = os.path.expanduser("~/immune_escape_project/results")

ver = pd.read_csv(f"{RES}/verify_step2_per_mutation.tsv", sep="\t")
ver["alt"] = ver.alt.astype(str)
ver = ver[(ver.alt.str.len() == 1) & (ver.alt != "-") & ver.landed]

def calls(path):
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
            for a in f[4].split(","):
                out.append((f[0], int(f[1]), f[3], a))
    return out

levels = sorted(glob.glob(f"{SW}/cov*"),
                key=lambda p: int(os.path.basename(p).replace("cov","").replace("x","")))

rows = []
for lv in levels:
    depth = int(os.path.basename(lv).replace("cov", "").replace("x", ""))
    # the sweep mirrors the cohort layout, so each sample keeps its
    # calls in a mutect2 subdirectory rather than beside the BAM
    for vcf in glob.glob(f"{lv}/*/mutect2/*.filtered.vcf.gz"):
        sid = os.path.basename(os.path.dirname(os.path.dirname(vcf)))
        g = ver[ver["sample"] == sid]
        if not len(g):
            continue
        c = calls(vcf)
        exact = {(x[0], x[1], x[3]) for x in c if len(x[2]) == 1 and len(x[3]) == 1}
        spans = [x for x in c if len(x[2]) == len(x[3]) and len(x[2]) > 1]

        found = 0
        for _, r in g.iterrows():
            key = (r.chrom, int(r.pos), r.alt)
            if key in exact:
                found += 1
                continue
            # a substitution can arrive inside a longer haplotype record
            hit = False
            for ch, pos, ref, alt in spans:
                off = int(r.pos) - pos
                if ch == r.chrom and 0 <= off < len(ref) \
                   and alt[off].upper() == str(r.alt).upper():
                    hit = True
                    break
            found += hit

        rows.append({"depth": depth, "sample": sid,
                     "verified": len(g), "found": found,
                     "recall": 100 * found / len(g),
                     "pass_calls": len(c)})

d = pd.DataFrame(rows)
if d.empty:
    sys.exit("no coverage-sweep VCFs matched")

W = 62
print("=" * W); print(" RECALL AGAINST DEPTH"); print("=" * W)
g = d.groupby("depth").agg(
    samples=("sample", "nunique"),
    verified=("verified", "sum"),
    found=("found", "sum")).reset_index()
g["recall"] = 100 * g.found / g.verified

print(f"\n  {'depth':<9}{'samples':>9}{'verified':>10}{'found':>8}{'recall':>9}")
for _, r in g.iterrows():
    print(f"  {int(r.depth):<9}{int(r.samples):>9}{int(r.verified):>10}"
          f"{int(r.found):>8}{r.recall:>8.1f}%")

base = g[g.depth == g.depth.max()].recall.iloc[0]
print(f"\n  relative to the deepest level:")
for _, r in g.iterrows():
    print(f"    {int(r.depth):>3}x   {r.recall - base:+6.1f} points")

print("\n" + "=" * W); print(" BY ALLELE FRACTION"); print("=" * W)
merged = []
for _, r in d.iterrows():
    sub = ver[ver["sample"] == r["sample"]]
    merged.append(sub.assign(depth=r.depth))
M = pd.concat(merged).drop_duplicates(subset=["sample","chrom","pos","depth"])
BINS=[0,.10,.20,.30,1.01]; LAB=["<10%","10-20%","20-30%",">30%"]
M["bin"]=pd.cut(M.observed_vaf,bins=BINS,labels=LAB,right=False)
print(f"\n  mutations available per bin (same at every depth):")
for l in LAB:
    print(f"    {l:<9}{int((M[M.depth==M.depth.min()].bin==l).sum())}")

out = f"{RES}/coverage_sweep.tsv"
d.to_csv(out, sep="\t", index=False)
g.to_csv(out.replace(".tsv", "_summary.tsv"), sep="\t", index=False)
print(f"\nwritten to {out}")
