#!/usr/bin/env python3
"""
Which five samples make the best demo.

A demo that runs every stage and produces nothing at one of them is worse
than one that skips the stage honestly, so the samples have to be chosen
against every stage rather than against mutation count alone. Six
requirements, and they pull in different directions.

  Enough substitutions to make recall meaningful, but not so many that
  injection takes an hour. Twenty to sixty.
  At least two indels, or step 2 has nothing to place and step 3 nothing
  to find.
  At least two heterozygous HLA loci, or LOHHLA has no allelic imbalance
  to test and returns nothing.
  All three HLA loci typed, since a missing HLA-C makes the OptiType
  output look broken to someone seeing it for the first time.
  Some peptides predicted to bind, or steps 5 and 7 produce empty tables.
  Coverage near the cohort median, so the numbers a viewer sees match the
  numbers on the slides.

Samples are scored on each and the ranking is printed with the components
visible, because a single score hides which requirement a sample fails.

Usage:
  python pick_demo_samples.py [workspace] [--n=5]
"""
import os
import sys
import glob
import numpy as np
import pandas as pd

args = [a for a in sys.argv[1:] if not a.startswith("--")]
WS = args[0] if args else os.environ.get(
    "WS", "/gpfs/bwfor/work/ws/fr_os136-immune_escape")
N_PICK = 5
for a in sys.argv[1:]:
    if a.startswith("--n="):
        N_PICK = int(a.split("=", 1)[1])

COHORT = f"{WS}/simulation/cohort"
RES = os.path.expanduser("~/immune_escape_project/results")

man = pd.read_csv(f"{COHORT}/manifest.tsv", sep="\t")
ind = pd.read_csv(f"{WS}/simulation/indels/indel_truth_all.tsv", sep="\t")
n_indel = ind.groupby("sample").size()

ver = pd.read_csv(f"{RES}/verify_step2_per_mutation.tsv", sep="\t")
ver["alt"] = ver.alt.astype(str)
ver = ver[(ver.alt.str.len() == 1) & (ver.alt != "-")]

v1 = pd.read_csv(f"{RES}/verify_step1.tsv", sep="\t") \
    if os.path.exists(f"{RES}/verify_step1.tsv") else None
v5 = pd.read_csv(f"{RES}/verify_step5.tsv", sep="\t") \
    if os.path.exists(f"{RES}/verify_step5.tsv") else None
v7 = pd.read_csv(f"{RES}/verify_step7.tsv", sep="\t") \
    if os.path.exists(f"{RES}/verify_step7.tsv") else None
irec = pd.read_csv(f"{RES}/indel_recall.tsv", sep="\t") \
    if os.path.exists(f"{RES}/indel_recall.tsv") else None
s3 = pd.read_csv(f"{RES}/step3_rescored_full.tsv", sep="\t") \
    if os.path.exists(f"{RES}/step3_rescored_full.tsv") else None

rows = []
for _, m in man.iterrows():
    s = m.sample_id
    r = {"sample": s, "cohort": m.cohort, "pop": m.get("superpopulation")}

    g = ver[ver["sample"] == s]
    r["snv"] = len(g)
    r["snv_landed"] = int(g.landed.sum())
    r["indels"] = int(n_indel.get(s, 0))

    if s3 is not None:
        h = s3[s3["sample"] == s]
        r["snv_recall"] = round(100 * h["pass"].mean(), 1) if len(h) else None
    if irec is not None:
        h = irec[irec["sample"] == s]
        r["indel_found"] = int(h.found.sum()) if len(h) else 0

    # HLA: how many loci were typed and how many of those are heterozygous,
    # since LOHHLA needs a locus where the two alleles differ
    p = f"{COHORT}/{s}/optitype/{s}_result.tsv"
    typed = het = 0
    if os.path.exists(p):
        t = pd.read_csv(p, sep="\t")
        if len(t):
            x = t.iloc[0]
            for l in "ABC":
                a1, a2 = str(x.get(f"{l}1")), str(x.get(f"{l}2"))
                if a1 not in ("nan", "None", ""):
                    typed += 1
                    if a1 != a2:
                        het += 1
    r["hla_typed"] = typed
    r["hla_het"] = het

    if v5 is not None:
        h = v5[v5["sample"] == s]
        if len(h):
            r["mf_strong"] = int(h.iloc[0].get("strong", 0))
            r["mf_peptides"] = int(h.iloc[0].get("peptides", 0))
    if v7 is not None:
        h = v7[v7["sample"] == s]
        if len(h):
            r["pv_strong"] = int(h.iloc[0].get("strong", 0))

    if v1 is not None:
        h = v1[v1["sample"] == s]
        if len(h):
            r["depth"] = round(float(h.iloc[0].panel_mean_depth), 1)
            r["size_MB"] = round(float(h.iloc[0].size_MB), 0)

    rows.append(r)

d = pd.DataFrame(rows)
med_depth = d.depth.median() if "depth" in d.columns else 34.5

# ------------------------------------------------------------------ score
def band(v, lo, hi, best_lo, best_hi):
    """1.0 inside the preferred band, tapering to 0 at the hard limits."""
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return 0.0
    if best_lo <= v <= best_hi:
        return 1.0
    if v < lo or v > hi:
        return 0.0
    if v < best_lo:
        return (v - lo) / (best_lo - lo)
    return (hi - v) / (hi - best_hi)


crit = {}
crit["substitutions"] = d.snv.apply(lambda v: band(v, 10, 120, 20, 60))
crit["indels"] = d.indels.apply(lambda v: min(v / 3.0, 1.0))
crit["HLA typed"] = d.hla_typed / 3.0
crit["HLA heterozygous"] = d.hla_het / 3.0
if "mf_strong" in d.columns:
    crit["binders"] = d.mf_strong.fillna(0).apply(lambda v: min(v / 20.0, 1.0))
if "depth" in d.columns:
    crit["typical depth"] = d.depth.apply(
        lambda v: max(0.0, 1.0 - abs(v - med_depth) / 12.0))
if "snv_recall" in d.columns:
    # a sample whose recall matches the cohort is representative; one far
    # from it would make the demo look better or worse than the work is
    crit["typical recall"] = d.snv_recall.apply(
        lambda v: max(0.0, 1.0 - abs((v or 0) - 74.5) / 30.0))

C = pd.DataFrame(crit)
d["score"] = C.mean(axis=1).round(3)
for k in C.columns:
    d[f"s_{k}"] = C[k].round(2)

pd.set_option("display.width", 250)
W = 84
print("=" * W)
print(" ALL SAMPLES, RANKED")
print("=" * W)
show = ["sample", "cohort", "snv", "indels", "hla_typed", "hla_het"]
for c in ("mf_strong", "pv_strong", "depth", "snv_recall", "indel_found"):
    if c in d.columns:
        show.append(c)
show.append("score")
print()
print(d.sort_values("score", ascending=False)[show].to_string(index=False))

print()
print("=" * W)
print(" HARD REQUIREMENTS")
print("=" * W)
ok = d[(d.indels >= 2) & (d.hla_typed == 3) & (d.hla_het >= 2) &
       (d.snv.between(15, 80))]
print(f"\n  2+ indels, all three HLA loci typed, 2+ heterozygous,")
print(f"  15-80 substitutions\n")
print(f"  {len(ok)} of {len(d)} samples qualify")
if "mf_strong" in d.columns:
    ok = ok[ok.mf_strong.fillna(0) > 0]
    print(f"  {len(ok)} of those also produce binders")

print()
print("=" * W)
print(f" THE {N_PICK} PICKED")
print("=" * W)
pick = ok.sort_values("score", ascending=False).head(N_PICK)
if len(pick) < N_PICK:
    print(f"\n  only {len(pick)} meet every requirement; filling from the")
    print(f"  ranking regardless\n")
    extra = d[~d["sample"].isin(pick["sample"])].sort_values(
        "score", ascending=False).head(N_PICK - len(pick))
    pick = pd.concat([pick, extra])

print()
print(pick[show].to_string(index=False))

print(f"\n  component scores:\n")
comp = ["sample"] + [f"s_{k}" for k in C.columns]
print(pick[comp].to_string(index=False))

print(f"\n  totals for the demo:")
print(f"    substitutions      {int(pick.snv.sum())}")
print(f"    indels             {int(pick.indels.sum())}")
if "mf_strong" in pick.columns:
    print(f"    binders            {int(pick.mf_strong.fillna(0).sum())}")
if "size_MB" in pick.columns:
    print(f"    BAM data           "
          f"{2 * pick.size_MB.sum() / 1024:.1f} GB  (normal + tumour)")
print(f"    cohorts            "
      f"{dict(pick.cohort.value_counts())}")
if "pop" in pick.columns:
    print(f"    populations        {sorted(set(pick['pop'].dropna()))}")

print(f"\n  SAMPLES=\"{' '.join(pick['sample'])}\"")

out = f"{RES}/demo_sample_choice.tsv"
d.sort_values("score", ascending=False).to_csv(out, sep="\t", index=False)
print(f"\nwritten to {out}")
