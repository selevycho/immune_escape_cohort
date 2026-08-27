#!/usr/bin/env python3
"""
What pVACseq and NetMHCpan produced, and how it compares.

This is the second route to an answer step 5 already gave. The two share
no code and start from different files: mhcflurry from the mutations that
were placed, pVACseq from the ones Mutect2 recovered, by way of VEP. Where
they agree, neither is carrying a systematic error. Where they differ, the
reason should be locatable.

Three questions.

  How many epitopes, and how many bind?
  Do the two routes agree on which mutations produce binders?
  Where they disagree, is it the input or the model?

The last one has an expected answer. pVACseq sees only what the caller
recovered, which at 74.5% recall is a quarter fewer mutations to work
from. If the difference is entirely that, the two predictors agree; if it
is more, they do not.

Usage:
  python check_pvacseq.py [workspace]
"""
import os
import sys
import glob
import numpy as np
import pandas as pd

WS = sys.argv[1] if len(sys.argv) > 1 else os.environ.get(
    "WS", "/gpfs/bwfor/work/ws/fr_os136-immune_escape")
COHORT = f"{WS}/simulation/cohort"
PVAC = f"{WS}/simulation/pvacseq/cohort"
RES = os.path.expanduser("~/immune_escape_project/results")

man = pd.read_csv(f"{COHORT}/manifest.tsv", sep="\t")

W = 74
rows, frames = [], []

for _, m in man.iterrows():
    sid = m.sample_id
    hits = glob.glob(f"{PVAC}/{sid}/pvacseq_out/MHC_Class_I/*.all_epitopes.tsv")
    if not hits:
        continue
    t = pd.read_csv(hits[0], sep="\t", low_memory=False)

    col = "Best MT Percentile"
    if col not in t.columns:
        col = next((c for c in t.columns if "Percentile" in c), None)
    if col is None:
        continue
    v = pd.to_numeric(t[col], errors="coerce")

    # the mutation a peptide came from, so the two routes can be compared
    # at the level of mutations rather than of peptide counts
    key = next((c for c in ("Mutation", "Protein Position", "Transcript")
                if c in t.columns), None)
    gene_col = next((c for c in ("Gene Name", "Gene") if c in t.columns), None)

    strong = t[v < 0.5]
    row = {"sample": sid, "cohort": m.cohort,
           "epitopes": len(t),
           "strong": int((v < 0.5).sum()),
           "weak": int(((v >= 0.5) & (v < 2.0)).sum())}
    if gene_col:
        row["genes_with_strong"] = strong[gene_col].nunique()
    rows.append(row)

    keep = [c for c in (gene_col, col, "HLA Allele", "MT Epitope Seq")
            if c and c in t.columns]
    if keep:
        f = t[keep].copy()
        f["sample"] = sid
        f["_p"] = v.values
        frames.append(f)

d = pd.DataFrame(rows)
if d.empty:
    sys.exit("no pVACseq results found")
P = pd.concat(frames) if frames else pd.DataFrame()

mf = pd.read_csv(f"{RES}/mhcflurry_summary.tsv", sep="\t") \
    if os.path.exists(f"{RES}/mhcflurry_summary.tsv") else None

pd.set_option("display.width", 200)

print("=" * W)
print(" 1. PER SAMPLE")
print("=" * W)
print()
print(d.to_string(index=False))

print()
print("=" * W)
print(" 2. TOTALS")
print("=" * W)
print(f"\n  samples                 {len(d)}")
print(f"  epitopes scored         {int(d.epitopes.sum()):,}")
print(f"  strong binders          {int(d.strong.sum())}")
print(f"  weak binders            {int(d.weak.sum())}")

print()
print("=" * W)
print(" 3. THE TWO ROUTES")
print("=" * W)
if mf is not None:
    j = d.merge(mf[["sample", "strong", "mutations"]], on="sample",
                suffixes=("_pv", "_mf"))
    print(f"\n  {'sample':<8}{'mhcflurry':>11}{'NetMHCpan':>11}{'ratio':>9}")
    for _, r in j.iterrows():
        ratio = r.strong_pv / r.strong_mf if r.strong_mf else np.nan
        print(f"  {r['sample']:<8}{int(r.strong_mf):>11}"
              f"{int(r.strong_pv):>11}{ratio:>9.2f}")
    print(f"  {'total':<8}{int(j.strong_mf.sum()):>11}"
          f"{int(j.strong_pv.sum()):>11}"
          f"{j.strong_pv.sum()/j.strong_mf.sum():>9.2f}")

    if len(j) > 3:
        r = float(np.corrcoef(j.strong_mf, j.strong_pv)[0, 1])
        print(f"\n  correlation across samples   r = {r:.3f}")
        print(f"\n  A correlation this high with a ratio away from one means")
        print(f"  the two rank samples the same way while counting")
        print(f"  different totals — which is what different inputs to the")
        print(f"  same question look like.")

print()
print("=" * W)
print(" 4. WHY THE TOTALS DIFFER")
print("=" * W)
print(f"\n  mhcflurry starts from the mutations that were placed;")
print(f"  pVACseq from the ones Mutect2 recovered.")
if mf is not None:
    placed = int(mf.mutations.sum())
    print(f"\n    missense placed            {placed}")
    print(f"    recall on substitutions    74.5%")
    print(f"    expected to reach pVACseq  {placed * 0.745:.0f}")
    print(f"\n  pVACseq also counts differently: one row per epitope per")
    print(f"  transcript, where mhcflurry scores one row per")
    print(f"  peptide-allele pair. The counts are not the same object.")

print()
print("=" * W)
print(" 5. WHICH ALLELES")
print("=" * W)
if len(P) and "HLA Allele" in P.columns:
    s = P[P._p < 0.5]
    by_locus = s["HLA Allele"].str.extract(r"HLA-([ABC])")[0].value_counts()
    print(f"\n  by locus:")
    for loc, k in by_locus.items():
        print(f"    HLA-{loc}   {k:>5}   ({100*k/len(s):.0f}%)")
    print(f"\n  most productive:")
    for a, k in s["HLA Allele"].value_counts().head(8).items():
        print(f"    {a:<16}{k:>5}")

out = f"{RES}/pvacseq_summary.tsv"
d.to_csv(out, sep="\t", index=False)
print(f"\nwritten to {out}")
