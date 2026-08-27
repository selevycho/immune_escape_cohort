#!/usr/bin/env python3
"""
What mhcflurry produced, and what shaped it.

Binder counts are easy to report and easy to misread. A sample with twice
the mutations will show twice the binders, and a sample carrying six
distinct HLA alleles will show more than one carrying four, for reasons
that have nothing to do with its mutations. Both effects are separated
here rather than folded into one number.

Four questions.

  How many peptides were tested, and how many bound?
  Does the count follow the mutations or the genotype?
  Which mutations produced nothing, and why?
  Which alleles do the binders belong to?

The last one matters more than it looks. If most binders sit on one or
two alleles, then losing those alleles removes most of the presented
neoantigens — which is the mechanism the whole project is about.

Usage:
  python check_mhcflurry.py [workspace]
"""
import os
import sys
import glob
import numpy as np
import pandas as pd

WS = sys.argv[1] if len(sys.argv) > 1 else os.environ.get(
    "WS", "/gpfs/bwfor/work/ws/fr_os136-immune_escape")
COHORT = f"{WS}/simulation/cohort"
RES = os.path.expanduser("~/immune_escape_project/results")

man = pd.read_csv(f"{COHORT}/manifest.tsv", sep="\t")

W = 74
per_sample, all_binders, all_peptides = [], [], []

for _, m in man.iterrows():
    sid = m.sample_id
    neo = f"{COHORT}/{sid}/neoantigens"
    per_mut = f"{neo}/neoantigens_per_mutation.tsv"
    if not os.path.exists(per_mut):
        continue

    t = pd.read_csv(per_mut, sep="\t")

    # the genotype, because the number of distinct alleles is half the
    # explanation for how many binders a sample yields
    hla = f"{COHORT}/{sid}/optitype/{sid}_result.tsv"
    alleles = []
    if os.path.exists(hla):
        x = pd.read_csv(hla, sep="\t")
        if len(x):
            r = x.iloc[0]
            alleles = [str(r.get(f"{l}{k}")) for l in "ABC" for k in (1, 2)]
            alleles = [a for a in alleles if a and a.lower() != "nan"]

    row = {"sample": sid, "cohort": m.cohort,
           "mutations": len(t),
           "peptides": int(t.peptides.sum()),
           "strong": int(t.strong.sum()),
           "weak": int(t.weak.sum()),
           "with_binder": int(((t.strong + t.weak) > 0).sum()),
           "n_alleles": len(set(alleles)),
           "best": t.best_percentile.min() if len(t) else np.nan}
    per_sample.append(row)

    b = f"{neo}/neoantigens_binders.tsv"
    if os.path.exists(b):
        bb = pd.read_csv(b, sep="\t")
        bb["sample"] = sid
        all_binders.append(bb)

    a = f"{neo}/neoantigens_all.tsv"
    if os.path.exists(a):
        aa = pd.read_csv(a, sep="\t", usecols=lambda c: c in
                         ("length", "binder", "allele", "mut_peptide"))
        aa["sample"] = sid
        all_peptides.append(aa)

d = pd.DataFrame(per_sample)
if d.empty:
    sys.exit("no mhcflurry results found")

B = pd.concat(all_binders) if all_binders else pd.DataFrame()
P = pd.concat(all_peptides) if all_peptides else pd.DataFrame()

pd.set_option("display.width", 200)

print("=" * W)
print(" 1. PER SAMPLE")
print("=" * W)
print()
cols = ["sample", "cohort", "mutations", "peptides", "strong", "weak",
        "with_binder", "n_alleles"]
print(d[cols].to_string(index=False))

print()
print("=" * W)
print(" 2. TOTALS")
print("=" * W)
print(f"\n  samples predicted        {len(d)}")
print(f"  missense mutations used  {int(d.mutations.sum())}")
print(f"  peptides tested          {int(d.peptides.sum()):,}")
print(f"  strong binders           {int(d.strong.sum())}")
print(f"  weak binders             {int(d.weak.sum())}")
print(f"\n  mutations with at least one binder  "
      f"{int(d.with_binder.sum())} of {int(d.mutations.sum())}"
      f"   ({100*d.with_binder.sum()/d.mutations.sum():.1f}%)")
print(f"  strong binders per mutation         "
      f"{d.strong.sum()/d.mutations.sum():.2f}")

print()
print("=" * W)
print(" 3. MUTATIONS OR GENOTYPE")
print("=" * W)
print(f"\n  A sample yields more binders when it has more mutations, and")
print(f"  also when it carries more distinct alleles. Separating the two:")
print()
if len(d) > 3:
    r_mut = float(np.corrcoef(d.mutations, d.strong)[0, 1])
    r_hla = float(np.corrcoef(d.n_alleles, d.strong)[0, 1])
    print(f"  strong binders against mutations     r = {r_mut:.3f}")
    print(f"  strong binders against allele count  r = {r_hla:.3f}")
    print(f"\n  per mutation, by how many alleles the sample carries:")
    print(f"    {'alleles':<10}{'samples':>9}{'binders/mutation':>20}")
    for n, g in d.groupby("n_alleles"):
        per = (g.strong / g.mutations).median()
        print(f"    {int(n):<10}{len(g):>9}{per:>20.2f}")
    print(f"\n  The count is not a property of the tumour alone. A")
    print(f"  homozygous patient presents fewer neoantigens from the same")
    print(f"  mutations — which is the argument for HLA loss as escape,")
    print(f"  visible here as a side effect of the measurement.")

print()
print("=" * W)
print(" 4. WHICH ALLELES CARRY THE BINDERS")
print("=" * W)
if len(B) and "allele" in B.columns:
    strong = B[B.binder == "STRONG"] if "binder" in B.columns else B
    by_locus = strong.allele.str.extract(r"HLA-([ABC])")[0].value_counts()
    print(f"\n  by locus:")
    for loc, k in by_locus.items():
        print(f"    HLA-{loc}   {k:>5}   ({100*k/len(strong):.0f}%)")

    print(f"\n  most productive alleles:")
    top = strong.allele.value_counts().head(10)
    for a, k in top.items():
        n_car = strong[strong.allele == a]["sample"].nunique()
        print(f"    {a:<16}{k:>5} binders   in {n_car} samples"
              f"   {k/n_car:.1f} each")

    print(f"\n  Concentration matters: if a patient's binders sit mostly on")
    print(f"  one allele, losing that allele removes most of what their")
    print(f"  tumour would have presented.")
    conc = []
    for sid, g in strong.groupby("sample"):
        vc = g.allele.value_counts()
        if len(vc):
            conc.append(100 * vc.iloc[0] / len(g))
    if conc:
        print(f"\n  share carried by each sample's single best allele:")
        print(f"    median {np.median(conc):.0f}%, "
              f"range {min(conc):.0f}-{max(conc):.0f}%")

print()
print("=" * W)
print(" 5. PEPTIDE LENGTH")
print("=" * W)
if len(P) and "length" in P.columns:
    print(f"\n  {'length':<9}{'tested':>10}{'strong':>9}{'rate':>9}")
    for L, g in P.groupby("length"):
        ns = int((g.binder == "STRONG").sum()) if "binder" in g else 0
        print(f"  {int(L):<9}{len(g):>10}{ns:>9}{100*ns/len(g):>8.2f}%")
    print(f"\n  Nine-mers are the canonical class I length and the")
    print(f"  predictor was trained accordingly; the others are included")
    print(f"  because real epitopes are not all nine residues long.")

print()
print("=" * W)
print(" 6. MUTATIONS THAT PRODUCED NOTHING")
print("=" * W)
skipped = []
for _, m in man.iterrows():
    p = f"{COHORT}/{m.sample_id}/neoantigens/skipped_mutations.tsv"
    if os.path.exists(p):
        t = pd.read_csv(p, sep="\t")
        t["sample"] = m.sample_id
        skipped.append(t)
if skipped:
    S = pd.concat(skipped)
    print(f"\n  {len(S)} missense mutations could not be turned into peptides")
    print(f"\n  why:")
    for reason, k in S.reason.str.split("(").str[0].value_counts().head(6).items():
        print(f"    {reason.strip():<48}{k}")
    print(f"\n  These are annotation mismatches, not injection failures:")
    print(f"  the residue TCGA numbered does not match any GENCODE isoform")
    print(f"  of that gene. Reported and skipped, never guessed at.")
else:
    print(f"\n  none")

out = f"{RES}/mhcflurry_summary.tsv"
d.to_csv(out, sep="\t", index=False)
print(f"\nwritten to {out}")
