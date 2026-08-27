#!/usr/bin/env python3
"""
Are the neoantigen results sound enough to put on a slide?

Steps 5 and 7 were verified before the injection check was found to be
wrong, and the lesson from that was not about mpileup. It was that a
verification script is software and can be measuring something other than
what it claims. So the predictor results get the same question asked of
them before they are used.

Four things are checked.

Do the peptides correspond to mutations that are actually in the files?
Both predictors were run on the truth set, and 31 of those substitutions
are not in the BAMs. Peptides derived from them are predictions about
mutations that do not exist.

Do the two predictors agree on the same objects? A correlation of 0.975
between per-sample binder counts says the totals track each other. It does
not say they are calling the same peptides, and two tools can produce
matching counts from disjoint sets.

Is the count of binders driven by the mutations or by the HLA genotype?
More alleles means more chances for any peptide to bind, so a sample with
six distinct alleles will show more binders than one with four regardless
of its mutations.

And what happened to the mutations that produced nothing? A peptide that
was never generated is different from one generated and predicted not to
bind.

Usage:
  python audit_predictors.py [workspace]
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

v5 = pd.read_csv(f"{RES}/verify_step5.tsv", sep="\t")
v7 = pd.read_csv(f"{RES}/verify_step7.tsv", sep="\t")
ver = pd.read_csv(f"{RES}/verify_step2_per_mutation.tsv", sep="\t")
ver["alt"] = ver.alt.astype(str)
snv = ver[(ver.alt.str.len() == 1) & (ver.alt != "-")]

W = 74
print("=" * W)
print(" 1. WHAT THESE TABLES CONTAIN")
print("=" * W)
print(f"\n  verify_step5.tsv   {len(v5)} rows")
print(f"    {', '.join(list(v5.columns)[:10])}")
print(f"\n  verify_step7.tsv   {len(v7)} rows")
print(f"    {', '.join(list(v7.columns)[:10])}")

for name, t in [("step5", v5), ("step7", v7)]:
    num = t.select_dtypes("number")
    if len(num.columns):
        print(f"\n  {name} totals:")
        for c in num.columns[:8]:
            print(f"    {c:<26}{num[c].sum():>10.0f}"
                  f"   median {num[c].median():>8.2f}")

print()
print("=" * W)
print(" 2. WERE THE PEPTIDES DERIVED FROM MUTATIONS THAT EXIST")
print("=" * W)
print(f"\n  31 of the 1 528 substitutions are not in the BAMs. If the")
print(f"  predictors were run on the truth set rather than on the verified")
print(f"  subset, some peptides describe mutations that were never placed.\n")

absent = snv[~snv.landed]
print(f"  substitutions in the truth set     {len(snv)}")
print(f"  verified present                   {int(snv.landed.sum())}")
print(f"  absent                             {len(absent)}")

pep_files = glob.glob(f"{COHORT}/*/mhcflurry/*.csv") + \
            glob.glob(f"{COHORT}/*/mhcflurry/*predictions*")
print(f"\n  mhcflurry output files found      {len(pep_files)}")

if pep_files:
    ex = pd.read_csv(pep_files[0])
    print(f"  columns in the first one:")
    print(f"    {', '.join(list(ex.columns)[:12])}")
    pos_col = next((c for c in ex.columns
                    if "pos" in c.lower() or "start" in c.lower()), None)
    print(f"  a position column: {pos_col or 'none found'}")
    if pos_col is None:
        print(f"\n  Without a coordinate in the output, a peptide cannot be")
        print(f"  traced back to the mutation that produced it, and whether")
        print(f"  absent mutations contributed cannot be settled from these")
        print(f"  files alone. The upper bound is {len(absent)} of {len(snv)},")
        print(f"  or {100*len(absent)/len(snv):.1f}%.")

print()
print("=" * W)
print(" 3. DO THE TWO PREDICTORS AGREE ON THE SAME PEPTIDES")
print("=" * W)

cmp_path = f"{RES}/predictor_comparison.tsv"
if os.path.exists(cmp_path):
    cmp = pd.read_csv(cmp_path, sep="\t")
    print(f"\n  predictor_comparison.tsv   {len(cmp)} rows")
    print(f"    {', '.join(list(cmp.columns)[:10])}")
    num = cmp.select_dtypes("number")
    if len(num.columns) >= 2:
        a, b = num.columns[0], num.columns[1]
        r = float(np.corrcoef(cmp[a], cmp[b])[0, 1])
        print(f"\n  {a} against {b}:  r = {r:.3f}")
        print(f"    totals {cmp[a].sum():.0f} and {cmp[b].sum():.0f}")
    overlap = [c for c in cmp.columns
               if "overlap" in c.lower() or "shared" in c.lower()
               or "common" in c.lower()]
    if overlap:
        print(f"\n  a column counting shared peptides: {overlap}")
        for c in overlap:
            print(f"    {c}: {cmp[c].sum()}")
    else:
        print(f"\n  No column records how many peptides both tools called.")
        print(f"  A correlation between counts is not agreement on")
        print(f"  identity: two tools can produce matching totals from")
        print(f"  disjoint peptide sets. Worth measuring before the slide")
        print(f"  claims they agree.")

print()
print("=" * W)
print(" 4. IS THE BINDER COUNT DRIVEN BY MUTATIONS OR BY HLA")
print("=" * W)

hla_rows = []
for _, m in pd.read_csv(f"{COHORT}/manifest.tsv", sep="\t").iterrows():
    p = f"{COHORT}/{m.sample_id}/optitype/{m.sample_id}_result.tsv"
    if not os.path.exists(p):
        continue
    t = pd.read_csv(p, sep="\t")
    if not len(t):
        continue
    x = t.iloc[0]
    alleles = {str(x.get(f"{l}{k}")) for l in "ABC" for k in (1, 2)}
    alleles = {a for a in alleles if a and a.lower() != "nan"}
    hla_rows.append({"sample": m.sample_id, "n_alleles": len(alleles)})
hla = pd.DataFrame(hla_rows)

key5 = next((c for c in v5.columns if c.lower() in ("sample", "sample_id")),
            None)
if key5:
    # verify_step5.tsv already carries an allele count, so a plain merge
    # produces n_alleles_x and n_alleles_y and neither name resolves
    left = v5.rename(columns={key5: "sample"})
    left = left.drop(columns=[c for c in left.columns if c == "n_alleles"])
    j = left.merge(hla, on="sample", how="left")
    mut = snv[snv.landed].groupby("sample").size().rename("n_verified")
    j = j.merge(mut, left_on="sample", right_index=True, how="left")
    binder = next((c for c in j.columns
                   if "strong" in c.lower() or "binder" in c.lower()), None)
    if binder and j[binder].notna().any():
        sub = j[j[binder].notna() & j.n_verified.notna() & j.n_alleles.notna()]
        if len(sub) > 3:
            r_mut = float(np.corrcoef(sub.n_verified, sub[binder])[0, 1])
            r_hla = float(np.corrcoef(sub.n_alleles, sub[binder])[0, 1])
            print(f"\n  {binder} against")
            print(f"    verified mutations   r = {r_mut:.3f}")
            print(f"    distinct HLA alleles r = {r_hla:.3f}")
            print(f"\n  alleles per sample: "
                  f"{int(sub.n_alleles.min())}-{int(sub.n_alleles.max())}, "
                  f"median {sub.n_alleles.median():.0f}")
            for n, g in sub.groupby("n_alleles"):
                per = (g[binder] / g.n_verified).median()
                print(f"    {int(n)} alleles: {len(g)} samples, "
                      f"{per:.2f} binders per mutation")
            print(f"\n  If binders per mutation rises with allele count, the")
            print(f"  totals partly measure the genotype rather than the")
            print(f"  mutations, and per-sample comparisons need care.")

print()
print("=" * W)
print(" WHAT TO SETTLE BEFORE THE SLIDE")
print("=" * W)
print(f"""
  - whether peptides trace back to verified mutations or to the whole
    truth set, upper bound {100*len(absent)/len(snv):.1f}% affected
  - whether the two predictors call the same peptides or merely the same
    number of them
  - how much of the binder count is the HLA genotype rather than the
    mutation load
""")
