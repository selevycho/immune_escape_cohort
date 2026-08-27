#!/usr/bin/env python3
"""
Both binding predictors, side by side, in one pass.

mhcflurry and NetMHCpan answer the same question from different starting
points: mhcflurry from the mutations that were placed, NetMHCpan by way of
VEP from the ones Mutect2 recovered. They share no code. Agreement is
evidence that neither carries a systematic error; disagreement should be
locatable in the input rather than in the models.

Everything the two slides need is printed here, in order, with nothing
cut off.

Usage:
  python check_predictors.py [workspace]
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
W = 76


def rule(t):
    print()
    print("=" * W)
    print(f" {t}")
    print("=" * W)


# ===================================================================== mf
mf_rows, mf_binders, mf_pep = [], [], []
for _, m in man.iterrows():
    sid = m.sample_id
    f = f"{COHORT}/{sid}/neoantigens/neoantigens_per_mutation.tsv"
    if not os.path.exists(f):
        continue
    t = pd.read_csv(f, sep="\t")

    alleles = []
    h = f"{COHORT}/{sid}/optitype/{sid}_result.tsv"
    if os.path.exists(h):
        x = pd.read_csv(h, sep="\t")
        if len(x):
            r = x.iloc[0]
            alleles = [str(r.get(f"{l}{k}")) for l in "ABC" for k in (1, 2)]
            alleles = [a for a in alleles if a and a.lower() != "nan"]

    mf_rows.append({"sample": sid, "cohort": m.cohort,
                    "mutations": len(t),
                    "peptides": int(t.peptides.sum()),
                    "strong": int(t.strong.sum()),
                    "weak": int(t.weak.sum()),
                    "with_binder": int(((t.strong + t.weak) > 0).sum()),
                    "n_alleles": len(set(alleles))})

    b = f"{COHORT}/{sid}/neoantigens/neoantigens_binders.tsv"
    if os.path.exists(b):
        bb = pd.read_csv(b, sep="\t",
                         usecols=lambda c: c in ("allele", "binder"))
        bb["sample"] = sid
        mf_binders.append(bb)

    a = f"{COHORT}/{sid}/neoantigens/neoantigens_all.tsv"
    if os.path.exists(a):
        aa = pd.read_csv(a, sep="\t",
                         usecols=lambda c: c in ("length", "binder"))
        mf_pep.append(aa)

MF = pd.DataFrame(mf_rows)
MFB = pd.concat(mf_binders) if mf_binders else pd.DataFrame()
MFP = pd.concat(mf_pep) if mf_pep else pd.DataFrame()

# ===================================================================== pv
pv_rows, pv_frames = [], []
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
    gene_col = next((c for c in ("Gene Name", "Gene") if c in t.columns), None)

    row = {"sample": sid, "cohort": m.cohort, "epitopes": len(t),
           "strong": int((v < 0.5).sum()),
           "weak": int(((v >= 0.5) & (v < 2.0)).sum())}
    if gene_col:
        row["genes"] = t[v < 0.5][gene_col].nunique()
    pv_rows.append(row)

    keep = [c for c in (gene_col, "HLA Allele") if c and c in t.columns]
    if keep:
        f2 = t[keep].copy()
        f2["sample"] = sid
        f2["_p"] = v.values
        pv_frames.append(f2)

PV = pd.DataFrame(pv_rows)
PVB = pd.concat(pv_frames) if pv_frames else pd.DataFrame()

pd.set_option("display.width", 220)

# ===================================================================== out
rule("1. mhcflurry, per sample")
print()
print(MF.to_string(index=False))

rule("2. mhcflurry, totals")
print(f"\n  samples                   {len(MF)}")
print(f"  missense mutations used   {int(MF.mutations.sum())}")
print(f"  peptides tested           {int(MF.peptides.sum()):,}")
print(f"  strong binders            {int(MF.strong.sum())}")
print(f"  weak binders              {int(MF.weak.sum())}")
print(f"  mutations with a binder   {int(MF.with_binder.sum())} of "
      f"{int(MF.mutations.sum())}"
      f"   ({100*MF.with_binder.sum()/MF.mutations.sum():.1f}%)")
print(f"  strong per mutation       "
      f"{MF.strong.sum()/MF.mutations.sum():.2f}")

rule("3. does the yield follow the mutations or the genotype")
if len(MF) > 3:
    print(f"\n  strong binders vs mutations      "
          f"r = {np.corrcoef(MF.mutations, MF.strong)[0,1]:.3f}")
    print(f"  strong binders vs allele count   "
          f"r = {np.corrcoef(MF.n_alleles, MF.strong)[0,1]:.3f}")
    print(f"\n  per mutation, by distinct alleles carried:")
    print(f"    {'alleles':<10}{'samples':>9}{'per mutation':>15}")
    for n, g in MF.groupby("n_alleles"):
        print(f"    {int(n):<10}{len(g):>9}"
              f"{(g.strong/g.mutations).median():>15.2f}")

rule("4. how concentrated each repertoire is")
if len(MFB):
    s = MFB[MFB.binder == "STRONG"]
    conc = []
    for sid, g in s.groupby("sample"):
        vc = g.allele.value_counts()
        if len(vc):
            conc.append(100 * vc.iloc[0] / len(g))
    print(f"\n  share on each sample's single best allele")
    print(f"    median {np.median(conc):.0f}%   "
          f"range {min(conc):.0f}-{max(conc):.0f}%")
    print(f"\n  binders by locus:")
    for loc, k in s.allele.str.extract(r"HLA-([ABC])")[0].value_counts().items():
        print(f"    HLA-{loc}   {k:>5}   ({100*k/len(s):.0f}%)")

rule("5. peptide length")
if len(MFP):
    print(f"\n  {'length':<9}{'tested':>10}{'strong':>9}{'rate':>9}")
    for L, g in MFP.groupby("length"):
        ns = int((g.binder == "STRONG").sum())
        print(f"  {int(L):<9}{len(g):>10}{ns:>9}{100*ns/len(g):>8.2f}%")

rule("6. NetMHCpan, per sample")
print()
print(PV.to_string(index=False))

rule("7. NetMHCpan, totals")
print(f"\n  samples                   {len(PV)}")
print(f"  epitopes scored           {int(PV.epitopes.sum()):,}")
print(f"  strong binders            {int(PV.strong.sum())}")
print(f"  weak binders              {int(PV.weak.sum())}")
if "genes" in PV.columns:
    print(f"  genes with a strong binder, median per sample  "
          f"{PV.genes.median():.0f}")

rule("8. the two routes")
if len(MF) and len(PV):
    j = PV.merge(MF[["sample", "strong", "mutations"]], on="sample",
                 suffixes=("_pv", "_mf"))
    print(f"\n  {'sample':<8}{'mhcflurry':>11}{'NetMHCpan':>11}{'ratio':>9}")
    for _, r in j.iterrows():
        ratio = r.strong_pv / r.strong_mf if r.strong_mf else np.nan
        print(f"  {r['sample']:<8}{int(r.strong_mf):>11}"
              f"{int(r.strong_pv):>11}{ratio:>9.2f}")
    print(f"  {'total':<8}{int(j.strong_mf.sum()):>11}"
          f"{int(j.strong_pv.sum()):>11}"
          f"{j.strong_pv.sum()/j.strong_mf.sum():>9.2f}")
    r = float(np.corrcoef(j.strong_mf, j.strong_pv)[0, 1])
    print(f"\n  correlation across samples   r = {r:.3f}")
    print(f"  median ratio                 "
          f"{(j.strong_pv/j.strong_mf).median():.2f}")

rule("9. why the totals differ")
placed = int(MF.mutations.sum())
print(f"\n  missense placed             {placed}")
print(f"  recall on substitutions     74.5%")
print(f"  expected to reach pVACseq   {placed*0.745:.0f}")
print(f"\n  and the two count different objects: mhcflurry scores one row")
print(f"  per peptide-allele pair, pVACseq one per epitope per transcript")

rule("10. NetMHCpan by allele")
if len(PVB) and "HLA Allele" in PVB.columns:
    s = PVB[PVB._p < 0.5]
    print(f"\n  by locus:")
    for loc, k in s["HLA Allele"].str.extract(r"HLA-([ABC])")[0].value_counts().items():
        print(f"    HLA-{loc}   {k:>5}   ({100*k/len(s):.0f}%)")
    print(f"\n  most productive:")
    for a, k in s["HLA Allele"].value_counts().head(8).items():
        print(f"    {a:<16}{k:>5}")

MF.to_csv(f"{RES}/mhcflurry_summary.tsv", sep="\t", index=False)
PV.to_csv(f"{RES}/pvacseq_summary.tsv", sep="\t", index=False)
print(f"\nwritten to {RES}/mhcflurry_summary.tsv and pvacseq_summary.tsv")
