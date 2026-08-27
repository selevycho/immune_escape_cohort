#!/usr/bin/env python3
"""
Compare the two neoantigen prediction routes on the same mutations.

Step 5 took mutations from the lifted MAF with HGVSp already annotated by
TCGA, cut peptides out of GENCODE protein sequences, and scored them with
mhcflurry. Step 7 sent the Mutect2 VCF through VEP, let pVACseq build the
peptides with its Wildtype and Frameshift plugins, and scored them with
NetMHCpan-4.1.

The routes share no code and start from different files. Where they agree,
neither is likely to carry a systematic error; where they disagree, the
disagreement localises to a step only one of them performs.

Three comparisons are made:
  counts     strong binders per sample, and whether they correlate
  mutations  which mutations each route found a binder for
  peptides   the actual 8-11mers, and how many are shared

Strong is defined identically for both: percentile rank below 0.5.

Usage:
  python compare_predictors.py [workspace]
"""
import sys, os, glob
import numpy as np
import pandas as pd
from scipy import stats

WS = sys.argv[1] if len(sys.argv) > 1 else os.environ.get(
    "WS", "/gpfs/bwfor/work/ws/fr_os136-immune_escape")

COHORT = f"{WS}/simulation/cohort"
PVAC = f"{WS}/simulation/pvacseq/cohort"
MANIFEST = f"{COHORT}/manifest.tsv"

STRONG = 0.5
WEAK = 2.0

man = pd.read_csv(MANIFEST, sep="\t")

rows = []
gene_overlap = []
peptide_overlap = []

for _, m in man.iterrows():
    sid = m.sample_id
    r = {"sample": sid, "cohort": m.cohort}

    # ---------------- mhcflurry ----------------
    f = f"{COHORT}/{sid}/neoantigens/neoantigens_per_mutation.tsv"
    mf = None
    if os.path.exists(f):
        mf = pd.read_csv(f, sep="\t")
        r["mf_mutations"] = len(mf)
        r["mf_peptides"] = int(mf.peptides.sum()) if "peptides" in mf else None
        r["mf_strong"] = int(mf.strong.sum())
        r["mf_weak"] = int(mf.weak.sum())
        r["mf_mut_with_binder"] = int((mf.strong > 0).sum())

    # ---------------- pVACseq / NetMHCpan ----------------
    files = glob.glob(f"{PVAC}/{sid}/pvacseq_out/MHC_Class_I/*.all_epitopes.tsv")
    pv = None
    if files:
        pv = pd.read_csv(files[0], sep="\t", low_memory=False)
        col = "Best MT Percentile"
        if col in pv.columns:
            v = pd.to_numeric(pv[col], errors="coerce")
            r["pv_epitopes"] = len(pv)
            r["pv_strong"] = int((v < STRONG).sum())
            r["pv_weak"] = int(((v >= STRONG) & (v < WEAK)).sum())
            gcol = "Gene Name" if "Gene Name" in pv.columns else None
            if gcol:
                r["pv_genes_with_strong"] = pv.loc[v < STRONG, gcol].nunique()

    # ---------------- overlap at the gene level ----------------
    if mf is not None and pv is not None and "Gene Name" in pv.columns:
        col = "Best MT Percentile"
        v = pd.to_numeric(pv[col], errors="coerce")
        pv_genes = set(pv.loc[v < STRONG, "Gene Name"].dropna())
        mf_genes = set(mf.loc[mf.strong > 0, "gene"].dropna())
        both = pv_genes & mf_genes
        gene_overlap.append({
            "sample": sid, "cohort": m.cohort,
            "mhcflurry_only": len(mf_genes - pv_genes),
            "shared": len(both),
            "netmhcpan_only": len(pv_genes - mf_genes),
            "jaccard": len(both) / len(pv_genes | mf_genes) if (pv_genes | mf_genes) else None,
        })

    # ---------------- overlap at the peptide level ----------------
    if pv is not None and "MT Epitope Seq" in pv.columns:
        col = "Best MT Percentile"
        v = pd.to_numeric(pv[col], errors="coerce")
        pv_pep = set(pv.loc[v < STRONG, "MT Epitope Seq"].dropna())
        peptide_overlap.append({"sample": sid, "n_pv_peptides": len(pv_pep)})

    rows.append(r)

a = pd.DataFrame(rows)
pd.set_option("display.width", 220)
pd.set_option("display.max_columns", 30)

print("=" * 104)
print(" PER SAMPLE")
print("=" * 104)
cols = [c for c in ["sample", "cohort", "mf_mutations", "mf_strong", "mf_weak",
                    "pv_epitopes", "pv_strong", "pv_weak"] if c in a.columns]
print(a[cols].to_string(index=False))

# ---------------- correlation ----------------
x = a.dropna(subset=["mf_strong", "pv_strong"])
print()
print("=" * 104)
print(" AGREEMENT")
print("=" * 104)
print(f"\n  samples with both routes : {len(x)} of {len(a)}")

if len(x) > 2:
    rp, pp = stats.pearsonr(x.mf_strong, x.pv_strong)
    rs, ps = stats.spearmanr(x.mf_strong, x.pv_strong)
    print(f"\n  strong binders per sample:")
    print(f"    Pearson  r = {rp:.3f}  (p = {pp:.3g})")
    print(f"    Spearman r = {rs:.3f}  (p = {ps:.3g})")
    print(f"\n    mhcflurry : total {int(x.mf_strong.sum()):>6}, "
          f"median {x.mf_strong.median():.0f} per sample")
    print(f"    NetMHCpan : total {int(x.pv_strong.sum()):>6}, "
          f"median {x.pv_strong.median():.0f} per sample")
    ratio = x.pv_strong.sum() / x.mf_strong.sum()
    print(f"\n    NetMHCpan calls {ratio:.1f}x as many strong binders as mhcflurry")

# ---------------- normalised by mutation count ----------------
if "mf_mutations" in x.columns:
    y = x[x.mf_mutations > 0].copy()
    y["mf_per_mut"] = y.mf_strong / y.mf_mutations
    y["pv_per_mut"] = y.pv_strong / y.mf_mutations
    print(f"\n  strong binders per injected mutation:")
    print(f"    mhcflurry : {y.mf_per_mut.median():.2f}  "
          f"(range {y.mf_per_mut.min():.2f} - {y.mf_per_mut.max():.2f})")
    print(f"    NetMHCpan : {y.pv_per_mut.median():.2f}  "
          f"(range {y.pv_per_mut.min():.2f} - {y.pv_per_mut.max():.2f})")

# ---------------- gene overlap ----------------
if gene_overlap:
    g = pd.DataFrame(gene_overlap)
    print()
    print("=" * 104)
    print(" WHICH GENES EACH ROUTE FLAGS")
    print("=" * 104)
    print(f"\n  {'sample':<8}{'mhcflurry only':>16}{'shared':>10}"
          f"{'NetMHCpan only':>16}{'Jaccard':>10}")
    for _, r in g.iterrows():
        j = f"{r.jaccard:.2f}" if pd.notna(r.jaccard) else "  - "
        print(f"  {r['sample']:<8}{r.mhcflurry_only:>16}{r.shared:>10}"
              f"{r.netmhcpan_only:>16}{j:>10}")
    print(f"\n  totals: {g.mhcflurry_only.sum()} mhcflurry-only, "
          f"{g.shared.sum()} shared, {g.netmhcpan_only.sum()} NetMHCpan-only")
    if g.jaccard.notna().any():
        print(f"  median Jaccard index: {g.jaccard.median():.3f}")

# ---------------- by cohort ----------------
print()
print("=" * 104)
print(" BY COHORT")
print("=" * 104)
for coh, gg in x.groupby("cohort"):
    print(f"\n  {coh.upper()}  n={len(gg)}")
    print(f"    mhcflurry strong : {int(gg.mf_strong.sum()):>6}"
          f"   median {gg.mf_strong.median():.0f}")
    print(f"    NetMHCpan strong : {int(gg.pv_strong.sum()):>6}"
          f"   median {gg.pv_strong.median():.0f}")
    if len(gg) > 2:
        rr, pp2 = stats.pearsonr(gg.mf_strong, gg.pv_strong)
        print(f"    correlation      : r = {rr:.3f} (p = {pp2:.3g})")

out = os.path.expanduser("~/immune_escape_project/results/predictor_comparison.tsv")
os.makedirs(os.path.dirname(out), exist_ok=True)
a.to_csv(out, sep="\t", index=False)
if gene_overlap:
    pd.DataFrame(gene_overlap).to_csv(
        out.replace(".tsv", "_genes.tsv"), sep="\t", index=False)
print(f"\nwritten to {out}")
