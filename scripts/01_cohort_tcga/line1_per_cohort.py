#!/usr/bin/env python3
"""
Line 1 recomputed per cohort, for the slides.

The hypothesis: a mutation in a gene the tumour does not transcribe yields
no neoantigen. Passenger genes can be silenced freely; driver genes cannot,
because the tumour depends on them. So the fraction of mutations sitting in
silenced genes should differ sharply between the two classes.

Silencing is per patient, not per gene in general. A gene counts as silenced
for a given mutation only if that patient's own tumour shows RSEM below the
threshold. A gene silenced in half the cohort and expressed in the other
half is a different claim from a gene silenced in everyone.

Genes are ranked only when they carry enough mutations to estimate a rate:
below the floor a single mutation would read as 0% or 100%.

Usage:
  python line1_per_cohort.py <tcga_root> <out_dir> [rsem_threshold] [min_mutations]
"""
import sys, os
import numpy as np
import pandas as pd

ROOT = sys.argv[1]
OUT = sys.argv[2]
RSEM_OFF = float(sys.argv[3]) if len(sys.argv) > 3 else 5.0
MIN_MUT_BRCA = int(sys.argv[4]) if len(sys.argv) > 4 else 10
MIN_MUT_OV = 5

NONSYN = {"Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
          "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
          "Splice_Site", "Nonstop_Mutation", "Translation_Start_Site"}

os.makedirs(OUT, exist_ok=True)


def norm(b):
    p = str(b).split("-")
    return "-".join(p[:4])[:15] if len(p) >= 4 else str(b)


summary = {}

for cohort, min_mut in [("brca", MIN_MUT_BRCA), ("ov", MIN_MUT_OV)]:
    print(f"\n{'='*72}\n {cohort.upper()}\n{'='*72}", flush=True)

    # ---------------- load ----------------
    m = pd.read_csv(f"{ROOT}/{cohort}/data_mutations.txt", sep="\t",
                    comment="#", low_memory=False)
    e = pd.read_csv(f"{ROOT}/{cohort}/data_mrna_seq_v2_rsem.txt",
                    sep="\t", low_memory=False)

    e = e[e.Hugo_Symbol.notna()].drop_duplicates("Hugo_Symbol")
    e = e.set_index("Hugo_Symbol").drop(columns=["Entrez_Gene_Id"], errors="ignore")
    e.columns = [norm(c) for c in e.columns]
    e = e.loc[:, ~e.columns.duplicated()].clip(lower=0)

    m["sample"] = m.Tumor_Sample_Barcode.map(norm)
    n_mut_all = len(m)
    m = m[m.Variant_Classification.isin(NONSYN)]
    n_nonsyn = len(m)

    # only mutations whose patient has expression, and whose gene is measured
    have_expr = set(e.columns)
    m = m[m["sample"].isin(have_expr)]
    m = m[m.Hugo_Symbol.isin(e.index)]

    n_pat = m["sample"].nunique()
    n_mut = len(m)
    n_genes = m.Hugo_Symbol.nunique()

    print(f"  mutations in file        {n_mut_all:>8,}")
    print(f"  non-synonymous           {n_nonsyn:>8,}")
    print(f"  usable (patient + gene)  {n_mut:>8,}")
    print(f"  patients                 {n_pat:>8}")
    print(f"  distinct mutated genes   {n_genes:>8,}")

    # ---------------- silencing per mutation ----------------
    print("  looking up expression for each mutation ...", flush=True)
    idx = pd.MultiIndex.from_arrays([m.Hugo_Symbol.values, m["sample"].values])
    stacked = e.stack()
    m = m.assign(rsem=stacked.reindex(idx).values)
    m = m[m.rsem.notna()].copy()
    m["silenced"] = m.rsem < RSEM_OFF

    pct = 100 * m.silenced.mean()
    print(f"\n  mutations in silenced genes (RSEM < {RSEM_OFF:g}): "
          f"{int(m.silenced.sum()):,} of {len(m):,}  ({pct:.1f}%)")

    # ---------------- per gene ----------------
    g = m.groupby("Hugo_Symbol").agg(
        n_mutations=("silenced", "size"),
        n_silenced=("silenced", "sum"),
        median_rsem=("rsem", "median"),
    ).reset_index()
    g["pct_silenced"] = (100 * g.n_silenced / g.n_mutations).round(1)
    g = g[g.n_mutations >= min_mut].copy()

    print(f"  genes with >= {min_mut} mutations: {len(g)}")

    kept = g.sort_values(["pct_silenced", "n_mutations"],
                         ascending=[True, False]).head(20)
    hidden = g.sort_values(["pct_silenced", "n_mutations"],
                           ascending=[False, False]).head(20)

    print(f"\n  --- KEPT EXPRESSED (lowest silencing) ---")
    print(f"  {'gene':<12}{'mutations':>10}{'silenced':>10}{'%':>8}{'median RSEM':>14}")
    for _, r in kept.iterrows():
        print(f"  {r.Hugo_Symbol:<12}{r.n_mutations:>10}{r.n_silenced:>10}"
              f"{r.pct_silenced:>8.1f}{r.median_rsem:>14.0f}")

    print(f"\n  --- SILENCED (highest silencing) ---")
    print(f"  {'gene':<12}{'mutations':>10}{'silenced':>10}{'%':>8}{'median RSEM':>14}")
    for _, r in hidden.iterrows():
        print(f"  {r.Hugo_Symbol:<12}{r.n_mutations:>10}{r.n_silenced:>10}"
              f"{r.pct_silenced:>8.1f}{r.median_rsem:>14.0f}")

    g.sort_values("pct_silenced").to_csv(
        f"{OUT}/line1_{cohort}_genes.tsv", sep="\t", index=False)
    summary[cohort] = {
        "patients": n_pat, "mutations": n_mut, "genes": n_genes,
        "pct_silenced": round(pct, 1),
        "n_ranked_genes": len(g),
        "kept": kept, "hidden": hidden, "table": g,
    }

# =====================================================================
print(f"\n{'='*72}\n CONVERGENCE BETWEEN COHORTS\n{'='*72}")

gb = summary["brca"]["table"].set_index("Hugo_Symbol")
go = summary["ov"]["table"].set_index("Hugo_Symbol")
shared = gb.index.intersection(go.index)
print(f"\n  genes ranked in both cohorts: {len(shared)}")

if len(shared) > 5:
    from scipy import stats
    x = gb.loc[shared, "pct_silenced"]
    y = go.loc[shared, "pct_silenced"]
    rho, p = stats.spearmanr(x, y)
    r, pr = stats.pearsonr(x, y)
    print(f"  Spearman rho = {rho:.3f}  (p = {p:.3g})")
    print(f"  Pearson  r   = {r:.3f}  (p = {pr:.3g})")

    both_kept = set(summary["brca"]["kept"].Hugo_Symbol) & \
                set(summary["ov"]["kept"].Hugo_Symbol)
    both_hidden = set(summary["brca"]["hidden"].Hugo_Symbol) & \
                  set(summary["ov"]["hidden"].Hugo_Symbol)
    print(f"\n  in the top-20 kept-expressed list of BOTH cohorts ({len(both_kept)}):")
    print(f"    {', '.join(sorted(both_kept)) or 'none'}")
    print(f"\n  in the top-20 silenced list of BOTH cohorts ({len(both_hidden)}):")
    print(f"    {', '.join(sorted(both_hidden)) or 'none'}")

    pd.DataFrame({"gene": shared,
                  "brca_pct": x.values, "ov_pct": y.values}).to_csv(
        f"{OUT}/line1_shared_genes.tsv", sep="\t", index=False)

print(f"\n{'='*72}\n NUMBERS FOR THE SLIDES\n{'='*72}")
for cohort in ["brca", "ov"]:
    s = summary[cohort]
    print(f"\n  {cohort.upper()}")
    print(f"    patients                {s['patients']:>8}")
    print(f"    mutations analysed      {s['mutations']:>8,}")
    print(f"    in silenced genes       {s['pct_silenced']:>7.1f}%")
    print(f"    genes ranked            {s['n_ranked_genes']:>8}")

print(f"\nwritten to {OUT}/")
