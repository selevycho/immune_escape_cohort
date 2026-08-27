#!/usr/bin/env python3
"""
LINE 2: MHC-I Transcriptional Downregulation (Software Control)

The tumour attempts to adaptively downregulate the core MHC class I machinery
(HLA-A, HLA-B, HLA-C, B2M) to compensate for mutations it cannot hide at the
transcript level.

Design (same as Line 1 logic):
  1. For every patient, compute the fraction of their mutations that are
     silenced (RSEM < threshold) - this is the Line 1 silencing ratio.
  2. Split patients into HIGH and LOW silencing groups at the median.
  3. Compare MHC-I component expression between the two groups.

If the tumour compensates, HIGH-silencing patients should show LOWER MHC-I.
If expression is flat, "software-level" RNA control is not the escape route.
"""
import sys, os
import numpy as np
import pandas as pd
from scipy import stats

TCGA_DIR = sys.argv[1]
OUT_DIR = sys.argv[2]
COHORT = sys.argv[3]
THRESHOLD = float(sys.argv[4]) if len(sys.argv) > 4 else 5.0

MHC = ["B2M", "HLA-A", "HLA-B", "HLA-C"]
EXTRA = ["TAP1", "TAP2", "TAPBP", "NLRC5", "PSMB8", "PSMB9"]

NONSYN = {"Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
          "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
          "Splice_Site", "Nonstop_Mutation", "Translation_Start_Site"}

os.makedirs(OUT_DIR, exist_ok=True)


def norm_id(s):
    p = str(s).split("-")
    return "-".join(p[:4])[:15] if len(p) >= 4 else str(s)


print("[1/5] reading expression (%s) ..." % COHORT, flush=True)
expr = pd.read_csv("%s/data_mrna_seq_v2_rsem.txt" % TCGA_DIR, sep="\t", low_memory=False)
expr = expr[expr.Hugo_Symbol.notna()].drop_duplicates("Hugo_Symbol")
expr = expr.set_index("Hugo_Symbol").drop(columns=["Entrez_Gene_Id"], errors="ignore")
expr.columns = [norm_id(c) for c in expr.columns]
expr = expr.loc[:, ~expr.columns.duplicated()].clip(lower=0)
print("      genes=%d  samples=%d" % (expr.shape[0], expr.shape[1]), flush=True)

print("[2/5] reading mutations and scoring silencing per patient ...", flush=True)
maf = pd.read_csv("%s/data_mutations.txt" % TCGA_DIR, sep="\t", low_memory=False,
                  usecols=["Hugo_Symbol", "Tumor_Sample_Barcode", "Variant_Classification"])
maf["sample"] = maf.Tumor_Sample_Barcode.map(norm_id)
maf = maf[maf.Variant_Classification.isin(NONSYN)]
maf = maf[maf["sample"].isin(expr.columns) & maf.Hugo_Symbol.isin(expr.index)]

E = expr.to_numpy(dtype=float)
gi = {g: i for i, g in enumerate(expr.index)}
si = {s: i for i, s in enumerate(expr.columns)}
maf["expression"] = E[maf.Hugo_Symbol.map(gi).to_numpy(),
                      maf["sample"].map(si).to_numpy()]
maf["silenced"] = maf.expression < THRESHOLD

pp = maf.groupby("sample").agg(
    n_mut=("silenced", "size"),
    n_sil=("silenced", "sum"),
).reset_index()
pp["silencing_ratio"] = pp.n_sil / pp.n_mut
pp = pp.set_index("sample")
print("      patients with mutations: %d" % len(pp), flush=True)

print("[3/5] splitting into HIGH / LOW silencing groups ...", flush=True)
med = pp.silencing_ratio.median()
pp["group"] = np.where(pp.silencing_ratio > med, "High Silencing", "Low Silencing")
n_hi = int((pp.group == "High Silencing").sum())
n_lo = int((pp.group == "Low Silencing").sum())
print("      median silencing ratio = %.3f" % med, flush=True)
print("      High Silencing n=%d   Low Silencing n=%d" % (n_hi, n_lo), flush=True)

print("[4/5] comparing MHC-I expression ...", flush=True)
genes = [g for g in MHC + EXTRA if g in expr.index]
print("      genes found: %s" % ", ".join(genes), flush=True)

L = np.log2(expr.loc[genes] + 1).T
common = L.index.intersection(pp.index)
L = L.loc[common]
grp = pp.loc[common, "group"]

rows = []
for g in genes:
    hi = L.loc[grp == "High Silencing", g].dropna()
    lo = L.loc[grp == "Low Silencing", g].dropna()
    u, p = stats.mannwhitneyu(hi, lo, alternative="two-sided")
    pooled = np.sqrt((hi.var() + lo.var()) / 2)
    rows.append({
        "Gene": g,
        "n_high": len(hi), "n_low": len(lo),
        "mean_high_log2": hi.mean(), "sem_high": hi.sem(),
        "mean_low_log2": lo.mean(), "sem_low": lo.sem(),
        "difference": hi.mean() - lo.mean(),
        "cohens_d": (hi.mean() - lo.mean()) / pooled if pooled > 0 else np.nan,
        "p_value": p,
    })

res = pd.DataFrame(rows)
res["p_bonferroni"] = (res.p_value * len(res)).clip(upper=1.0)
res["cohort"] = COHORT

print("[5/5] writing ...", flush=True)
out = L.copy()
out["silencing_ratio"] = pp.loc[common, "silencing_ratio"]
out["group"] = grp
out["cohort"] = COHORT
out.index.name = "sample"
out.to_csv("%s/%s_line2_expression.tsv" % (OUT_DIR, COHORT), sep="\t")
res.to_csv("%s/%s_line2_stats.csv" % (OUT_DIR, COHORT), index=False)
pp.to_csv("%s/%s_line2_patient_groups.tsv" % (OUT_DIR, COHORT), sep="\t")

pd.set_option("display.width", 220)
print()
print("===== %s - LINE 2 =====" % COHORT.upper())
print("silencing threshold : RSEM < %s" % THRESHOLD)
print("median split at     : %.3f" % med)
print("High Silencing n=%d   Low Silencing n=%d" % (n_hi, n_lo))
print()
core = res[res.Gene.isin(MHC)]
print("--- CORE MHC-I COMPLEX ---")
print(core[["Gene", "mean_high_log2", "mean_low_log2", "difference",
            "cohens_d", "p_value", "p_bonferroni"]].to_string(
    index=False, float_format=lambda x: "%.4f" % x))
print()
print("--- WIDER ANTIGEN PRESENTATION MACHINERY ---")
print(res[~res.Gene.isin(MHC)][["Gene", "mean_high_log2", "mean_low_log2",
                                "difference", "cohens_d", "p_value",
                                "p_bonferroni"]].to_string(
    index=False, float_format=lambda x: "%.4f" % x))
print()
sig = res[res.p_bonferroni < 0.05]
if len(sig) == 0:
    print("RESULT: no significant difference after correction.")
    print("        MHC-I expression is stable regardless of silencing strategy")
    print("        -> regulatory stalemate; RNA-level control is not the escape route.")
else:
    print("RESULT: significant genes: %s" % ", ".join(sig.Gene))
print()
print("Note on effect size: with n=%d vs %d this design detects" % (n_hi, n_lo))
print("Cohen's d of about %.2f at 80%% power. Smaller effects would be missed."
      % (2.8 * np.sqrt(1.0 / n_hi + 1.0 / n_lo)))
