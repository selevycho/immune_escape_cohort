#!/usr/bin/env python3
"""
Do the three escape mechanisms converge on the same patients, or are they
alternatives?

For every patient with mutations, expression and copy number, three scores
are computed - one per line of defence - and then tested against each other.

  LINE 1  transcriptomic silencing
          fraction of that patient's non-synonymous mutations sitting in
          genes the same tumour does not transcribe (RSEM < 5)

  LINE 2  MHC-I downregulation
          mean z-score of B2M, HLA-A/B/C, TAP1/2, TAPBP and NLRC5 relative
          to the rest of the cohort, sign-flipped so that a high score
          means low expression

  LINE 3  genomic disruption
          number of antigen-presentation genes with a GISTIC loss call,
          reported both raw and relative to that patient's genome-wide
          loss rate, since a tumour that has lost half its genome will hit
          these genes by chance

A fourth axis is added for context, because the single-patient profile
suggested it matters more than any of the three:

  CHECKPOINT  mean z-score of CD274, PDCD1, CTLA4, LAG3, IDO1, HAVCR2

Convergence is tested three ways: correlation between the scores, overlap
of the top quartiles against what independence would predict, and how many
patients are extreme on zero, one, two or three lines.

Usage:
  python convergence_three_lines.py <tcga_dir> <out_dir> <cohort>
"""
import sys, os
import numpy as np
import pandas as pd
from scipy import stats

TCGA_DIR = sys.argv[1]
OUT_DIR = sys.argv[2]
COHORT = sys.argv[3]

RSEM_OFF = 5.0
TOP_FRAC = 0.25          # "extreme" means top quartile on a line

MHC = ["B2M", "HLA-A", "HLA-B", "HLA-C", "TAP1", "TAP2", "TAPBP", "NLRC5"]
CHECKPOINT = ["CD274", "PDCD1", "CTLA4", "LAG3", "IDO1", "HAVCR2"]
CYT = ["GZMA", "PRF1"]
TCELL = ["CD8A", "CD3E", "GZMB", "IFNG"]

NONSYN = {"Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
          "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
          "Splice_Site", "Nonstop_Mutation", "Translation_Start_Site"}

os.makedirs(OUT_DIR, exist_ok=True)


def norm_id(s):
    p = str(s).split("-")
    return "-".join(p[:4])[:15] if len(p) >= 4 else str(s)


def zscore(df, genes):
    """Row-wise z-score across the cohort, averaged over the gene set."""
    present = [g for g in genes if g in df.index]
    if not present:
        return None, []
    sub = np.log2(df.loc[present] + 1)
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1) + 1e-9, axis=0)
    return z.mean(axis=0), present


print("[1/6] reading expression (%s) ..." % COHORT, flush=True)
expr = pd.read_csv("%s/data_mrna_seq_v2_rsem.txt" % TCGA_DIR,
                   sep="\t", low_memory=False)
expr = expr[expr.Hugo_Symbol.notna()].drop_duplicates("Hugo_Symbol")
expr = expr.set_index("Hugo_Symbol").drop(columns=["Entrez_Gene_Id"],
                                          errors="ignore")
expr.columns = [norm_id(c) for c in expr.columns]
expr = expr.loc[:, ~expr.columns.duplicated()].clip(lower=0)
print("      genes=%d  samples=%d" % expr.shape, flush=True)

print("[2/6] LINE 1 - transcriptomic silencing ...", flush=True)
maf = pd.read_csv("%s/data_mutations.txt" % TCGA_DIR, sep="\t",
                  low_memory=False,
                  usecols=["Hugo_Symbol", "Tumor_Sample_Barcode",
                           "Variant_Classification"])
maf["sample"] = maf.Tumor_Sample_Barcode.map(norm_id)
maf = maf[maf.Variant_Classification.isin(NONSYN)]
maf = maf[maf["sample"].isin(expr.columns) & maf.Hugo_Symbol.isin(expr.index)]

E = expr.to_numpy(dtype=float)
gi = {g: i for i, g in enumerate(expr.index)}
si = {s: i for i, s in enumerate(expr.columns)}
maf["rsem"] = E[maf.Hugo_Symbol.map(gi).to_numpy(),
                maf["sample"].map(si).to_numpy()]
maf["silenced"] = maf.rsem < RSEM_OFF

line1 = maf.groupby("sample").agg(
    n_mut=("silenced", "size"),
    n_silenced=("silenced", "sum")).reset_index()
line1["line1_silencing"] = line1.n_silenced / line1.n_mut
line1 = line1[line1.n_mut >= 10].set_index("sample")
print("      patients with >=10 usable mutations: %d" % len(line1), flush=True)
print("      median silencing fraction: %.3f"
      % line1.line1_silencing.median(), flush=True)

print("[3/6] LINE 2 - MHC-I expression ...", flush=True)
mhc_z, mhc_used = zscore(expr, MHC)
print("      genes used: %s" % ", ".join(mhc_used), flush=True)
line2 = pd.DataFrame({"mhc_z": mhc_z})
line2["line2_downreg"] = -line2.mhc_z     # high score = low MHC-I

print("[4/6] LINE 3 - copy number of the machinery ...", flush=True)
cna = pd.read_csv("%s/data_cna.txt" % TCGA_DIR, sep="\t", low_memory=False)
cna = cna[cna.Hugo_Symbol.notna()].drop_duplicates("Hugo_Symbol")
cna = cna.set_index("Hugo_Symbol").drop(columns=["Entrez_Gene_Id"],
                                        errors="ignore")
cna.columns = [norm_id(c) for c in cna.columns]
cna = cna.loc[:, ~cna.columns.duplicated()]

mhc_cna = [g for g in MHC if g in cna.index]
lost = (cna.loc[mhc_cna] < 0).sum(axis=0)
background = (cna < 0).mean(axis=0)          # genome-wide loss rate
expected = background * len(mhc_cna)
line3 = pd.DataFrame({
    "line3_genes_lost": lost,
    "background_loss_rate": background,
    "expected_lost": expected,
})
line3["line3_enrichment"] = line3.line3_genes_lost / (line3.expected_lost + 1e-9)
print("      machinery genes scored: %s" % ", ".join(mhc_cna), flush=True)
print("      median lost: %.1f of %d, median background rate %.3f"
      % (line3.line3_genes_lost.median(), len(mhc_cna),
         line3.background_loss_rate.median()), flush=True)

print("[5/6] context - checkpoints and immune infiltration ...", flush=True)
ckpt_z, ckpt_used = zscore(expr, CHECKPOINT)
cyt_z, cyt_used = zscore(expr, CYT)
tcell_z, tcell_used = zscore(expr, TCELL)
print("      checkpoints: %s" % ", ".join(ckpt_used), flush=True)

ctx = pd.DataFrame({"checkpoint_z": ckpt_z, "cytolytic_z": cyt_z,
                    "tcell_z": tcell_z})

print("[6/6] merging and testing ...", flush=True)
df = line1[["n_mut", "n_silenced", "line1_silencing"]] \
    .join(line2[["line2_downreg"]], how="inner") \
    .join(line3[["line3_genes_lost", "background_loss_rate",
                 "line3_enrichment"]], how="inner") \
    .join(ctx, how="inner")
df["TMB"] = df.n_mut
df["cohort"] = COHORT
print("      patients with all three layers: %d" % len(df), flush=True)

LINES = {"line1_silencing": "Line 1 silencing",
         "line2_downreg": "Line 2 MHC-I down",
         "line3_genes_lost": "Line 3 CN loss"}

for col in LINES:
    thr = df[col].quantile(1 - TOP_FRAC)
    df[col + "_hit"] = df[col] >= thr
df["n_lines_hit"] = df[[c + "_hit" for c in LINES]].sum(axis=1)

df.to_csv("%s/%s_three_lines.tsv" % (OUT_DIR, COHORT), sep="\t")

pd.set_option("display.width", 220)
print()
print("=" * 68)
print("THREE LINES OF DEFENCE - %s  (n = %d)" % (COHORT.upper(), len(df)))
print("=" * 68)

print()
print("--- SCORE DISTRIBUTIONS ---")
desc = df[list(LINES) + ["checkpoint_z", "cytolytic_z", "TMB"]].describe(
    percentiles=[0.25, 0.5, 0.75]).T
print(desc[["mean", "std", "25%", "50%", "75%"]].to_string(
    float_format=lambda x: "%.3f" % x))

print()
print("--- CORRELATION BETWEEN THE LINES ---")
print("  (if the mechanisms converge, these should be positive)")
cols = list(LINES)
rows = []
for i in range(len(cols)):
    for j in range(i + 1, len(cols)):
        a, b = cols[i], cols[j]
        r, p = stats.spearmanr(df[a], df[b])
        rows.append({"pair": "%s vs %s" % (LINES[a], LINES[b]),
                     "spearman_rho": r, "p_value": p})
corr = pd.DataFrame(rows)
print(corr.to_string(index=False, float_format=lambda x: "%.4f" % x))

print()
print("--- OVERLAP OF THE TOP QUARTILES ---")
print("  observed vs expected under independence")
rows = []
n = len(df)
for i in range(len(cols)):
    for j in range(i + 1, len(cols)):
        a, b = cols[i] + "_hit", cols[j] + "_hit"
        obs = int((df[a] & df[b]).sum())
        exp = df[a].sum() * df[b].sum() / n
        tab = pd.crosstab(df[a], df[b])
        try:
            odds, pf = stats.fisher_exact(tab)
        except Exception:
            odds, pf = np.nan, np.nan
        rows.append({"pair": "%s + %s" % (LINES[cols[i]], LINES[cols[j]]),
                     "observed": obs, "expected": exp,
                     "ratio": obs / exp if exp else np.nan,
                     "odds_ratio": odds, "fisher_p": pf})
ov = pd.DataFrame(rows)
print(ov.to_string(index=False, float_format=lambda x: "%.3f" % x))

print()
print("--- HOW MANY LINES PER PATIENT ---")
obs_counts = df.n_lines_hit.value_counts().sort_index()
p_hit = TOP_FRAC
exp_counts = pd.Series(
    [stats.binom.pmf(k, 3, p_hit) * n for k in range(4)], index=range(4))
tbl = pd.DataFrame({"observed": obs_counts, "expected_if_independent":
                    exp_counts}).fillna(0)
tbl["ratio"] = tbl.observed / tbl.expected_if_independent
print(tbl.to_string(float_format=lambda x: "%.1f" % x))
chi2, pchi = stats.chisquare(
    tbl.observed.values,
    f_exp=tbl.expected_if_independent.values * tbl.observed.sum()
    / tbl.expected_if_independent.sum())
print("  chi-square = %.2f   p = %.3g" % (chi2, pchi))

print()
print("--- CONTEXT: what distinguishes patients hitting many lines ---")
ctx_tbl = df.groupby("n_lines_hit").agg(
    patients=("TMB", "size"),
    median_TMB=("TMB", "median"),
    mean_checkpoint_z=("checkpoint_z", "mean"),
    mean_cytolytic_z=("cytolytic_z", "mean"),
    mean_background_loss=("background_loss_rate", "mean"))
print(ctx_tbl.to_string(float_format=lambda x: "%.3f" % x))

print()
print("--- THE CONFOUNDER CHECK ---")
r1, p1 = stats.spearmanr(df.background_loss_rate, df.line3_genes_lost)
r2, p2 = stats.spearmanr(df.background_loss_rate, df.line1_silencing)
r3, p3 = stats.spearmanr(df.TMB, df.line1_silencing)
print("  genome-wide loss rate vs Line 3 : rho = %.3f  p = %.3g" % (r1, p1))
print("  genome-wide loss rate vs Line 1 : rho = %.3f  p = %.3g" % (r2, p2))
print("  mutation count vs Line 1        : rho = %.3f  p = %.3g" % (r3, p3))
print("  (a strong first correlation means Line 3 is mostly aneuploidy)")

print()
print("--- CHECKPOINTS AS A FOURTH ROUTE ---")
thr = df.checkpoint_z.quantile(1 - TOP_FRAC)
df["checkpoint_hit"] = df.checkpoint_z >= thr
none_of_three = df[df.n_lines_hit == 0]
print("  patients hitting none of the three lines: %d" % len(none_of_three))
print("  of those, in the top checkpoint quartile: %d  (%.1f%%)"
      % (int(none_of_three.checkpoint_hit.sum()),
         100 * none_of_three.checkpoint_hit.mean()))
print("  expected by chance: %.1f%%" % (100 * TOP_FRAC))
r, p = stats.spearmanr(df.checkpoint_z, df.cytolytic_z)
print("  checkpoint vs cytolytic activity: rho = %.3f  p = %.3g" % (r, p))

corr.to_csv("%s/%s_line_correlations.tsv" % (OUT_DIR, COHORT),
            sep="\t", index=False)
ov.to_csv("%s/%s_line_overlaps.tsv" % (OUT_DIR, COHORT),
          sep="\t", index=False)
tbl.to_csv("%s/%s_lines_per_patient.tsv" % (OUT_DIR, COHORT), sep="\t")

print()
print("files written to %s" % OUT_DIR)
