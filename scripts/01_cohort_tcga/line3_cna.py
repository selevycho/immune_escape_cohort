#!/usr/bin/env python3
"""
LINE 3: Genomic Disruption of Antigen Presentation (Hardware Kill-Switch)

The tumour shifts from regulatory control to physical destruction: copy number
losses remove the antigen-presentation machinery from the genome, producing a
permanent shutdown that no amount of transcriptional regulation can reverse.

B2M sits on chromosome 15 and is the light chain of every MHC-I molecule, so
its physical loss collapses the entire display system regardless of how much
HLA the tumour still transcribes.

This script:
  1. reads GISTIC calls (-2 deep deletion, -1 shallow, 0 neutral, +1/+2 gain)
  2. classifies every patient as B2M deleted or intact
  3. joins that to the Line 1 silencing ratio
  4. tests whether the two escape routes co-occur or act as alternatives
  5. repeats for the other antigen-presentation genes
"""
import sys, os
import numpy as np
import pandas as pd
from scipy import stats

TCGA_DIR = sys.argv[1]
OUT_DIR = sys.argv[2]
COHORT = sys.argv[3]

TARGETS = ["B2M", "HLA-A", "HLA-B", "HLA-C", "TAP1", "TAP2", "TAPBP", "NLRC5"]

os.makedirs(OUT_DIR, exist_ok=True)


def norm_id(s):
    p = str(s).split("-")
    return "-".join(p[:4])[:15] if len(p) >= 4 else str(s)


print("[1/6] reading copy number (%s) ..." % COHORT, flush=True)
cna = pd.read_csv("%s/data_cna.txt" % TCGA_DIR, sep="\t", low_memory=False)
cna = cna[cna.Hugo_Symbol.notna()].drop_duplicates("Hugo_Symbol")
cna = cna.set_index("Hugo_Symbol").drop(columns=["Entrez_Gene_Id"], errors="ignore")
cna.columns = [norm_id(c) for c in cna.columns]
cna = cna.loc[:, ~cna.columns.duplicated()]
print("      genes=%d  samples=%d" % (cna.shape[0], cna.shape[1]), flush=True)

genes = [g for g in TARGETS if g in cna.index]
print("      targets found: %s" % ", ".join(genes), flush=True)

print("[2/6] frequency of loss per gene ...", flush=True)
rows = []
for g in genes:
    v = cna.loc[g].dropna()
    rows.append({
        "Gene": g,
        "n_samples": len(v),
        "deep_deletion_-2": int((v == -2).sum()),
        "shallow_deletion_-1": int((v == -1).sum()),
        "neutral_0": int((v == 0).sum()),
        "gain_1": int((v == 1).sum()),
        "amplification_2": int((v == 2).sum()),
        "pct_any_loss": round(100 * (v < 0).mean(), 1),
        "pct_deep_deletion": round(100 * (v == -2).mean(), 1),
    })
freq = pd.DataFrame(rows)
freq["cohort"] = COHORT

print("[3/6] joining to Line 1 silencing ratio ...", flush=True)
pg = pd.read_csv("%s/%s_line2_patient_groups.tsv" % (OUT_DIR, COHORT),
                 sep="\t").set_index("sample")

common = cna.columns.intersection(pg.index)
sub = cna.loc[genes, common].T
sr = pg.loc[common, "silencing_ratio"]
print("      patients with CNA + silencing ratio: %d" % len(common), flush=True)

print("[4/6] B2M status vs silencing ...", flush=True)
b2m_lost = sub["B2M"] < 0
hi = sr[b2m_lost].dropna()
lo = sr[~b2m_lost].dropna()
if len(hi) > 2 and len(lo) > 2:
    u, p_b2m = stats.mannwhitneyu(hi, lo, alternative="two-sided")
    pooled = np.sqrt((hi.var() + lo.var()) / 2)
    d_b2m = (hi.mean() - lo.mean()) / pooled if pooled > 0 else np.nan
else:
    p_b2m, d_b2m = np.nan, np.nan

print("[5/6] testing every target gene ...", flush=True)
assoc = []
for g in genes:
    lost = sub[g] < 0
    a = sr[lost].dropna()
    b = sr[~lost].dropna()
    if len(a) < 3 or len(b) < 3:
        continue
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    pooled = np.sqrt((a.var() + b.var()) / 2)
    assoc.append({
        "Gene": g,
        "n_lost": len(a), "n_intact": len(b),
        "mean_silencing_lost": a.mean(),
        "mean_silencing_intact": b.mean(),
        "difference": a.mean() - b.mean(),
        "cohens_d": (a.mean() - b.mean()) / pooled if pooled > 0 else np.nan,
        "p_value": p,
    })
assoc = pd.DataFrame(assoc)
assoc["p_bonferroni"] = (assoc.p_value * len(assoc)).clip(upper=1.0)
assoc["cohort"] = COHORT

print("[6/6] writing ...", flush=True)
land = pd.DataFrame({
    "sample": common,
    "silencing_ratio": sr.values,
    "B2M_cna": sub["B2M"].values,
    "B2M_deleted": b2m_lost.values,
})
for g in genes:
    land["%s_cna" % g] = sub[g].values
land["cohort"] = COHORT
land = land.sort_values("silencing_ratio", ascending=False).reset_index(drop=True)

land.to_csv("%s/%s_line3_landscape.tsv" % (OUT_DIR, COHORT), sep="\t", index=False)
freq.to_csv("%s/%s_line3_cna_frequency.csv" % (OUT_DIR, COHORT), index=False)
assoc.to_csv("%s/%s_line3_association.csv" % (OUT_DIR, COHORT), index=False)

pd.set_option("display.width", 220)
print()
print("===== %s - LINE 3 =====" % COHORT.upper())
print("patients analysed: %d" % len(land))
print()
print("--- COPY NUMBER STATUS OF ANTIGEN PRESENTATION GENES ---")
print(freq[["Gene", "n_samples", "deep_deletion_-2", "shallow_deletion_-1",
            "neutral_0", "gain_1", "amplification_2",
            "pct_any_loss", "pct_deep_deletion"]].to_string(index=False))
print()
n_del = int(b2m_lost.sum())
print("--- B2M, THE KILL-SWITCH ---")
print("B2M lost (any level): %d / %d  (%.1f%%)"
      % (n_del, len(land), 100.0 * n_del / len(land)))
print("  mean silencing ratio, B2M deleted : %.3f  (n=%d)" % (hi.mean(), len(hi)))
print("  mean silencing ratio, B2M intact  : %.3f  (n=%d)" % (lo.mean(), len(lo)))
print("  difference = %.3f   Cohen's d = %.3f   p = %.4g"
      % (hi.mean() - lo.mean(), d_b2m, p_b2m))
print()
print("--- SILENCING RATIO BY LOSS STATUS, ALL TARGETS ---")
print(assoc[["Gene", "n_lost", "n_intact", "mean_silencing_lost",
             "mean_silencing_intact", "difference", "cohens_d",
             "p_value", "p_bonferroni"]].to_string(
    index=False, float_format=lambda x: "%.4f" % x))
print()
print("How to read this:")
print("  large positive Cohen's d -> the two escape routes COMBINE:")
print("     tumours that delete the machinery also hide their mutations")
print("  d near zero -> the routes are INDEPENDENT or ALTERNATIVE:")
print("     deleting the hardware makes hiding at the RNA level unnecessary")
