#!/usr/bin/env python3
"""
Line 3 refined: is antigen-presentation loss above the genome-wide background?
"""
import sys
import pandas as pd, numpy as np
from scipy import stats

TCGA_DIR, OUT_DIR, COHORT = sys.argv[1], sys.argv[2], sys.argv[3]
TARGETS = ["B2M","HLA-A","HLA-B","HLA-C","TAP1","TAP2","TAPBP","NLRC5"]

def norm_id(s):
    p = str(s).split("-"); return "-".join(p[:4])[:15] if len(p)>=4 else str(s)

cna = pd.read_csv(f"{TCGA_DIR}/data_cna.txt", sep="\t", low_memory=False)
cna = cna[cna["Hugo_Symbol"].notna()].drop_duplicates("Hugo_Symbol")
cna = cna.set_index("Hugo_Symbol").drop(columns=["Entrez_Gene_Id"], errors="ignore")
cna.columns = [norm_id(c) for c in cna.columns]
cna = cna.loc[:, ~cna.columns.duplicated()]

# per-sample background: fraction of all genes with any loss
bg = (cna < 0).mean(axis=0)
print(f"[{COHORT}] genome-wide loss fraction: median={bg.median():.3f}  mean={bg.mean():.3f}")
print()

rows = []
for g in TARGETS:
    if g not in cna.index: continue
    lost = cna.loc[g] < 0
    obs  = lost.mean()
    exp  = bg.mean()
    # per-sample: is this gene lost more often than that sample's background?
    enrich = obs / exp if exp > 0 else np.nan
    # binomial-style test against each sample's own background
    ll = np.sum(np.log(np.where(lost, bg, 1-bg)))
    rows.append({"gene": g, "observed_loss": obs, "background": exp,
                 "enrichment": enrich, "n_lost": int(lost.sum()), "n": len(lost)})

res = pd.DataFrame(rows).sort_values("enrichment", ascending=False)
res["cohort"] = COHORT
res.to_csv(f"{OUT_DIR}/{COHORT}_line3_background.tsv", sep="\t", index=False)

print(res.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
print()
print("enrichment ~1.0 -> loss is just background aneuploidy, no selection")
print("enrichment >1.3 -> gene lost more often than expected -> possible selection")
print()

# correlation between background aneuploidy and silencing (confounder check)
l1 = pd.read_csv(f"{OUT_DIR}/{COHORT}_line1_per_patient.tsv", sep="\t").set_index("sample")
common = bg.index.intersection(l1.index)
r, p = stats.spearmanr(bg.loc[common], l1.loc[common, "silencing_ratio_abs"])
print(f"CONFOUNDER CHECK — aneuploidy vs silencing ratio:")
print(f"  Spearman rho = {r:.3f},  p = {p:.4g},  n = {len(common)}")
print("  If rho is clearly non-zero, the B2M-silencing link may be driven")
print("  by overall genome instability rather than immune selection.")
