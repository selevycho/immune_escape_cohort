#!/usr/bin/env python3
"""
LINE 1 - convergence between breast and ovarian cancer.

If the stealth strategy is a general mechanism rather than a tissue quirk,
the same genes should be silenced in both cancers, and the same genes should
be forced to stay expressed in both.

Three views:
  1. Overlap of the high-silencing and low-silencing gene lists
  2. Correlation of Silenced_Percent across all genes tested in both cohorts
  3. The genes that diverge most between the two cancers
"""
import sys, os
import numpy as np
import pandas as pd
from scipy import stats

RES_DIR = sys.argv[1]
TOP_N = int(sys.argv[2]) if len(sys.argv) > 2 else 25

brca = pd.read_csv("%s/brca_all_genes.csv" % RES_DIR)
ov = pd.read_csv("%s/ov_all_genes.csv" % RES_DIR)
print("[1/4] genes passing filter: BRCA %d, OV %d" % (len(brca), len(ov)), flush=True)

m = brca.merge(ov, on="Hugo_Symbol", suffixes=("_brca", "_ov"))
print("      tested in BOTH cohorts   : %d" % len(m), flush=True)

print("[2/4] correlating silencing between cancers ...", flush=True)
r_s, p_s = stats.spearmanr(m.Silenced_Percent_brca, m.Silenced_Percent_ov)
r_p, p_p = stats.pearsonr(m.Silenced_Percent_brca, m.Silenced_Percent_ov)

print("[3/4] overlap of the extreme lists ...", flush=True)


def top_low(df, n):
    return set(df.sort_values(["Silenced_Percent", "Total_Mutations"],
                              ascending=[True, False]).head(n).Hugo_Symbol)


def top_high(df, n):
    return set(df.sort_values(["Silenced_Percent", "Total_Mutations"],
                              ascending=[False, False]).head(n).Hugo_Symbol)


universe = set(m.Hugo_Symbol)
bl = top_low(brca, TOP_N) & universe
ol = top_low(ov, TOP_N) & universe
bh = top_high(brca, TOP_N) & universe
oh = top_high(ov, TOP_N) & universe


def hyper(a, b, N):
    k = len(a & b)
    p = stats.hypergeom.sf(k - 1, N, len(a), len(b)) if k > 0 else 1.0
    exp = len(a) * len(b) / N if N > 0 else np.nan
    return k, exp, p


N = len(universe)
k_low, e_low, p_low = hyper(bl, ol, N)
k_high, e_high, p_high = hyper(bh, oh, N)

print("[4/4] writing ...", flush=True)
shared_low = sorted(bl & ol)
shared_high = sorted(bh & oh)

m["delta_silencing"] = m.Silenced_Percent_brca - m.Silenced_Percent_ov
out = m[["Hugo_Symbol",
         "Total_Mutations_brca", "Silenced_Percent_brca",
         "Total_Mutations_ov", "Silenced_Percent_ov", "delta_silencing"]]
out.sort_values("delta_silencing").to_csv(
    "%s/convergence_all_genes.csv" % RES_DIR, index=False)
pd.DataFrame({"Hugo_Symbol": shared_low}).to_csv(
    "%s/convergence_shared_LOW_silencing.csv" % RES_DIR, index=False)
pd.DataFrame({"Hugo_Symbol": shared_high}).to_csv(
    "%s/convergence_shared_HIGH_silencing.csv" % RES_DIR, index=False)

pd.set_option("display.width", 200)
print()
print("========== LINE 1 CONVERGENCE : BRCA vs OV ==========")
print("genes tested in both cohorts: %d" % N)
print()
print("--- Correlation of Silenced_Percent across shared genes ---")
print("  Spearman rho = %.3f   p = %.3g" % (r_s, p_s))
print("  Pearson  r   = %.3f   p = %.3g" % (r_p, p_p))
print()
print("--- Overlap of the top-%d LOW-silencing lists (forced-expressed drivers) ---" % TOP_N)
print("  shared: %d   expected by chance: %.1f   hypergeometric p = %.3g" % (k_low, e_low, p_low))
print("  genes : %s" % (", ".join(shared_low) if shared_low else "(none)"))
print()
print("--- Overlap of the top-%d HIGH-silencing lists (hidden passengers) ---" % TOP_N)
print("  shared: %d   expected by chance: %.1f   hypergeometric p = %.3g" % (k_high, e_high, p_high))
print("  genes : %s" % (", ".join(shared_high) if shared_high else "(none)"))
print()
print("--- 15 genes MORE silenced in OVARIAN than in BREAST ---")
print(out.sort_values("delta_silencing").head(15).to_string(index=False))
print()
print("--- 15 genes MORE silenced in BREAST than in OVARIAN ---")
print(out.sort_values("delta_silencing", ascending=False).head(15).to_string(index=False))
