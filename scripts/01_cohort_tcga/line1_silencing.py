#!/usr/bin/env python3
"""
LINE 1: Transcriptomic Silencing of Mutations (Stealth Mode)

The tumour suppresses genes harbouring somatic mutations to prevent neoantigen
production. A mutation in a gene that is not transcribed yields no peptide and
is therefore invisible to the immune system.

Threshold: RSEM < 5 counts as "not expressed".
Rationale: in TCGA RSEM data, genes below 5 sit at the detection limit and are
reliably OFF, 5-50 is an ambiguous zone, and above 50 genes are reliably ON
(Sci Rep 2025). This is the RSEM equivalent of the <1 FPKM rule.

Outputs (CPTAC-style columns):
  Hugo_Symbol, Total_Mutations, Silenced_Count, Expressed_Count, Silenced_Percent
"""
import sys, os
import numpy as np
import pandas as pd

TCGA_DIR  = sys.argv[1]
OUT_DIR   = sys.argv[2]
COHORT    = sys.argv[3]
THRESHOLD = float(sys.argv[4]) if len(sys.argv) > 4 else 5.0
MIN_MUT   = int(sys.argv[5])   if len(sys.argv) > 5 else 10
TOP_N     = 25

SWEEP = [1.0, 5.0, 10.0, 50.0]

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

print("[2/5] reading mutations ...", flush=True)
maf = pd.read_csv("%s/data_mutations.txt" % TCGA_DIR, sep="\t", low_memory=False,
                  usecols=["Hugo_Symbol", "Tumor_Sample_Barcode", "Variant_Classification"])
maf["sample"] = maf.Tumor_Sample_Barcode.map(norm_id)
maf = maf[maf.Variant_Classification.isin(NONSYN)]
maf = maf[maf["sample"].isin(expr.columns) & maf.Hugo_Symbol.isin(expr.index)]
n_pat = maf["sample"].nunique()
print("      non-synonymous mutations: %d" % len(maf), flush=True)
print("      patients                : %d" % n_pat, flush=True)

print("[3/5] reading expression of each mutated gene in its own tumour ...", flush=True)
E = expr.to_numpy(dtype=float)
gi = {g: i for i, g in enumerate(expr.index)}
si = {s: i for i, s in enumerate(expr.columns)}
maf["expression"] = E[maf.Hugo_Symbol.map(gi).to_numpy(),
                      maf["sample"].map(si).to_numpy()]
maf["silenced"] = maf.expression < THRESHOLD

print("[4/5] threshold sweep ...", flush=True)
sweep = pd.DataFrame([{"RSEM_threshold": t,
                       "silenced_pct_of_all_mutations": round(100 * (maf.expression < t).mean(), 1)}
                      for t in SWEEP])

print("[5/5] building gene tables ...", flush=True)
tab = maf.groupby("Hugo_Symbol").agg(
    Total_Mutations=("silenced", "size"),
    Silenced_Count=("silenced", "sum"),
    Median_RSEM=("expression", "median"),
).reset_index()
tab["Silenced_Count"] = tab.Silenced_Count.astype(int)
tab["Expressed_Count"] = tab.Total_Mutations - tab.Silenced_Count
tab["Silenced_Percent"] = (100 * tab.Silenced_Count / tab.Total_Mutations).round(1)
tab["Median_RSEM"] = tab.Median_RSEM.round(2)
tab["cohort"] = COHORT

tab = tab[tab.Total_Mutations >= MIN_MUT].copy()
COLS = ["Hugo_Symbol", "Total_Mutations", "Silenced_Count",
        "Expressed_Count", "Silenced_Percent"]

low = (tab.sort_values(["Silenced_Percent", "Total_Mutations"], ascending=[True, False])
          .head(TOP_N)[COLS].reset_index(drop=True))
high = (tab.sort_values(["Silenced_Percent", "Total_Mutations"], ascending=[False, False])
           .head(TOP_N)[COLS].reset_index(drop=True))

full = tab.sort_values("Total_Mutations", ascending=False)[COLS + ["Median_RSEM", "cohort"]]
full.to_csv("%s/%s_all_genes.csv" % (OUT_DIR, COHORT), index=False)
low.to_csv("%s/detailed_expressed_low_%s.csv" % (OUT_DIR, COHORT), index=False)
high.to_csv("%s/detailed_silenced_high_%s.csv" % (OUT_DIR, COHORT), index=False)
sweep.to_csv("%s/%s_threshold_sweep.csv" % (OUT_DIR, COHORT), index=False)
maf.to_csv("%s/%s_line1_mutations.tsv" % (OUT_DIR, COHORT), sep="\t", index=False)

print()
print("===== %s - LINE 1 =====" % COHORT.upper())
print("threshold             : RSEM < %s" % THRESHOLD)
print("min mutations per gene: %d" % MIN_MUT)
print("patients              : %d" % n_pat)
print("mutations analysed    : %d" % len(maf))
print("genes passing filter  : %d" % len(tab))
print("SILENCED FRACTION OF ALL MUTATIONS: %.1f%%" % (100 * maf.silenced.mean()))
print()
print("Threshold sensitivity:")
print(sweep.to_string(index=False))
print()
print("### LOW SILENCING (drivers the tumour must keep expressing) ###")
print(low.to_csv(index=False))
print("### HIGH SILENCING (passengers the tumour hides) ###")
print(high.to_csv(index=False))
