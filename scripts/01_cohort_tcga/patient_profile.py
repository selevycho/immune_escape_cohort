#!/usr/bin/env python3
"""
Assemble a complete profile for the TCGA patient behind a simulated sample.

The simulation reproduces one real case, so every layer of that patient's
data can be pulled from the same cBioPortal tables and placed alongside the
pipeline output. This joins them into one view:

  clinical    - subtype, stage, survival
  mutations   - what the tumour carries, genome-wide and inside the panel
  expression  - RSEM for every mutated gene, and for the antigen machinery
  copy number - GISTIC calls, including the presentation genes
  simulation  - what was injected, what Mutect2 recovered, which peptides
                mhcflurry called binders

The point is to see whether the pipeline's neoantigen predictions survive
contact with the patient's own expression data: a peptide from a gene the
tumour does not transcribe cannot exist.

Usage:
  python patient_profile.py <tcga_dir> <sim_dir> <barcode> <out_dir>
"""
import sys, os
import numpy as np
import pandas as pd

TCGA_DIR = sys.argv[1]
SIM_DIR = sys.argv[2]
BARCODE = sys.argv[3]
OUT_DIR = sys.argv[4]

RSEM_OFF = 5.0
MACHINERY = ["B2M", "HLA-A", "HLA-B", "HLA-C", "TAP1", "TAP2", "TAPBP",
             "NLRC5", "PSMB8", "PSMB9", "CALR", "ERAP1", "ERAP2"]
IMMUNE = ["GZMA", "PRF1", "CD8A", "CD3E", "IFNG", "GZMB", "CD274",
          "PDCD1", "CTLA4", "LAG3", "IDO1"]

NONSYN = {"Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
          "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
          "Splice_Site", "Nonstop_Mutation", "Translation_Start_Site"}

os.makedirs(OUT_DIR, exist_ok=True)


def norm_id(s):
    p = str(s).split("-")
    return "-".join(p[:4])[:15] if len(p) >= 4 else str(s)


SAMPLE = norm_id(BARCODE)
PATIENT = "-".join(SAMPLE.split("-")[:3])

print("=" * 70)
print("PATIENT PROFILE: %s" % SAMPLE)
print("=" * 70)

print("\n[1/6] clinical ...", flush=True)
for fn, key in [("data_clinical_patient.txt", PATIENT),
                ("data_clinical_sample.txt", SAMPLE)]:
    path = "%s/%s" % (TCGA_DIR, fn)
    if not os.path.exists(path):
        continue
    cl = pd.read_csv(path, sep="\t", comment="#", low_memory=False)
    idcol = cl.columns[0]
    row = cl[cl[idcol].astype(str).str.startswith(key)]
    if len(row) == 0:
        continue
    print("  --- %s ---" % fn)
    for c in row.columns:
        v = row.iloc[0][c]
        if pd.notna(v) and str(v).strip() not in ("", "[Not Available]",
                                                  "[Not Applicable]"):
            print("    %-32s %s" % (c, v))

print("\n[2/6] expression ...", flush=True)
expr = pd.read_csv("%s/data_mrna_seq_v2_rsem.txt" % TCGA_DIR,
                   sep="\t", low_memory=False)
expr = expr[expr.Hugo_Symbol.notna()].drop_duplicates("Hugo_Symbol")
expr = expr.set_index("Hugo_Symbol").drop(columns=["Entrez_Gene_Id"],
                                          errors="ignore")
expr.columns = [norm_id(c) for c in expr.columns]
expr = expr.loc[:, ~expr.columns.duplicated()].clip(lower=0)

if SAMPLE not in expr.columns:
    raise SystemExit("sample not found in expression matrix")

e = expr[SAMPLE]
cohort_median = expr.median(axis=1)
print("  genes measured: %d" % len(e))
print("  median RSEM in this tumour: %.1f" % e.median())

print("\n  --- antigen presentation machinery ---")
mach = pd.DataFrame({
    "RSEM": e.reindex(MACHINERY),
    "cohort_median": cohort_median.reindex(MACHINERY),
})
mach["ratio_to_cohort"] = (mach.RSEM / mach.cohort_median).round(2)
mach["percentile"] = [
    round(100 * (expr.loc[g] < e[g]).mean(), 1) if g in expr.index else np.nan
    for g in MACHINERY]
print(mach.round(1).to_string())

print("\n  --- immune activity ---")
imm = pd.DataFrame({
    "RSEM": e.reindex(IMMUNE),
    "cohort_median": cohort_median.reindex(IMMUNE),
})
imm["percentile"] = [
    round(100 * (expr.loc[g] < e[g]).mean(), 1) if g in expr.index else np.nan
    for g in IMMUNE]
print(imm.round(1).to_string())

cyt = np.sqrt(max(e.get("GZMA", 0), 0.01) * max(e.get("PRF1", 0), 0.01))
cyt_all = np.sqrt(expr.loc["GZMA"].clip(lower=0.01) *
                  expr.loc["PRF1"].clip(lower=0.01))
print("\n  cytolytic activity (geometric mean of GZMA and PRF1): %.1f" % cyt)
print("  percentile within the cohort: %.1f" % (100 * (cyt_all < cyt).mean()))

print("\n[3/6] copy number ...", flush=True)
cna_path = "%s/data_cna.txt" % TCGA_DIR
if os.path.exists(cna_path):
    cna = pd.read_csv(cna_path, sep="\t", low_memory=False)
    cna = cna[cna.Hugo_Symbol.notna()].drop_duplicates("Hugo_Symbol")
    cna = cna.set_index("Hugo_Symbol").drop(columns=["Entrez_Gene_Id"],
                                            errors="ignore")
    cna.columns = [norm_id(c) for c in cna.columns]
    cna = cna.loc[:, ~cna.columns.duplicated()]
    if SAMPLE in cna.columns:
        c = cna[SAMPLE]
        print("  genome-wide: %d genes lost, %d gained, %.1f%% altered"
              % (int((c < 0).sum()), int((c > 0).sum()),
                 100 * (c != 0).mean()))
        print("\n  --- presentation machinery ---")
        cm = pd.DataFrame({"GISTIC": c.reindex(MACHINERY)})
        cm["call"] = cm.GISTIC.map({-2: "deep deletion", -1: "shallow loss",
                                    0: "neutral", 1: "gain",
                                    2: "amplification"})
        print(cm.to_string())
    else:
        print("  sample not in the CNA table")
else:
    print("  no CNA file")

print("\n[4/6] mutations ...", flush=True)
maf = pd.read_csv("%s/data_mutations.txt" % TCGA_DIR, sep="\t",
                  low_memory=False)
maf["sample"] = maf.Tumor_Sample_Barcode.map(norm_id)
mine = maf[maf["sample"] == SAMPLE].copy()
mine_ns = mine[mine.Variant_Classification.isin(NONSYN)]
print("  total mutations       : %d" % len(mine))
print("  non-synonymous        : %d" % len(mine_ns))
print("  TMB (non-syn, rough)  : %.1f per Mb of exome" % (len(mine_ns) / 38.0))

den = mine_ns.t_ref_count.fillna(0) + mine_ns.t_alt_count.fillna(0)
mine_ns = mine_ns.assign(VAF=np.where(den > 0, mine_ns.t_alt_count / den,
                                      np.nan))
print("  median VAF            : %.3f" % mine_ns.VAF.median())

in_expr = mine_ns[mine_ns.Hugo_Symbol.isin(e.index)].copy()
in_expr["RSEM"] = in_expr.Hugo_Symbol.map(e)
in_expr["silenced"] = in_expr.RSEM < RSEM_OFF
print("  mutations in measured genes: %d" % len(in_expr))
print("  SILENCED (RSEM < %.0f)      : %d  (%.1f%%)"
      % (RSEM_OFF, int(in_expr.silenced.sum()),
         100 * in_expr.silenced.mean()))

print("\n[5/6] pipeline output ...", flush=True)
truth_path = "%s/truth_set.tsv" % SIM_DIR
neo_path = "%s/neoantigens/neoantigens_per_mutation.tsv" % SIM_DIR
cmp_path = "%s/comparison/truth_vs_calls.tsv" % SIM_DIR

summary = {}
if os.path.exists(truth_path):
    t = pd.read_csv(truth_path, sep="\t")
    summary["injected"] = len(t)
if os.path.exists(cmp_path):
    c = pd.read_csv(cmp_path, sep="\t")
    summary["recovered_by_mutect2"] = int(c.detected_pass.sum())
if os.path.exists(neo_path):
    n = pd.read_csv(neo_path, sep="\t")
    summary["mutations_with_peptides"] = len(n)
    summary["strong_binders"] = int(n.strong.sum())
for k, v in summary.items():
    print("  %-28s %s" % (k, v))

print("\n[6/6] do the predicted neoantigens actually exist? ...", flush=True)
if os.path.exists(neo_path):
    n = pd.read_csv(neo_path, sep="\t")
    n["RSEM_in_this_tumour"] = n.gene.map(e)
    n["expressed"] = n.RSEM_in_this_tumour >= RSEM_OFF
    n = n.sort_values("best_percentile")

    tot_strong = int(n.strong.sum())
    ok_strong = int(n[n.expressed].strong.sum())
    dead_strong = tot_strong - ok_strong

    print("  strong binders predicted        : %d" % tot_strong)
    print("  from genes the tumour expresses : %d" % ok_strong)
    print("  from genes that are OFF here    : %d  (%.1f%% wasted)"
          % (dead_strong,
             100.0 * dead_strong / tot_strong if tot_strong else 0))
    print()
    print("  --- REAL targets: strong binder AND gene expressed ---")
    real = n[(n.expressed) & (n.strong > 0)]
    print(real[["gene", "hgvsp", "VAF", "RSEM_in_this_tumour",
                "strong", "best_percentile"]].head(15).to_string(
        index=False, float_format=lambda x: "%.2f" % x))
    print()
    print("  --- PHANTOM targets: strong binder but gene silent ---")
    dead = n[(~n.expressed) & (n.strong > 0)]
    print(dead[["gene", "hgvsp", "VAF", "RSEM_in_this_tumour",
                "strong", "best_percentile"]].head(15).to_string(
        index=False, float_format=lambda x: "%.2f" % x))

    n.to_csv("%s/%s_neoantigens_with_expression.tsv" % (OUT_DIR, SAMPLE),
             sep="\t", index=False)
else:
    print("  no neoantigen table found")

in_expr.to_csv("%s/%s_mutations_with_expression.tsv" % (OUT_DIR, SAMPLE),
               sep="\t", index=False)
mach.to_csv("%s/%s_machinery.tsv" % (OUT_DIR, SAMPLE), sep="\t")

print()
print("=" * 70)
print("files written to %s" % OUT_DIR)
