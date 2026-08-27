#!/usr/bin/env python3
"""
Count what is actually available in each TCGA study, layer by layer.

The number of patients in a study is not the number usable for analysis.
A patient contributes to Line 1 only if they have both mutations and
expression; to Line 3 only if they also have copy number. Reporting the
study size when the analysis used a smaller intersection would overstate
the work.

Sample barcodes are truncated to the first four fields (TCGA-XX-YYYY-NN)
before comparison. cBioPortal appends different suffixes in different
files - -01 for primary tumour in one, longer aliquot identifiers in
another - and untruncated barcodes do not match across files.

Usage:
  python count_cohorts.py <tcga_root>
"""
import sys, os
import pandas as pd

ROOT = sys.argv[1] if len(sys.argv) > 1 else os.environ.get(
    "WS", "/gpfs/bwfor/work/ws/fr_os136-immune_escape") + "/tcga"

FILES = {
    "mutations": "data_mutations.txt",
    "expression": "data_mrna_seq_v2_rsem.txt",
    "cna": "data_cna.txt",
    "clinical_patient": "data_clinical_patient.txt",
    "clinical_sample": "data_clinical_sample.txt",
}


def norm(barcode):
    """TCGA-AC-A23H-01A-11D-A159-09 -> TCGA-AC-A23H-01"""
    p = str(barcode).split("-")
    return "-".join(p[:4])[:15] if len(p) >= 4 else str(barcode)


def patient(barcode):
    """TCGA-AC-A23H-01 -> TCGA-AC-A23H"""
    return "-".join(str(barcode).split("-")[:3])


results = {}

for cohort in ["brca", "ov"]:
    d = f"{ROOT}/{cohort}"
    print(f"reading {cohort.upper()} ...", flush=True)
    r = {}

    # ---- mutations ----
    p = f"{d}/{FILES['mutations']}"
    if os.path.exists(p):
        m = pd.read_csv(p, sep="\t", comment="#", low_memory=False)
        ids = {norm(x) for x in m.Tumor_Sample_Barcode.dropna()}
        r["mutations_samples"] = len(ids)
        r["mutations_rows"] = len(m)
        r["mut_ids"] = ids
        nonsyn = {"Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
                  "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
                  "Splice_Site", "Nonstop_Mutation", "Translation_Start_Site"}
        r["mutations_nonsyn"] = int(m.Variant_Classification.isin(nonsyn).sum())

    # ---- expression ----
    p = f"{d}/{FILES['expression']}"
    if os.path.exists(p):
        e = pd.read_csv(p, sep="\t", nrows=1)
        cols = [c for c in e.columns if c.startswith("TCGA")]
        ids = {norm(c) for c in cols}
        r["expression_samples"] = len(ids)
        r["expr_ids"] = ids
        full = pd.read_csv(p, sep="\t", usecols=[0], low_memory=False)
        r["expression_genes"] = full.iloc[:, 0].notna().sum()

    # ---- copy number ----
    p = f"{d}/{FILES['cna']}"
    if os.path.exists(p):
        c = pd.read_csv(p, sep="\t", nrows=1)
        cols = [x for x in c.columns if x.startswith("TCGA")]
        ids = {norm(x) for x in cols}
        r["cna_samples"] = len(ids)
        r["cna_ids"] = ids

    # ---- clinical ----
    p = f"{d}/{FILES['clinical_patient']}"
    if os.path.exists(p):
        cl = pd.read_csv(p, sep="\t", comment="#", low_memory=False)
        r["clinical_patients"] = len(cl)
        idc = cl.columns[0]
        r["clin_ids"] = set(cl[idc].dropna())

    results[cohort] = r

# =====================================================================
print()
print("=" * 78)
print(" WHAT EACH FILE CONTAINS")
print("=" * 78)
print(f"\n {'':<24}{'BRCA':>12}{'OV':>12}")
for label, key in [("mutation rows", "mutations_rows"),
                   ("  non-synonymous", "mutations_nonsyn"),
                   ("samples with mutations", "mutations_samples"),
                   ("samples with expression", "expression_samples"),
                   ("samples with copy number", "cna_samples"),
                   ("patients in clinical", "clinical_patients"),
                   ("genes in expression", "expression_genes")]:
    b = results["brca"].get(key, "-")
    o = results["ov"].get(key, "-")
    b = f"{b:,}" if isinstance(b, (int,)) else b
    o = f"{o:,}" if isinstance(o, (int,)) else o
    print(f" {label:<24}{b:>12}{o:>12}")

print()
print("=" * 78)
print(" INTERSECTIONS - what the analysis could actually use")
print("=" * 78)

for cohort in ["brca", "ov"]:
    r = results[cohort]
    mut = r.get("mut_ids", set())
    expr = r.get("expr_ids", set())
    cna = r.get("cna_ids", set())

    both = mut & expr
    all3 = mut & expr & cna

    print(f"\n  {cohort.upper()}")
    print(f"    mutations only              {len(mut):>6}")
    print(f"    expression only             {len(expr):>6}")
    print(f"    copy number only            {len(cna):>6}")
    print(f"    mutations AND expression    {len(both):>6}"
          f"   <- Line 1, Line 2")
    print(f"    all three                   {len(all3):>6}"
          f"   <- Line 3, convergence")
    print(f"    lost to missing expression  {len(mut - expr):>6}")

    r["usable_line12"] = len(both)
    r["usable_all"] = len(all3)

print()
print("=" * 78)
print(" NUMBERS FOR THE SLIDE")
print("=" * 78)
b, o = results["brca"], results["ov"]
print(f"""
  TCGA PanCancer Atlas 2018, downloaded from cBioPortal

                        BRCA        OV
  patients analysed   {b['usable_line12']:>6}    {o['usable_line12']:>6}
  somatic mutations   {b['mutations_rows']:>6,}    {o['mutations_rows']:>6,}
    non-synonymous    {b['mutations_nonsyn']:>6,}    {o['mutations_nonsyn']:>6,}

  combined: {b['usable_line12'] + o['usable_line12']} patients,
            {b['mutations_rows'] + o['mutations_rows']:,} mutations
""")

out = os.path.expanduser("~/immune_escape_project/results/cohort_counts.tsv")
os.makedirs(os.path.dirname(out), exist_ok=True)
rows = []
for cohort in ["brca", "ov"]:
    r = {k: v for k, v in results[cohort].items() if not k.endswith("_ids")}
    r["cohort"] = cohort
    rows.append(r)
pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
print(f"written to {out}")
