#!/usr/bin/env python3
"""
Census of variant types in the two cBioPortal cohorts.

The simulation currently injects SNVs only. Before deciding whether to add
indels or structural variants, it is worth knowing what fraction of the
data each type actually represents - both across the whole cohort and
inside the panel, which is what the simulation can reach at all.

Three views are produced:
  1. all mutations in the MAF, by Variant_Type and Variant_Classification
  2. the same restricted to non-synonymous consequences
  3. the same restricted to the panel, since anything outside it is out of
     scope for the simulated samples regardless of type

Structural variants are counted from data_sv.txt if it is present. That
file is separate from the MAF: fusions, translocations and large
rearrangements are not represented as point mutations and cBioPortal
stores them apart.

Usage:
  python mutation_type_census.py <tcga_dir_brca> <tcga_dir_ov> \
                                 <liftover_dir> <panel_bed> <out_dir>
"""
import sys, os
import numpy as np
import pandas as pd

TCGA_BRCA = sys.argv[1]
TCGA_OV = sys.argv[2]
LIFT_DIR = sys.argv[3]
PANEL = sys.argv[4]
OUT_DIR = sys.argv[5]

MIN_VAF = 0.05
MHC = ("chr6", 29600000, 33100000)

NONSYN = {"Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
          "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
          "Splice_Site", "Nonstop_Mutation", "Translation_Start_Site"}

INDEL_CLASS = {"Frame_Shift_Del", "Frame_Shift_Ins",
               "In_Frame_Del", "In_Frame_Ins"}

os.makedirs(OUT_DIR, exist_ok=True)

print("[1/4] loading the panel ...", flush=True)
panel = pd.read_csv(PANEL, sep="\t", header=None,
                    names=["chrom", "start", "end", "name"])
by_chrom = {c: g[["start", "end"]].to_numpy()
            for c, g in panel.groupby("chrom")}
panel_mb = (panel.end - panel.start).sum() / 1e6
print("      %d intervals, %.1f Mb" % (len(panel), panel_mb), flush=True)


def in_panel_mask(df):
    out = np.zeros(len(df), dtype=bool)
    for i, (c, p) in enumerate(zip(df.Chromosome_hg38.values,
                                   df.Start_Position_hg38.values)):
        iv = by_chrom.get(c)
        if iv is None:
            continue
        pp = int(p) - 1
        out[i] = bool(((iv[:, 0] <= pp) & (pp < iv[:, 1])).any())
    return out


summary_rows = []
detail = {}

for cohort, tcga_dir in [("brca", TCGA_BRCA), ("ov", TCGA_OV)]:
    print("\n[2/4] %s - reading the lifted MAF ..." % cohort.upper(), flush=True)
    maf = pd.read_csv("%s/%s.hg38.maf.tsv" % (LIFT_DIR, cohort),
                      sep="\t", low_memory=False)
    n_all = len(maf)
    n_pat = maf.Tumor_Sample_Barcode.nunique()
    print("      %d mutations, %d patients" % (n_all, n_pat), flush=True)

    maf["is_indel_class"] = maf.Variant_Classification.isin(INDEL_CLASS)
    maf["is_nonsyn"] = maf.Variant_Classification.isin(NONSYN)

    den = maf.t_ref_count.fillna(0) + maf.t_alt_count.fillna(0)
    maf["VAF"] = np.where(den > 0, maf.t_alt_count / den, np.nan)

    # ---- level 1: everything ----
    vt_all = maf.Variant_Type.value_counts()

    # ---- level 2: non-synonymous only ----
    ns = maf[maf.is_nonsyn].copy()
    vt_ns = ns.Variant_Type.value_counts()

    # ---- level 3: inside the panel, VAF filtered ----
    ns_v = ns[ns.VAF >= MIN_VAF].copy()
    mask = in_panel_mask(ns_v)
    inp = ns_v[mask].copy()
    inp["in_mhc"] = ((inp.Chromosome_hg38 == MHC[0]) &
                     (inp.Start_Position_hg38 >= MHC[1]) &
                     (inp.Start_Position_hg38 < MHC[2]))
    vt_panel = inp.Variant_Type.value_counts()

    def pct(series, key):
        tot = series.sum()
        return 100.0 * series.get(key, 0) / tot if tot else 0.0

    n_snp_p = int(vt_panel.get("SNP", 0))
    n_ind_p = int(len(inp) - n_snp_p)

    summary_rows.append({
        "cohort": cohort,
        "patients": n_pat,
        "mutations_all": n_all,
        "SNP_pct_all": round(pct(vt_all, "SNP"), 1),
        "indel_pct_all": round(100 - pct(vt_all, "SNP"), 1),
        "nonsyn": len(ns),
        "SNP_pct_nonsyn": round(pct(vt_ns, "SNP"), 1),
        "indel_pct_nonsyn": round(100 - pct(vt_ns, "SNP"), 1),
        "in_panel": len(inp),
        "SNP_in_panel": n_snp_p,
        "indel_in_panel": n_ind_p,
        "indel_pct_panel": round(100.0 * n_ind_p / len(inp), 1) if len(inp) else 0,
        "indel_in_mhc": int(inp[inp.in_mhc & (inp.Variant_Type != "SNP")].shape[0]),
    })

    detail[cohort] = {
        "vt_all": vt_all, "vt_ns": vt_ns, "vt_panel": vt_panel,
        "vc_ns": ns.Variant_Classification.value_counts(),
        "vc_panel": inp.Variant_Classification.value_counts(),
        "indel_vaf": inp[inp.Variant_Type != "SNP"].VAF,
        "snp_vaf": inp[inp.Variant_Type == "SNP"].VAF,
        "indel_per_patient": inp[inp.Variant_Type != "SNP"]
                             .groupby("Tumor_Sample_Barcode").size(),
    }

    inp.to_csv("%s/%s_panel_mutations.tsv" % (OUT_DIR, cohort),
               sep="\t", index=False)

print("\n[3/4] structural variants ...", flush=True)
sv_rows = []
for cohort, tcga_dir in [("brca", TCGA_BRCA), ("ov", TCGA_OV)]:
    sv_path = "%s/data_sv.txt" % tcga_dir
    if not os.path.exists(sv_path):
        print("      %s: data_sv.txt not present in %s"
              % (cohort.upper(), tcga_dir), flush=True)
        sv_rows.append({"cohort": cohort, "sv_events": 0,
                        "patients_with_sv": 0, "note": "file absent"})
        continue
    sv = pd.read_csv(sv_path, sep="\t", low_memory=False)
    idcol = [c for c in sv.columns if "Sample_Id" in c or c == "Sample_ID"]
    idcol = idcol[0] if idcol else sv.columns[0]
    print("      %s: %d events, %d samples"
          % (cohort.upper(), len(sv), sv[idcol].nunique()), flush=True)
    sv_rows.append({"cohort": cohort, "sv_events": len(sv),
                    "patients_with_sv": int(sv[idcol].nunique()),
                    "note": ""})
    detail[cohort]["sv"] = sv
    detail[cohort]["sv_idcol"] = idcol

print("\n[4/4] writing ...", flush=True)
summ = pd.DataFrame(summary_rows)
summ.to_csv("%s/variant_type_summary.tsv" % OUT_DIR, sep="\t", index=False)
pd.DataFrame(sv_rows).to_csv("%s/sv_summary.tsv" % OUT_DIR,
                             sep="\t", index=False)

pd.set_option("display.width", 240)
print()
print("=" * 72)
print(" VARIANT TYPE CENSUS")
print("=" * 72)

for cohort in ["brca", "ov"]:
    d = detail[cohort]
    r = summ[summ.cohort == cohort].iloc[0]
    print()
    print("--- %s ---" % cohort.upper())
    print("  patients: %d   mutations in the MAF: %d"
          % (r.patients, r.mutations_all))
    print()
    print("  all mutations:")
    print(d["vt_all"].to_string().replace("\n", "\n    ").rjust(4))
    print("    SNP %.1f%%, indel %.1f%%" % (r.SNP_pct_all, r.indel_pct_all))
    print()
    print("  non-synonymous only (%d):" % r.nonsyn)
    print("    SNP %.1f%%, indel %.1f%%"
          % (r.SNP_pct_nonsyn, r.indel_pct_nonsyn))
    print()
    print("  inside the panel, VAF >= %.2f (%d mutations):" % (MIN_VAF, r.in_panel))
    print("    SNP   %6d  (%.1f%%)" % (r.SNP_in_panel, 100 - r.indel_pct_panel))
    print("    indel %6d  (%.1f%%)" % (r.indel_in_panel, r.indel_pct_panel))
    print("    of those indels, %d sit in the MHC where injection fails"
          % r.indel_in_mhc)
    print()
    iv = d["indel_vaf"].dropna()
    sv_ = d["snp_vaf"].dropna()
    if len(iv):
        print("  VAF, indels : median %.3f, %.0f%% below 0.10"
              % (iv.median(), 100 * (iv < 0.10).mean()))
    if len(sv_):
        print("  VAF, SNVs   : median %.3f, %.0f%% below 0.10"
              % (sv_.median(), 100 * (sv_ < 0.10).mean()))
    ipp = d["indel_per_patient"]
    if len(ipp):
        print("  patients carrying at least one panel indel: %d"
              % (ipp > 0).sum())
        print("  indels per such patient: median %d, max %d"
              % (ipp.median(), ipp.max()))

print()
print("=" * 72)
print(" SUMMARY TABLE")
print("=" * 72)
print(summ[["cohort", "patients", "mutations_all", "indel_pct_all",
            "indel_pct_nonsyn", "in_panel", "indel_in_panel",
            "indel_pct_panel"]].to_string(index=False))

print()
print("=" * 72)
print(" STRUCTURAL VARIANTS")
print("=" * 72)
for cohort in ["brca", "ov"]:
    d = detail[cohort]
    if "sv" not in d:
        print("  %s: no data_sv.txt in the study folder" % cohort.upper())
        continue
    sv = d["sv"]
    print("  %s: %d events across %d samples"
          % (cohort.upper(), len(sv), sv[d["sv_idcol"]].nunique()))
    for col in ["Class", "SV_Status", "Event_Info", "Site1_Hugo_Symbol"]:
        if col in sv.columns:
            print("    %s:" % col)
            print(sv[col].value_counts().head(8).to_string()
                  .replace("\n", "\n      ").rjust(6))
    per = sv.groupby(d["sv_idcol"]).size()
    print("    events per sample: median %d, max %d" % (per.median(), per.max()))

print()
print("files written to %s" % OUT_DIR)
