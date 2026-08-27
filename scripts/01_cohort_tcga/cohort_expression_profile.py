#!/usr/bin/env python3
"""
Join each simulated sample back to its donor's real RNA and copy number.

Every synthetic patient carries one real TCGA case's mutations, so that
case's expression can be pulled from the same cBioPortal tables and placed
beside the pipeline output. This produces, per sample:

  Line 1   what fraction of the donor's panel mutations sit in genes the
           donor's own tumour does not transcribe
  Line 2   MHC-I machinery expression, as a z-score against the cohort
  Line 3   copy number status of the same machinery genes
  context  tumour mutational burden, cytolytic activity, checkpoints

And the question the simulation exists to ask: of the peptides mhcflurry
called strong binders, how many come from genes that are switched off in
the donor's tumour and therefore cannot produce a protein at all.

Note on the HLA genotype. Peptide binding was predicted against the
backbone genome's HLA, not the donor's - the two are different people.
That is fine for measuring what the pipeline detects, but it means the
neoantigen lists are not the donor's real neoantigens.

Usage:
  python cohort_expression_profile.py <manifest> <cohort_dir> <tcga_root> <out_dir>
"""
import sys, os, glob
import numpy as np
import pandas as pd

MANIFEST = sys.argv[1]
COHORT_DIR = sys.argv[2]
TCGA_ROOT = sys.argv[3]
OUT_DIR = sys.argv[4]

RSEM_OFF = 5.0
MHC = ["B2M", "HLA-A", "HLA-B", "HLA-C", "TAP1", "TAP2", "TAPBP", "NLRC5"]
CHECKPOINT = ["CD274", "PDCD1", "CTLA4", "LAG3", "IDO1", "HAVCR2"]
CYT = ["GZMA", "PRF1"]
TCELL = ["CD8A", "CD3E", "GZMB", "IFNG"]

NONSYN = {"Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
          "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Splice_Site"}

os.makedirs(OUT_DIR, exist_ok=True)


def norm_id(s):
    p = str(s).split("-")
    return "-".join(p[:4])[:15] if len(p) >= 4 else str(s)


def zrow(df, genes):
    """Cohort-wide z-score per gene, averaged over the set."""
    present = [g for g in genes if g in df.index]
    if not present:
        return None
    sub = np.log2(df.loc[present] + 1)
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1) + 1e-9, axis=0)
    return z.mean(axis=0)


print("[1/5] reading the manifest ...", flush=True)
man = pd.read_csv(MANIFEST, sep="\t")
print("      samples: %d" % len(man), flush=True)

data = {}
for cohort in ["brca", "ov"]:
    print("[2/5] loading %s expression and copy number ..." % cohort.upper(),
          flush=True)
    expr = pd.read_csv("%s/%s/data_mrna_seq_v2_rsem.txt" % (TCGA_ROOT, cohort),
                       sep="\t", low_memory=False)
    expr = expr[expr.Hugo_Symbol.notna()].drop_duplicates("Hugo_Symbol")
    expr = expr.set_index("Hugo_Symbol").drop(columns=["Entrez_Gene_Id"],
                                              errors="ignore")
    expr.columns = [norm_id(c) for c in expr.columns]
    expr = expr.loc[:, ~expr.columns.duplicated()].clip(lower=0)

    cna = None
    p = "%s/%s/data_cna.txt" % (TCGA_ROOT, cohort)
    if os.path.exists(p):
        cna = pd.read_csv(p, sep="\t", low_memory=False)
        cna = cna[cna.Hugo_Symbol.notna()].drop_duplicates("Hugo_Symbol")
        cna = cna.set_index("Hugo_Symbol").drop(columns=["Entrez_Gene_Id"],
                                                errors="ignore")
        cna.columns = [norm_id(c) for c in cna.columns]
        cna = cna.loc[:, ~cna.columns.duplicated()]

    cl = None
    p = "%s/%s/data_clinical_patient.txt" % (TCGA_ROOT, cohort)
    if os.path.exists(p):
        cl = pd.read_csv(p, sep="\t", comment="#", low_memory=False)

    data[cohort] = {
        "expr": expr, "cna": cna, "clinical": cl,
        "mhc_z": zrow(expr, MHC),
        "ckpt_z": zrow(expr, CHECKPOINT),
        "cyt_z": zrow(expr, CYT),
        "tcell_z": zrow(expr, TCELL),
        "cohort_median": expr.median(axis=1),
    }
    print("      genes %d, samples %d" % expr.shape, flush=True)

print("[3/5] building the per-sample profile ...", flush=True)
rows = []
for _, m in man.iterrows():
    sid = m.sample_id
    cohort = m.cohort
    sample = norm_id(m.tcga_barcode)
    d = data[cohort]
    expr = d["expr"]

    r = {"sample_id": sid, "cohort": cohort, "donor": sample,
         "backbone": m.backbone, "superpop": m.superpopulation,
         "panel_mutations": m.n_panel_mut, "total_mutations": m.n_total_mut}

    if sample not in expr.columns:
        r["in_expression"] = False
        rows.append(r)
        continue
    r["in_expression"] = True
    e = expr[sample]

    # ---- clinical ----
    if d["clinical"] is not None:
        pid = "-".join(sample.split("-")[:3])
        cl = d["clinical"]
        idc = cl.columns[0]
        row = cl[cl[idc].astype(str).str.startswith(pid)]
        if len(row):
            for col, key in [("SUBTYPE", "subtype"), ("AGE", "age"),
                             ("AJCC_PATHOLOGIC_TUMOR_STAGE", "stage"),
                             ("GENETIC_ANCESTRY_LABEL", "donor_ancestry")]:
                if col in row.columns:
                    r[key] = row.iloc[0][col]

    # ---- Line 1: silencing in the donor's own tumour ----
    truth = "%s/%s/truth_set.tsv" % (COHORT_DIR, sid)
    if os.path.exists(truth):
        t = pd.read_csv(truth, sep="\t")
        t = t[t.Hugo_Symbol.isin(expr.index)]
        if len(t):
            rsem = t.Hugo_Symbol.map(e)
            r["injected_mutations"] = len(t)
            r["silenced_in_donor"] = int((rsem < RSEM_OFF).sum())
            r["silencing_pct"] = round(100 * (rsem < RSEM_OFF).mean(), 1)
            r["median_VAF"] = round(t.VAF.median(), 3)

    # ---- Line 2: MHC-I expression ----
    for g in MHC:
        if g in expr.index:
            r["RSEM_" + g] = round(float(e[g]), 1)
    for key, z in [("mhc_z", "mhc_z"), ("checkpoint_z", "ckpt_z"),
                   ("cytolytic_z", "cyt_z"), ("tcell_z", "tcell_z")]:
        s = d[z]
        if s is not None and sample in s.index:
            r[key] = round(float(s[sample]), 3)

    # ---- Line 3: copy number of the machinery ----
    if d["cna"] is not None and sample in d["cna"].columns:
        c = d["cna"][sample]
        present = [g for g in MHC if g in c.index]
        r["machinery_genes_lost"] = int((c.reindex(present) < 0).sum())
        r["genome_loss_fraction"] = round(float((c < 0).mean()), 3)
        for g in ["B2M", "HLA-A", "HLA-B", "HLA-C"]:
            if g in c.index:
                r["CN_" + g] = int(c[g])

    # ---- pipeline: what Mutect2 recovered ----
    comp = "%s/%s/comparison/truth_vs_calls.tsv" % (COHORT_DIR, sid)
    if os.path.exists(comp):
        cc = pd.read_csv(comp, sep="\t")
        r["mutect2_recall_pct"] = round(100 * cc.detected_pass.mean(), 1)

    # ---- pipeline: neoantigens, and how many are phantoms ----
    neo = "%s/%s/neoantigens/neoantigens_per_mutation.tsv" % (COHORT_DIR, sid)
    if os.path.exists(neo):
        n = pd.read_csv(neo, sep="\t")
        n["donor_rsem"] = n.gene.map(e)
        n["expressed"] = n.donor_rsem >= RSEM_OFF
        tot = int(n.strong.sum())
        real = int(n[n.expressed].strong.sum())
        r["strong_binders"] = tot
        r["strong_expressed"] = real
        r["strong_phantom"] = tot - real
        r["phantom_pct"] = round(100 * (tot - real) / tot, 1) if tot else 0.0

    # ---- HLA from the backbone ----
    hla = "%s/%s/optitype/%s_result.tsv" % (COHORT_DIR, sid, sid)
    if os.path.exists(hla):
        h = pd.read_csv(hla, sep="\t")
        if len(h):
            a1, a2, b1, b2, c1, c2 = [str(h.iloc[0][k])
                                      for k in ["A1", "A2", "B1", "B2", "C1", "C2"]]
            r["backbone_HLA"] = "%s/%s %s/%s %s/%s" % (a1, a2, b1, b2, c1, c2)
            het = ""
            for locus, (p_, q_) in zip("ABC", [(a1, a2), (b1, b2), (c1, c2)]):
                if str(p_) != str(q_):
                    het += locus
            r["het_loci"] = het if het else "-"

    rows.append(r)

prof = pd.DataFrame(rows)
print("      rows: %d" % len(prof), flush=True)

print("[4/5] writing ...", flush=True)
prof.to_csv("%s/cohort_profile.tsv" % OUT_DIR, sep="\t", index=False)

print("[5/5] summary", flush=True)
pd.set_option("display.width", 250)
print()
print("=" * 78)
print(" COHORT PROFILE")
print("=" * 78)

cols = [c for c in ["sample_id", "cohort", "donor", "subtype", "backbone",
                    "superpop", "injected_mutations", "silencing_pct",
                    "median_VAF", "mutect2_recall_pct", "strong_binders",
                    "phantom_pct", "het_loci"] if c in prof.columns]
print(prof[cols].to_string(index=False))

print()
print("--- LINE 1: silencing in the donor's own tumour ---")
if "silencing_pct" in prof.columns:
    for c, g in prof.groupby("cohort"):
        v = g.silencing_pct.dropna()
        if len(v):
            print("  %-5s n=%2d   median %.1f%%   range %.1f - %.1f"
                  % (c.upper(), len(v), v.median(), v.min(), v.max()))

print()
print("--- LINE 2: MHC-I expression, z against the cohort ---")
mz = [c for c in ["mhc_z", "cytolytic_z", "checkpoint_z", "tcell_z"]
      if c in prof.columns]
if mz:
    print(prof.groupby("cohort")[mz].median().round(3).to_string())

print()
print("--- LINE 3: copy number of the machinery ---")
if "machinery_genes_lost" in prof.columns:
    print(prof.groupby("cohort")[["machinery_genes_lost",
                                  "genome_loss_fraction"]]
          .median().round(3).to_string())

print()
print("--- PIPELINE ---")
pp = [c for c in ["mutect2_recall_pct", "strong_binders", "strong_expressed",
                  "strong_phantom", "phantom_pct"] if c in prof.columns]
if pp:
    print(prof.groupby("cohort")[pp].median().round(2).to_string())

if "strong_binders" in prof.columns:
    tot = prof.strong_binders.sum()
    ph = prof.strong_phantom.sum()
    print()
    print("  strong binders across the cohort : %d" % tot)
    print("  from genes silent in the donor   : %d  (%.1f%%)"
          % (ph, 100 * ph / tot if tot else 0))
    print("  those peptides cannot exist - the gene makes no protein")

print()
print("--- MHC-I EXPRESSION, RAW RSEM ---")
rc = [c for c in prof.columns if c.startswith("RSEM_")]
if rc:
    print(prof.groupby("cohort")[rc].median().round(0).to_string())

print()
print("written to %s/cohort_profile.tsv" % OUT_DIR)
