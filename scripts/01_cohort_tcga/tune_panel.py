#!/usr/bin/env python3
"""
Explore how panel size trades against mutation coverage.

Widening the panel pulls in more patients with enough mutations to be worth
simulating, but the return diminishes: the genes with the most mutations are
already in, and everything added after them is rarer. This sweeps TOP_N and
reports, for each setting, how large the panel would be and how many
patients would clear the selection threshold.

Nothing is written; this only informs the choice of TOP_N for make_panel.py.

Usage:
  python tune_panel.py <results_dir> <gencode_gtf_gz> <liftover_dir> [thresholds]
"""
import sys, os, gzip, re
import pandas as pd
import numpy as np

RES_DIR = sys.argv[1]
GTF = sys.argv[2]
LIFT_DIR = sys.argv[3]
TOP_N_VALUES = [int(x) for x in sys.argv[4].split(",")] \
    if len(sys.argv) > 4 else [60, 100, 150, 200]

MIN_MUT_BRCA = 10
MIN_MUT_OV = 5
FLANK = 1000
SEL_THRESHOLD = 15          # mutations per patient needed to be usable

SYNONYMS = {"GPR112": "ADGRG4", "KIAA1109": "BLTP1", "DPCR1": "MUC21",
            "C6orf10": "TSBP1", "GPR98": "ADGRV1", "EMR3": "ADGRE3",
            "KIAA0754": "MACF1"}

MACHINERY = ["B2M", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G",
             "TAP1", "TAP2", "TAPBP", "NLRC5", "PSMB8", "PSMB9",
             "CALR", "CANX", "PDIA3", "ERAP1", "ERAP2"]
HLA_BLOCK_MB = 3.5
MAIN_CHROMS = set(["chr%s" % c for c in list(range(1, 23)) + ["X", "Y"]])
NONSYN = {"Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
          "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Splice_Site"}

print("[1/4] reading Line 1 tables ...", flush=True)
brca = pd.read_csv("%s/brca_all_genes.csv" % RES_DIR)
ov = pd.read_csv("%s/ov_all_genes.csv" % RES_DIR)
brca = brca[brca.Total_Mutations >= MIN_MUT_BRCA]
ov = ov[ov.Total_Mutations >= MIN_MUT_OV]
print("      BRCA %d genes, OV %d genes" % (len(brca), len(ov)), flush=True)


def pick(df, high, n):
    return list(df.sort_values(["Silenced_Percent", "Total_Mutations"],
                               ascending=[not high, False]).head(n).Hugo_Symbol)


print("[2/4] indexing GENCODE exons (one pass) ...", flush=True)
# collect exon spans for every gene we might ever want
all_wanted = set(MACHINERY)
for n in TOP_N_VALUES:
    for df in (brca, ov):
        all_wanted |= set(pick(df, False, n)) | set(pick(df, True, n))
gencode_names = {SYNONYMS.get(g, g): g for g in all_wanted}

exons = {}
name_re = re.compile(r'gene_name "([^"]+)"')
with gzip.open(GTF, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if f[2] != "exon":
            continue
        m = name_re.search(f[8])
        if not m:
            continue
        gn = m.group(1)
        if gn not in gencode_names or f[0] not in MAIN_CHROMS:
            continue
        exons.setdefault(gencode_names[gn], []).append(
            (f[0], int(f[3]), int(f[4])))
print("      genes with exons: %d / %d" % (len(exons), len(all_wanted)),
      flush=True)


def gene_span_bp(gene):
    """Merged exon length with flanks, in bp."""
    iv = exons.get(gene)
    if not iv:
        return 0, []
    out = []
    for c, g in pd.DataFrame(iv, columns=["c", "s", "e"]).groupby("c"):
        g = g.sort_values("s")
        cs = ce = None
        for s, e in zip(g.s, g.e):
            s0, e0 = max(0, s - 1 - FLANK), e + FLANK
            if cs is None:
                cs, ce = s0, e0
            elif s0 <= ce:
                ce = max(ce, e0)
            else:
                out.append((c, cs, ce))
                cs, ce = s0, e0
        if cs is not None:
            out.append((c, cs, ce))
    return sum(e - s for _, s, e in out), out


print("[3/4] loading mutation coordinates ...", flush=True)
mafs = {}
for coh in ["brca", "ov"]:
    m = pd.read_csv("%s/%s.hg38.maf.tsv" % (LIFT_DIR, coh), sep="\t",
                    low_memory=False)
    m = m[m.Variant_Classification.isin(NONSYN) & (m.Variant_Type == "SNP")]
    den = m.t_ref_count.fillna(0) + m.t_alt_count.fillna(0)
    m = m[(den > 0) & (m.t_alt_count / den >= 0.05)]
    mafs[coh] = m[["Tumor_Sample_Barcode", "Chromosome_hg38",
                   "Start_Position_hg38"]]
    print("      %s: %d usable mutations, %d patients"
          % (coh.upper(), len(m), m.Tumor_Sample_Barcode.nunique()), flush=True)

print("[4/4] sweeping TOP_N ...", flush=True)
rows = []
for n in TOP_N_VALUES:
    genes = set(MACHINERY)
    for df in (brca, ov):
        genes |= set(pick(df, False, n)) | set(pick(df, True, n))

    intervals = []
    total_bp = 0
    for g in genes:
        bp, iv = gene_span_bp(g)
        total_bp += bp
        intervals.extend(iv)
    total_mb = total_bp / 1e6 + HLA_BLOCK_MB

    by = {}
    for c, s, e in intervals:
        by.setdefault(c, []).append((s, e))
    by = {c: np.array(v) for c, v in by.items()}
    # the HLA block
    by.setdefault("chr6", np.empty((0, 2)))
    by["chr6"] = np.vstack([by["chr6"], [[29600000, 33100000]]])

    row = {"TOP_N": n, "genes": len(genes), "panel_Mb": round(total_mb, 1)}
    for coh in ["brca", "ov"]:
        m = mafs[coh]
        hit = np.zeros(len(m), dtype=bool)
        for i, (c, p) in enumerate(zip(m.Chromosome_hg38.values,
                                       m.Start_Position_hg38.values)):
            iv = by.get(c)
            if iv is None or len(iv) == 0:
                continue
            pp = int(p) - 1
            if ((iv[:, 0] <= pp) & (pp < iv[:, 1])).any():
                hit[i] = True
        cnt = m[hit].groupby("Tumor_Sample_Barcode").size()
        row["%s_mutations" % coh] = int(hit.sum())
        row["%s_patients_ge15" % coh] = int((cnt >= SEL_THRESHOLD).sum())
        row["%s_patients_ge20" % coh] = int((cnt >= 20).sum())
        row["%s_median" % coh] = float(cnt.median()) if len(cnt) else 0.0
    rows.append(row)
    print("      TOP_N=%d done" % n, flush=True)

res = pd.DataFrame(rows)
pd.set_option("display.width", 220)
print()
print("=================== PANEL SIZE SWEEP ===================")
print(res.to_string(index=False))
print()
print("Reading this: the useful column is patients_ge15 - how many cases")
print("carry enough panel mutations to give a meaningful recall estimate.")
print("If it stops climbing while panel_Mb keeps growing, the extra genes")
print("are only costing disk and runtime.")
