#!/usr/bin/env python3
"""
Build the simulation gene panel from the Line 1 cohort results.

The panel defines which slice of the genome the simulated patients carry.
It is assembled from measured data rather than picked by hand:

  DRIVERS    - genes with the LOWEST silencing percentage. The tumour is
               forced to keep expressing these, so mutations in them can
               produce neoantigens.
  PASSENGERS - genes with the HIGHEST silencing percentage. Mutations here
               are transcriptionally invisible.
  MACHINERY  - the antigen presentation apparatus: the HLA region on chr6,
               B2M on chr15, TAP/PSMB/NLRC5. Needed for OptiType, LOHHLA
               and for Line 2 / Line 3.

Only EXONS are included, not whole gene bodies. Giant genes such as TTN,
DMD and PCDH15 are mostly intron, which no mutation in the MAF touches;
including them would inflate the panel several-fold for nothing. The HLA
locus is the exception and goes in as one continuous block, because HLA
typing and LOH detection need the whole region.

TCGA uses gene symbols that GENCODE has since renamed, so a synonym table
is applied before lookup. Without it four genes are silently dropped.

Coordinates come from GENCODE (GRCh38), matching the reference and the
lifted-over mutation coordinates.

Usage:
  python make_panel.py <results_dir> <gencode_gtf_gz> <out_dir> [top_n] [flank]
"""
import sys, os, gzip, re
import pandas as pd
import numpy as np

RES_DIR = sys.argv[1]
GTF = sys.argv[2]
OUT_DIR = sys.argv[3]
TOP_N = int(sys.argv[4]) if len(sys.argv) > 4 else 20
FLANK = int(sys.argv[5]) if len(sys.argv) > 5 else 1000

MIN_MUT_BRCA = 10
MIN_MUT_OV = 5

# TCGA symbol -> current GENCODE symbol
SYNONYMS = {
    "GPR112": "ADGRG4",
    "KIAA1109": "BLTP1",
    "DPCR1": "MUC21",
    "C6orf10": "TSBP1",
    "GPR98": "ADGRV1",
    "EMR3": "ADGRE3",
    "CCDC168": "CFAP418",
    "FAM47C": "FAM47C",
    "RP1L1": "RP1L1",
    "KIAA0754": "MACF1",
}

MACHINERY = ["B2M", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G",
             "TAP1", "TAP2", "TAPBP", "NLRC5", "PSMB8", "PSMB9",
             "CALR", "CANX", "PDIA3", "ERAP1", "ERAP2"]

HLA_REGION = ("chr6", 29600000, 33100000, "HLA_REGION_BLOCK")

MAIN_CHROMS = set(["chr%s" % c for c in list(range(1, 23)) + ["X", "Y"]])

os.makedirs(OUT_DIR, exist_ok=True)

print("[1/5] reading Line 1 results ...", flush=True)
brca = pd.read_csv("%s/brca_all_genes.csv" % RES_DIR)
ov = pd.read_csv("%s/ov_all_genes.csv" % RES_DIR)
brca = brca[brca.Total_Mutations >= MIN_MUT_BRCA]
ov = ov[ov.Total_Mutations >= MIN_MUT_OV]
print("      BRCA genes: %d   OV genes: %d" % (len(brca), len(ov)), flush=True)


def pick(df, high, n):
    return list(df.sort_values(["Silenced_Percent", "Total_Mutations"],
                               ascending=[not high, False]).head(n).Hugo_Symbol)


drv_b = pick(brca, False, TOP_N)
drv_o = pick(ov, False, TOP_N)
pas_b = pick(brca, True, TOP_N)
pas_o = pick(ov, True, TOP_N)

drivers = sorted(set(drv_b) | set(drv_o))
passengers = sorted(set(pas_b) | set(pas_o))
shared_drv = sorted(set(drv_b) & set(drv_o))
shared_pas = sorted(set(pas_b) & set(pas_o))
passengers = [g for g in passengers if g not in set(drivers)]

print("      drivers    : %d  (shared: %d)"
      % (len(drivers), len(shared_drv)), flush=True)
print("      passengers : %d  (shared: %d)"
      % (len(passengers), len(shared_pas)), flush=True)
print("      machinery  : %d" % len(MACHINERY), flush=True)

category = {}
for g in drivers:
    category[g] = "DRIVER"
for g in passengers:
    category[g] = "PASSENGER"
for g in MACHINERY:
    category[g] = "MACHINERY"

shared = set(shared_drv) | set(shared_pas)

# map TCGA symbols to the names GENCODE actually uses
lookup = {}          # gencode symbol -> original TCGA symbol
for g in category:
    lookup[SYNONYMS.get(g, g)] = g
wanted = set(lookup)

renamed = {g: SYNONYMS[g] for g in category if g in SYNONYMS}
if renamed:
    print("      applying %d synonyms:" % len(renamed), flush=True)
    for k, v in sorted(renamed.items()):
        print("        %-12s -> %s" % (k, v), flush=True)

print("[2/5] scanning GENCODE exons for %d genes ..." % len(wanted), flush=True)
found = {}
exons = {}
name_re = re.compile(r'gene_name "([^"]+)"')
type_re = re.compile(r'gene_type "([^"]+)"')

n_ex = 0
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
        gname = m.group(1)
        if gname not in wanted or f[0] not in MAIN_CHROMS:
            continue
        n_ex += 1
        gt = type_re.search(f[8])
        gt = gt.group(1) if gt else "NA"
        found.setdefault(gname, (f[0], f[6], gt))
        exons.setdefault(gname, []).append((int(f[3]), int(f[4])))

print("      exon records kept: %d" % n_ex, flush=True)
print("      matched genes    : %d / %d" % (len(found), len(wanted)), flush=True)

missing = sorted(lookup[g] for g in (wanted - set(found)))
if missing:
    print("      STILL NOT FOUND (%d): %s"
          % (len(missing), ", ".join(missing)), flush=True)

print("[3/5] building exon intervals ...", flush=True)
rows = []
for gname, (c, strand, gt) in found.items():
    orig = lookup[gname]
    iv = sorted(exons[gname])
    cs = ce = None
    merged_ex = []
    for s, e in iv:
        s0 = max(0, s - 1 - FLANK)
        e0 = e + FLANK
        if cs is None:
            cs, ce = s0, e0
        elif s0 <= ce:
            ce = max(ce, e0)
        else:
            merged_ex.append((cs, ce))
            cs, ce = s0, e0
    if cs is not None:
        merged_ex.append((cs, ce))
    for s0, e0 in merged_ex:
        rows.append({
            "chrom": c, "start": s0, "end": e0, "name": orig,
            "gencode_name": gname,
            "category": category[orig],
            "shared_both_cohorts": orig in shared,
            "gene_type": gt, "strand": strand,
            "length_bp": e0 - s0,
        })

c, s, e, nm = HLA_REGION
rows.append({"chrom": c, "start": s, "end": e, "name": nm,
             "gencode_name": nm, "category": "MACHINERY",
             "shared_both_cohorts": False, "gene_type": "region",
             "strand": ".", "length_bp": e - s})

panel = pd.DataFrame(rows)


def chrom_key(c):
    v = c.replace("chr", "")
    return (int(v) if v.isdigit() else {"X": 23, "Y": 24}.get(v, 99))


panel["_k"] = panel.chrom.map(chrom_key)
panel = panel.sort_values(["_k", "start"]).drop(columns="_k").reset_index(drop=True)

print("[4/5] merging overlaps across genes ...", flush=True)
merged = []
for c, g in panel.groupby("chrom", sort=False):
    g = g.sort_values("start")
    cs = ce = None
    names = []
    for _, r in g.iterrows():
        if cs is None:
            cs, ce, names = r.start, r.end, [r["name"]]
        elif r.start <= ce:
            ce = max(ce, r.end)
            if r["name"] not in names:
                names.append(r["name"])
        else:
            merged.append((c, cs, ce, "+".join(names)))
            cs, ce, names = r.start, r.end, [r["name"]]
    if cs is not None:
        merged.append((c, cs, ce, "+".join(names)))

mdf = pd.DataFrame(merged, columns=["chrom", "start", "end", "name"])
mdf["_k"] = mdf.chrom.map(chrom_key)
mdf = mdf.sort_values(["_k", "start"]).drop(columns="_k")

print("[5/5] writing ...", flush=True)
mdf.to_csv("%s/panel.bed" % OUT_DIR, sep="\t", header=False, index=False)
panel.to_csv("%s/panel_annotated.tsv" % OUT_DIR, sep="\t", index=False)

pd.DataFrame({"Hugo_Symbol": drivers,
              "shared": [g in set(shared_drv) for g in drivers]}
             ).to_csv("%s/panel_drivers.csv" % OUT_DIR, index=False)
pd.DataFrame({"Hugo_Symbol": passengers,
              "shared": [g in set(shared_pas) for g in passengers]}
             ).to_csv("%s/panel_passengers.csv" % OUT_DIR, index=False)

total_bp = int((mdf.end - mdf.start).sum())
by_cat = panel.groupby("category").agg(
    intervals=("name", "size"), bp=("length_bp", "sum"))
by_cat["Mb"] = (by_cat.bp / 1e6).round(2)
by_cat = by_cat.drop(columns="bp")

pd.set_option("display.width", 200)
print()
print("============ PANEL SUMMARY ============")
print("BED intervals    : %d" % len(mdf))
print("total size       : %.2f Mb" % (total_bp / 1e6))
print("genes included   : %d" % len(found))
print()
print(by_cat.to_string())
print()
print("--- DRIVERS shared by both cancers ---")
print("  %s" % (", ".join(shared_drv) if shared_drv else "(none)"))
print()
print("--- PASSENGERS shared by both cancers ---")
print("  %s" % (", ".join(shared_pas) if shared_pas else "(none)"))
print()
print("output: %s/panel.bed" % OUT_DIR)
