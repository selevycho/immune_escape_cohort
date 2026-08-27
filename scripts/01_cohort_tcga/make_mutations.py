#!/usr/bin/env python3
"""
Turn one real TCGA patient's mutations into BAMSurgeon input files.

The simulated patient reproduces a real TCGA case: same genes, same base
changes, same variant allele frequencies. Nothing is invented.

Output formats follow BAMSurgeon's variant file convention:
  SNV   : chrom  pos  pos  VAF  altbase
  INDEL : chrom  start  end  VAF  INS|DEL  [bases]

Also writes a truth set, which is what Mutect2 output gets compared against.

Usage:
  python make_mutations.py <lifted_maf> <panel_bed> <barcode> <out_dir> [min_vaf]
"""
import sys, os
import numpy as np
import pandas as pd

MAF_IN = sys.argv[1]
PANEL = sys.argv[2]
BARCODE = sys.argv[3]
OUT_DIR = sys.argv[4]
MIN_VAF = float(sys.argv[5]) if len(sys.argv) > 5 else 0.05

# The HLA locus is extremely polymorphic; read reassembly there is unreliable.
# Set to True to skip mutations inside it.
SKIP_HLA = False
HLA = ("chr6", 29600000, 33100000)

NONSYN = {"Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
          "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
          "Splice_Site", "Nonstop_Mutation", "Translation_Start_Site"}

BASES = {"A", "C", "G", "T"}

os.makedirs(OUT_DIR, exist_ok=True)

print("[1/6] reading lifted MAF ...", flush=True)
maf = pd.read_csv(MAF_IN, sep="\t", low_memory=False)
print("      rows: %d" % len(maf), flush=True)

print("[2/6] selecting patient %s ..." % BARCODE, flush=True)
sub = maf[maf.Tumor_Sample_Barcode == BARCODE].copy()
if len(sub) == 0:
    cands = maf.Tumor_Sample_Barcode.unique()[:5]
    raise SystemExit("barcode not found. examples: %s" % list(cands))
print("      mutations for this patient: %d" % len(sub), flush=True)

sub = sub[sub.Variant_Classification.isin(NONSYN)]
print("      non-synonymous            : %d" % len(sub), flush=True)

print("[3/6] restricting to the panel ...", flush=True)
panel = pd.read_csv(PANEL, sep="\t", header=None,
                    names=["chrom", "start", "end", "name"])
by_chrom = {}
for c, g in panel.groupby("chrom"):
    by_chrom[c] = g[["start", "end"]].to_numpy()


def in_panel(c, pos):
    iv = by_chrom.get(c)
    if iv is None:
        return False
    p = pos - 1                       # BED is 0-based half-open
    return bool(((iv[:, 0] <= p) & (p < iv[:, 1])).any())


keep = [in_panel(r.Chromosome_hg38, int(r.Start_Position_hg38))
        for _, r in sub.iterrows()]
sub = sub[keep].copy()
print("      inside panel              : %d" % len(sub), flush=True)

if SKIP_HLA:
    c, s, e = HLA
    inhla = ((sub.Chromosome_hg38 == c) &
             (sub.Start_Position_hg38 >= s) &
             (sub.Start_Position_hg38 < e))
    print("      dropping %d HLA-region mutations" % int(inhla.sum()), flush=True)
    sub = sub[~inhla].copy()

print("[4/6] computing VAF ...", flush=True)
den = sub.t_ref_count.fillna(0) + sub.t_alt_count.fillna(0)
sub["VAF"] = np.where(den > 0, sub.t_alt_count / den, np.nan)
sub = sub[sub.VAF.notna()].copy()

n_before = len(sub)
sub = sub[sub.VAF >= MIN_VAF].copy()
print("      dropped %d below VAF %.2f" % (n_before - len(sub), MIN_VAF), flush=True)
sub["VAF"] = sub.VAF.clip(upper=0.95).round(4)

print("[5/6] splitting SNV and indel ...", flush=True)
is_snv = (sub.Reference_Allele.isin(BASES) &
          sub.Tumor_Seq_Allele2.isin(BASES) &
          (sub.Variant_Type == "SNP"))
snv = sub[is_snv].copy()
ind = sub[~is_snv].copy()
print("      SNV   : %d" % len(snv), flush=True)
print("      indel : %d" % len(ind), flush=True)

snv_out = pd.DataFrame({
    "chrom": snv.Chromosome_hg38,
    "start": snv.Start_Position_hg38.astype(int),
    "end": snv.Start_Position_hg38.astype(int),
    "vaf": snv.VAF,
    "alt": snv.Tumor_Seq_Allele2,
}).sort_values(["chrom", "start"])

rows = []
for _, r in ind.iterrows():
    ref = str(r.Reference_Allele)
    alt = str(r.Tumor_Seq_Allele2)
    s = int(r.Start_Position_hg38)
    e = int(r.End_Position_hg38)
    if ref == "-" or r.Variant_Type in ("INS",):
        rows.append((r.Chromosome_hg38, s, s + 1, r.VAF, "INS", alt))
    elif alt == "-" or r.Variant_Type in ("DEL",):
        rows.append((r.Chromosome_hg38, s, e, r.VAF, "DEL", ""))
ind_out = pd.DataFrame(rows, columns=["chrom", "start", "end",
                                      "vaf", "type", "seq"])
if len(ind_out):
    ind_out = ind_out.sort_values(["chrom", "start"])

print("[6/6] writing ...", flush=True)
tag = BARCODE.replace("-", "_")

snv_path = "%s/snv_mutations.txt" % OUT_DIR
ind_path = "%s/indel_mutations.txt" % OUT_DIR
snv_out.to_csv(snv_path, sep="\t", header=False, index=False)
ind_out.to_csv(ind_path, sep="\t", header=False, index=False)

truth = sub[["Hugo_Symbol", "Chromosome_hg38", "Start_Position_hg38",
             "End_Position_hg38", "Reference_Allele", "Tumor_Seq_Allele2",
             "Variant_Type", "Variant_Classification", "VAF",
             "t_ref_count", "t_alt_count"]].copy()
if "HGVSp_Short" in sub.columns:
    truth["HGVSp_Short"] = sub.HGVSp_Short
truth["source_patient"] = BARCODE
truth = truth.sort_values(["Chromosome_hg38", "Start_Position_hg38"])
truth.to_csv("%s/truth_set.tsv" % OUT_DIR, sep="\t", index=False)

pd.set_option("display.width", 220)
print()
print("========== MUTATION SET: %s ==========" % BARCODE)
print("SNVs to inject   : %d" % len(snv_out))
print("indels to inject : %d" % len(ind_out))
print("total            : %d" % (len(snv_out) + len(ind_out)))
print()
print("VAF distribution:")
q = sub.VAF.describe(percentiles=[0.1, 0.25, 0.5, 0.75, 0.9])
print(q.to_string())
print()
bins = [0, 0.1, 0.2, 0.3, 0.5, 1.01]
lab = ["<10%", "10-20%", "20-30%", "30-50%", ">50%"]
vc = pd.cut(sub.VAF, bins=bins, labels=lab, right=False).value_counts().sort_index()
print("VAF bins (these become the sensitivity curve):")
print(vc.to_string())
print()
print("Genes covered (top 15 by mutation count):")
print(sub.Hugo_Symbol.value_counts().head(15).to_string())
print()
print("Chromosomes: %s" % ", ".join(sorted(sub.Chromosome_hg38.unique(),
                                           key=lambda c: c)))
print()
print("files:")
print("  %s" % snv_path)
print("  %s" % ind_path)
print("  %s/truth_set.tsv" % OUT_DIR)
