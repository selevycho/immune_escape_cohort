#!/usr/bin/env python3
"""
Lift TCGA MAF coordinates from GRCh37/hg19 to GRCh38/hg38.

Why this is needed: cBioPortal MAF files use GRCh37 positions, while the
1000 Genomes 30x alignments and the reference on this cluster are GRCh38.
Injecting hg19 coordinates into an hg38 BAM places every mutation in the
wrong place, and no downstream tool reports an error.

Uses pyliftover (pure Python) rather than the UCSC binary, which requires
an OpenSSL version not present on RHEL 8.

The final step verifies each lifted SNV against the hg38 FASTA: if the base
at the new position does not match Reference_Allele, the lift is wrong.
This is the only way to catch a silent failure.

Usage:
  python liftover_maf.py <maf_in> <chain> <hg38_fasta> <out_prefix>
"""
import sys, os
import pandas as pd
import numpy as np

MAF_IN = sys.argv[1]
CHAIN = sys.argv[2]
FASTA = sys.argv[3]
OUT_PREFIX = sys.argv[4]

BASES = {"A", "C", "G", "T"}

print("[1/6] reading MAF ...", flush=True)
maf = pd.read_csv(MAF_IN, sep="\t", low_memory=False, comment="#")
print("      rows: %d" % len(maf), flush=True)

need = ["Hugo_Symbol", "Chromosome", "Start_Position", "End_Position",
        "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Classification",
        "Variant_Type", "Tumor_Sample_Barcode"]
missing = [c for c in need if c not in maf.columns]
if missing:
    raise SystemExit("MAF is missing columns: %s" % missing)

build_col = [c for c in maf.columns if c.upper() == "NCBI_BUILD"]
if build_col:
    print("      NCBI_Build: %s"
          % maf[build_col[0]].astype(str).unique()[:5], flush=True)

maf = maf.reset_index(drop=True)

print("[2/6] loading chain file ...", flush=True)
from pyliftover import LiftOver
lo = LiftOver(CHAIN)
print("      chain loaded", flush=True)

print("[3/6] converting coordinates ...", flush=True)
chrom19 = ("chr" + maf.Chromosome.astype(str)
           .str.replace("^chr", "", regex=True)).to_numpy()
start1 = maf.Start_Position.astype(np.int64).to_numpy()
end1 = maf.End_Position.astype(np.int64).to_numpy()
span = end1 - start1

new_chrom = np.empty(len(maf), dtype=object)
new_start = np.full(len(maf), -1, dtype=np.int64)
n_unmapped = 0

for i in range(len(maf)):
    # pyliftover works in 0-based coordinates
    res = lo.convert_coordinate(chrom19[i], int(start1[i]) - 1)
    if not res:
        n_unmapped += 1
        continue
    c, p = res[0][0], res[0][1]
    new_chrom[i] = c
    new_start[i] = p + 1                     # back to 1-based
    if (i + 1) % 20000 == 0:
        print("      %d / %d" % (i + 1, len(maf)), flush=True)

maf["Chromosome_hg38"] = new_chrom
maf["Start_Position_hg38"] = new_start
maf["End_Position_hg38"] = new_start + span

print("      unmapped: %d  (%.2f%%)"
      % (n_unmapped, 100.0 * n_unmapped / len(maf)), flush=True)

print("[4/6] filtering ...", flush=True)
out = maf[maf.Start_Position_hg38 > 0].copy()

same_chr = (out.Chromosome_hg38.str.replace("^chr", "", regex=True) ==
            out.Chromosome.astype(str).str.replace("^chr", "", regex=True))
n_jumped = int((~same_chr).sum())
if n_jumped:
    print("      WARNING: %d mutations changed chromosome - dropped"
          % n_jumped, flush=True)
out = out[same_chr].copy()
print("      kept: %d" % len(out), flush=True)

print("[5/6] verifying reference bases against hg38 FASTA ...", flush=True)
import pysam
fa = pysam.FastaFile(FASTA)
if "chr1" not in set(fa.references):
    print("      FASTA has no 'chr' prefix - stripping", flush=True)
    out["Chromosome_hg38"] = out.Chromosome_hg38.str.replace("^chr", "",
                                                             regex=True)

snv = out[out.Reference_Allele.isin(BASES) &
          out.Tumor_Seq_Allele2.isin(BASES)]
check = snv.sample(min(len(snv), 5000), random_state=1) if len(snv) else snv

ok, bad, err = 0, 0, 0
bad_rows = []
for _, r in check.iterrows():
    try:
        b = fa.fetch(r.Chromosome_hg38,
                     int(r.Start_Position_hg38) - 1,
                     int(r.Start_Position_hg38)).upper()
    except (KeyError, ValueError):
        err += 1
        continue
    if b == r.Reference_Allele:
        ok += 1
    else:
        bad += 1
        if len(bad_rows) < 10:
            bad_rows.append((r.Hugo_Symbol, r.Chromosome_hg38,
                             int(r.Start_Position_hg38),
                             r.Reference_Allele, b))
fa.close()

tot = ok + bad
pct = 100.0 * ok / tot if tot else float("nan")
print("      SNVs checked: %d   match: %d   mismatch: %d   errors: %d"
      % (tot, ok, bad, err), flush=True)
print("      concordance : %.2f%%" % pct, flush=True)
for row in bad_rows:
    print("        mismatch: %s" % (row,), flush=True)

print("[6/6] writing ...", flush=True)
cols = ["Hugo_Symbol", "Tumor_Sample_Barcode", "Variant_Classification",
        "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele2",
        "Chromosome", "Start_Position", "End_Position",
        "Chromosome_hg38", "Start_Position_hg38", "End_Position_hg38"]
extra = [c for c in ["t_ref_count", "t_alt_count", "t_depth", "HGVSp_Short"]
         if c in out.columns]
out[cols + extra].to_csv(OUT_PREFIX + ".hg38.maf.tsv", sep="\t", index=False)

print()
print("=========== LIFTOVER SUMMARY ===========")
print("input mutations   : %d" % len(maf))
print("lifted to hg38    : %d  (%.2f%%)" % (len(out), 100.0 * len(out) / len(maf)))
print("unmapped          : %d" % n_unmapped)
print("changed chromosome: %d" % n_jumped)
print("reference concordance: %.2f%%" % pct)
print()
if pct >= 99.0:
    print("PASS - coordinates are consistent with the hg38 reference.")
elif pct >= 95.0:
    print("MARGINAL - inspect the mismatches before injecting anything.")
else:
    print("FAIL - do NOT use this output.")
    print("Likely causes: wrong chain file, or a coordinate convention error.")
print()
print("output: %s.hg38.maf.tsv" % OUT_PREFIX)
