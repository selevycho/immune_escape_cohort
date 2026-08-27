#!/usr/bin/env python3
"""
Verify that every mutation in the truth set is actually present in the
tumour BAM and absent from the normal.

This checks the simulation itself rather than any downstream caller. For
each injected position it counts, in both BAMs: total depth, reads carrying
the alternate base, and the resulting observed VAF.

A mutation is only usable as ground truth if the alternate base is present
in the tumour and absent in the normal. Anything else is a failed injection
and must be excluded before recall is computed, otherwise the caller gets
blamed for something it never had a chance to find.

Reads are counted the way GATK counts them - duplicates, secondary and
supplementary alignments excluded - and low mapping quality is reported
separately, since that is what silently removes evidence.

Usage:
  python verify_injection.py <truth_tsv> <normal_bam> <tumour_bam> <out_dir>
"""
import sys, os
import numpy as np
import pandas as pd
import pysam

TRUTH = sys.argv[1]
NORMAL_BAM = sys.argv[2]
TUMOR_BAM = sys.argv[3]
OUT_DIR = sys.argv[4]

MIN_MAPQ = 20
MIN_BASEQ = 10
MHC = ("chr6", 29600000, 33100000)

os.makedirs(OUT_DIR, exist_ok=True)

print("[1/4] reading truth set ...", flush=True)
truth = pd.read_csv(TRUTH, sep="\t", low_memory=False)
truth = truth[truth.Variant_Type == "SNP"].copy()
print("      SNVs to verify: %d" % len(truth), flush=True)


def count_at(bam, chrom, pos1, ref, alt):
    out = {"depth_hq": 0, "alt_hq": 0, "ref_hq": 0,
           "depth_all": 0, "alt_all": 0, "lowmapq": 0}
    for col in bam.pileup(chrom, pos1 - 1, pos1,
                          truncate=True, stepper="nofilter",
                          min_base_quality=0, max_depth=100000):
        if col.reference_pos != pos1 - 1:
            continue
        for read in col.pileups:
            if read.is_del or read.is_refskip:
                continue
            aln = read.alignment
            if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
                continue
            if aln.is_duplicate:
                continue
            base = aln.query_sequence[read.query_position].upper()
            bq = aln.query_qualities[read.query_position]
            out["depth_all"] += 1
            if base == alt:
                out["alt_all"] += 1
            if aln.mapping_quality < MIN_MAPQ:
                out["lowmapq"] += 1
                continue
            if bq < MIN_BASEQ:
                continue
            out["depth_hq"] += 1
            if base == alt:
                out["alt_hq"] += 1
            elif base == ref:
                out["ref_hq"] += 1
    return out


print("[2/4] counting reads at each position ...", flush=True)
nb = pysam.AlignmentFile(NORMAL_BAM, "rb")
tb = pysam.AlignmentFile(TUMOR_BAM, "rb")

rows = []
for i, (_, r) in enumerate(truth.iterrows()):
    c = r.Chromosome_hg38
    p = int(r.Start_Position_hg38)
    ref = str(r.Reference_Allele).upper()
    alt = str(r.Tumor_Seq_Allele2).upper()

    n = count_at(nb, c, p, ref, alt)
    t = count_at(tb, c, p, ref, alt)

    rows.append({
        "gene": r.Hugo_Symbol, "chrom": c, "pos": p,
        "ref": ref, "alt": alt, "target_vaf": r.VAF,
        "normal_depth": n["depth_hq"], "normal_alt": n["alt_hq"],
        "tumour_depth": t["depth_hq"], "tumour_alt": t["alt_hq"],
        "tumour_depth_all": t["depth_all"], "tumour_alt_all": t["alt_all"],
        "tumour_lowmapq": t["lowmapq"],
        "normal_depth_all": n["depth_all"], "normal_lowmapq": n["lowmapq"],
    })
    if (i + 1) % 40 == 0:
        print("      %d / %d" % (i + 1, len(truth)), flush=True)

nb.close()
tb.close()

df = pd.DataFrame(rows)
df["obs_vaf"] = np.where(df.tumour_depth > 0,
                         df.tumour_alt / df.tumour_depth, np.nan)
df["obs_vaf_all"] = np.where(df.tumour_depth_all > 0,
                             df.tumour_alt_all / df.tumour_depth_all, np.nan)
df["normal_vaf"] = np.where(df.normal_depth > 0,
                            df.normal_alt / df.normal_depth, np.nan)
df["vaf_error"] = df.obs_vaf - df.target_vaf
df["in_mhc"] = ((df.chrom == MHC[0]) &
                (df.pos >= MHC[1]) & (df.pos < MHC[2]))
df["lowmapq_frac"] = np.where(df.tumour_depth_all > 0,
                              df.tumour_lowmapq / df.tumour_depth_all, np.nan)

print("[3/4] classifying ...", flush=True)


def verdict(r):
    if r.tumour_depth_all == 0:
        return "NO_READS"
    if r.tumour_alt_all == 0:
        return "NOT_INJECTED"
    if r.tumour_depth < 5:
        return "LOW_USABLE_DEPTH"
    if r.tumour_alt == 0:
        return "INJECTED_BUT_UNUSABLE"
    if r.normal_alt > 0 and r.normal_vaf > 0.05:
        return "PRESENT_IN_NORMAL"
    return "OK"


df["status"] = df.apply(verdict, axis=1)

print("[4/4] writing ...", flush=True)
df.to_csv("%s/injection_verification.tsv" % OUT_DIR, sep="\t", index=False)
df[df.status != "OK"].to_csv("%s/injection_problems.tsv" % OUT_DIR,
                             sep="\t", index=False)

pd.set_option("display.width", 240)
print()
print("========== INJECTION VERIFICATION ==========")
print("positions checked : %d" % len(df))
print()
print("--- STATUS ---")
print(df.status.value_counts().to_string())
print()
n_ok = int((df.status == "OK").sum())
print("usable as ground truth: %d / %d  (%.1f%%)"
      % (n_ok, len(df), 100.0 * n_ok / len(df)))
print()

print("--- INSIDE vs OUTSIDE THE MHC REGION ---")
g = df.groupby("in_mhc").agg(
    n=("gene", "size"),
    ok=("status", lambda s: int((s == "OK").sum())),
    mean_depth_all=("tumour_depth_all", "mean"),
    mean_depth_usable=("tumour_depth", "mean"),
    mean_lowmapq_frac=("lowmapq_frac", "mean"),
).reset_index()
g["ok_pct"] = (100 * g.ok / g.n).round(1)
print(g.to_string(index=False, float_format=lambda x: "%.3f" % x))
print()

ok = df[df.status == "OK"]
if len(ok) > 2:
    from scipy import stats
    r, p = stats.pearsonr(ok.target_vaf, ok.obs_vaf)
    print("--- VAF ACCURACY (usable positions only) ---")
    print("  n = %d" % len(ok))
    print("  Pearson r        = %.4f  (p = %.3g)" % (r, p))
    print("  mean error       = %+.4f" % ok.vaf_error.mean())
    print("  median abs error = %.4f" % ok.vaf_error.abs().median())
    print("  within +/-0.05   = %.1f%%"
          % (100.0 * (ok.vaf_error.abs() <= 0.05).mean()))
    print("  mean usable depth= %.1f" % ok.tumour_depth.mean())
    print()

print("--- NORMAL SAMPLE CONTAMINATION CHECK ---")
print("  positions with any alt read in normal: %d"
      % int((df.normal_alt > 0).sum()))
print("  positions with normal VAF above 5 pct: %d"
      % int((df.normal_vaf > 0.05).sum()))
print("  (both should be near zero - the normal is the un-mutated backbone)")
print()

bad = df[df.status != "OK"]
if len(bad):
    print("--- PROBLEM POSITIONS (%d) ---" % len(bad))
    cols = ["gene", "chrom", "pos", "target_vaf", "in_mhc",
            "tumour_depth_all", "tumour_depth", "tumour_alt",
            "tumour_lowmapq", "status"]
    print(bad[cols].sort_values(["in_mhc", "pos"])
          .to_string(index=False, float_format=lambda x: "%.3f" % x))
    print()

print("--- 10 LARGEST VAF DEVIATIONS AMONG USABLE POSITIONS ---")
if len(ok):
    cols = ["gene", "chrom", "pos", "target_vaf", "obs_vaf",
            "vaf_error", "tumour_depth", "tumour_alt"]
    print(ok.reindex(ok.vaf_error.abs().sort_values(ascending=False).index)
          .head(10)[cols].to_string(index=False, float_format=lambda x: "%.3f" % x))
print()
print("files written to %s" % OUT_DIR)
