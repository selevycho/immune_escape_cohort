#!/usr/bin/env python3
"""
Compare Mutect2 calls against the injected truth set.

This is the measurement the whole simulation exists for: how many of the
mutations that were actually put into the BAM does the caller recover, and
how does that depend on variant allele frequency?

Every injected mutation is known exactly - position, base change and target
VAF - so recall, precision and the VAF detection curve are all measurable
rather than assumed.

Outputs:
  - per-mutation table with detected / missed status and observed VAF
  - recall broken down by VAF bin (the sensitivity curve)
  - false positives, and what filters Mutect2 applied
  - target VAF vs observed VAF agreement

Usage:
  python compare_mutect2_truth.py <truth_tsv> <filtered_vcf_gz> <out_dir>
"""
import sys, os
import numpy as np
import pandas as pd

TRUTH = sys.argv[1]
VCF = sys.argv[2]
OUT_DIR = sys.argv[3]

VAF_BINS = [0.0, 0.10, 0.20, 0.30, 0.50, 1.01]
VAF_LABELS = ["<10%", "10-20%", "20-30%", "30-50%", ">50%"]

os.makedirs(OUT_DIR, exist_ok=True)

print("[1/6] reading truth set ...", flush=True)
truth = pd.read_csv(TRUTH, sep="\t", low_memory=False)
truth = truth[truth.Variant_Type == "SNP"].copy()
truth["key"] = (truth.Chromosome_hg38.astype(str) + ":" +
                truth.Start_Position_hg38.astype(int).astype(str) + ":" +
                truth.Reference_Allele.astype(str) + ">" +
                truth.Tumor_Seq_Allele2.astype(str))
print("      injected SNVs: %d" % len(truth), flush=True)

print("[2/6] reading Mutect2 VCF ...", flush=True)
import pysam
vcf = pysam.VariantFile(VCF)

samples = list(vcf.header.samples)
print("      samples in VCF: %s" % samples, flush=True)
tumor_sample = None
for s in samples:
    if "TUMOR" in s.upper():
        tumor_sample = s
if tumor_sample is None and samples:
    tumor_sample = samples[0]
print("      using as tumour: %s" % tumor_sample, flush=True)

calls = []
for rec in vcf:
    if rec.alts is None:
        continue
    for alt in rec.alts:
        if len(rec.ref) != 1 or len(alt) != 1:
            continue                      # SNVs only
        filt = list(rec.filter.keys())
        passed = (len(filt) == 0 or filt == ["PASS"])
        af = dp = np.nan
        if tumor_sample and tumor_sample in rec.samples:
            smp = rec.samples[tumor_sample]
            v = smp.get("AF")
            if v is not None:
                af = float(v[0]) if isinstance(v, tuple) else float(v)
            d = smp.get("DP")
            if d is not None:
                dp = float(d)
        calls.append({
            "key": "%s:%d:%s>%s" % (rec.chrom, rec.pos, rec.ref, alt),
            "chrom": rec.chrom, "pos": rec.pos,
            "ref": rec.ref, "alt": alt,
            "PASS": passed,
            "filters": ",".join(filt) if filt else "PASS",
            "obs_vaf": af, "obs_depth": dp,
        })
vcf.close()

called = pd.DataFrame(calls)
print("      total SNV records: %d" % len(called), flush=True)
if len(called) == 0:
    raise SystemExit("no SNV records in the VCF - check the Mutect2 run")
print("      PASS records     : %d" % int(called.PASS.sum()), flush=True)

print("[3/6] matching truth against calls ...", flush=True)
pass_map = called[called.PASS].set_index("key")
any_map = called.set_index("key")
pass_map = pass_map[~pass_map.index.duplicated()]
any_map = any_map[~any_map.index.duplicated()]

truth["detected_pass"] = truth.key.isin(pass_map.index)
truth["detected_any"] = truth.key.isin(any_map.index)
truth["obs_vaf"] = truth.key.map(any_map.obs_vaf)
truth["obs_depth"] = truth.key.map(any_map.obs_depth)
truth["filters"] = truth.key.map(any_map.filters).fillna("NOT_CALLED")

n = len(truth)
n_pass = int(truth.detected_pass.sum())
n_any = int(truth.detected_any.sum())
print("      recovered (PASS)      : %d / %d" % (n_pass, n), flush=True)
print("      recovered (any filter): %d / %d" % (n_any, n), flush=True)

print("[4/6] recall by VAF bin ...", flush=True)
truth["vaf_bin"] = pd.cut(truth.VAF, bins=VAF_BINS, labels=VAF_LABELS, right=False)
curve = truth.groupby("vaf_bin", observed=False).agg(
    injected=("key", "size"),
    found_pass=("detected_pass", "sum"),
    found_any=("detected_any", "sum"),
    mean_target_vaf=("VAF", "mean"),
    mean_obs_vaf=("obs_vaf", "mean"),
).reset_index()
curve["recall_pass"] = (100 * curve.found_pass / curve.injected).round(1)
curve["recall_any"] = (100 * curve.found_any / curve.injected).round(1)

print("[5/6] false positives ...", flush=True)
truth_keys = set(truth.key)
fp = called[called.PASS & ~called.key.isin(truth_keys)].copy()
print("      PASS calls not in truth: %d" % len(fp), flush=True)

precision = 100.0 * n_pass / max(int(called.PASS.sum()), 1)

print("[6/6] writing ...", flush=True)
cols = ["Hugo_Symbol", "Chromosome_hg38", "Start_Position_hg38",
        "Reference_Allele", "Tumor_Seq_Allele2", "VAF",
        "detected_pass", "detected_any", "obs_vaf", "obs_depth", "filters"]
if "HGVSp_Short" in truth.columns:
    cols.insert(1, "HGVSp_Short")
truth[cols].sort_values("VAF").to_csv(
    "%s/truth_vs_calls.tsv" % OUT_DIR, sep="\t", index=False)
curve.to_csv("%s/sensitivity_curve.tsv" % OUT_DIR, sep="\t", index=False)
fp.to_csv("%s/false_positives.tsv" % OUT_DIR, sep="\t", index=False)

pd.set_option("display.width", 220)
print()
print("========== MUTECT2 vs TRUTH ==========")
print("injected SNVs        : %d" % n)
print("recovered, PASS      : %d  (recall %.1f%%)" % (n_pass, 100.0 * n_pass / n))
print("recovered, any filter: %d  (recall %.1f%%)" % (n_any, 100.0 * n_any / n))
print("PASS calls total     : %d" % int(called.PASS.sum()))
print("false positives      : %d  (precision %.1f%%)" % (len(fp), precision))
print()
print("--- SENSITIVITY CURVE BY VAF ---")
print(curve[["vaf_bin", "injected", "found_pass", "recall_pass",
             "found_any", "recall_any", "mean_target_vaf", "mean_obs_vaf"]]
      .to_string(index=False, float_format=lambda x: "%.3f" % x))
print()

miss = truth[~truth.detected_any]
if len(miss):
    print("--- COMPLETELY MISSED (%d) ---" % len(miss))
    mc = ["Hugo_Symbol", "Chromosome_hg38", "Start_Position_hg38", "VAF"]
    print(miss.sort_values("VAF")[mc].head(20)
          .to_string(index=False, float_format=lambda x: "%.3f" % x))
    print()

fo = truth[truth.detected_any & ~truth.detected_pass]
if len(fo):
    print("--- CALLED BUT FILTERED OUT (%d) ---" % len(fo))
    print(fo.filters.value_counts().head(10).to_string())
    print()

det = truth[truth.detected_any & truth.obs_vaf.notna()]
if len(det) > 2:
    from scipy import stats
    r, p = stats.pearsonr(det.VAF, det.obs_vaf)
    bias = (det.obs_vaf - det.VAF).mean()
    print("--- VAF RECOVERY ---")
    print("  n = %d   Pearson r = %.3f   p = %.3g" % (len(det), r, p))
    print("  mean observed - target = %+.4f" % bias)
    print("  a negative bias means BAMSurgeon under-injected, or the caller")
    print("  under-estimates the fraction")
    print()

if len(fp):
    print("--- FALSE POSITIVES, top 10 by depth ---")
    print(fp.sort_values("obs_depth", ascending=False).head(10)[
        ["chrom", "pos", "ref", "alt", "obs_vaf", "obs_depth"]]
        .to_string(index=False, float_format=lambda x: "%.3f" % x))
    print()

print("files written to %s" % OUT_DIR)
