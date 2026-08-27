#!/usr/bin/env python3
"""
What exists for the indels, before running anything against it.

Step 8 wrote its output into a separate BAM per sample and Mutect2 was
never pointed at it, so the pipeline has been shown to place indels and
not to find them. Closing that gap needs a calling run, and a calling run
needs to know what is actually on disk: which samples have an indel BAM,
whether it is indexed, whether a matched normal exists, and whether the
truth table lines up with the files.

Eight indels are also unaccounted for — the cohort truth sets carry 83
while indel_truth_all.tsv carries 75 — so this reports where they went
rather than leaving the discrepancy to be noticed later.

Nothing is modified. The output is a plan.

Usage:
  python survey_indel_files.py [workspace]
"""
import os
import sys
import subprocess
import pandas as pd

WS = sys.argv[1] if len(sys.argv) > 1 else os.environ.get(
    "WS", "/gpfs/bwfor/work/ws/fr_os136-immune_escape")

COHORT = f"{WS}/simulation/cohort"
INDELS = f"{WS}/simulation/indels"
TRUTH = f"{INDELS}/indel_truth_all.tsv"
MANIFEST = f"{COHORT}/manifest.tsv"

SAMTOOLS = None
for cand in ["samtools",
             os.path.expanduser("~/miniconda3/envs/bio_work/bin/samtools"),
             "/home/fr/fr_fr/fr_os136/miniconda3/envs/bio_work/bin/samtools"]:
    if subprocess.run(f"{cand} --version", shell=True,
                      capture_output=True).returncode == 0:
        SAMTOOLS = cand
        break

man = pd.read_csv(MANIFEST, sep="\t")
truth = pd.read_csv(TRUTH, sep="\t") if os.path.exists(TRUTH) else pd.DataFrame()

W = 74
print("=" * W)
print(" 1. FILES ON DISK")
print("=" * W)

rows = []
for _, m in man.iterrows():
    sid = m.sample_id
    d = f"{INDELS}/{sid}"
    bam = f"{d}/{sid}_tumor_snv_indel.bam"
    rows.append({
        "sample": sid, "cohort": m.cohort,
        "dir": os.path.isdir(d),
        "bam": os.path.exists(bam),
        "bai": os.path.exists(bam + ".bai"),
        "size_MB": round(os.path.getsize(bam) / 1e6, 1) if os.path.exists(bam) else 0,
        "log": os.path.exists(f"{d}/addindel.log"),
        "vcf": os.path.exists(f"{d}/raw.addindel.indels.vcf"),
        "normal": os.path.exists(f"{COHORT}/{sid}/{sid}_normal.bam"),
        "n_truth": int((truth["sample"] == sid).sum()) if len(truth) else 0,
    })

f = pd.DataFrame(rows)
pd.set_option("display.width", 200)
print(f"\n  samples in the manifest        {len(f)}")
print(f"  with an indel directory        {int(f.dir.sum())}")
print(f"  with an indel BAM              {int(f.bam.sum())}")
print(f"  indexed                        {int(f.bai.sum())}")
print(f"  with a matched normal          {int(f.normal.sum())}")
print(f"  total size                     {f.size_MB.sum()/1024:.1f} GB")

ready = f[f.bam & f.bai & f.normal & (f.n_truth > 0)]
print(f"\n  ready to call                  {len(ready)}")

problems = f[(f.n_truth > 0) & ~(f.bam & f.bai & f.normal)]
if len(problems):
    print(f"\n  have indels but cannot be called:")
    print(problems[["sample", "bam", "bai", "normal", "n_truth"]]
          .to_string(index=False))

none_ = f[f.n_truth == 0]
if len(none_):
    print(f"\n  no indel inside the panel      {len(none_)} samples")
    print(f"    {' '.join(none_['sample'])}")

print()
print("=" * W)
print(" 2. THE EIGHT MISSING INDELS")
print("=" * W)

cohort_ind = []
for _, m in man.iterrows():
    p = f"{COHORT}/{m.sample_id}/truth_set.tsv"
    if not os.path.exists(p):
        continue
    t = pd.read_csv(p, sep="\t")
    if "Variant_Type" not in t.columns:
        continue
    g = t[t.Variant_Type != "SNP"].copy()
    g["sample"] = m.sample_id
    cohort_ind.append(g)

if cohort_ind:
    C = pd.concat(cohort_ind, ignore_index=True)
    print(f"\n  indels in the cohort truth sets   {len(C)}")
    print(f"  indels in indel_truth_all.tsv     {len(truth)}")
    print(f"  difference                        {len(C) - len(truth)}")

    if len(truth):
        pos_col = "pos" if "pos" in truth.columns else "Start_Position_hg38"
        chr_col = "chrom" if "chrom" in truth.columns else "Chromosome_hg38"
        have = set(zip(truth["sample"], truth[chr_col], truth[pos_col]))
        C["key"] = list(zip(C["sample"], C.Chromosome_hg38,
                            C.Start_Position_hg38))
        gone = C[~C.key.isin(have)]
        if len(gone):
            print(f"\n  present in a cohort truth set but not in the indel run:\n")
            cols = [c for c in ["sample", "Hugo_Symbol", "Chromosome_hg38",
                                "Start_Position_hg38", "Variant_Type",
                                "Variant_Classification", "VAF"]
                    if c in gone.columns]
            print(gone[cols].to_string(index=False))
            print(f"\n  by type:")
            for t_, k in gone.Variant_Type.value_counts().items():
                print(f"    {t_:<12}{k}")
            if "VAF" in gone.columns:
                print(f"\n  their requested fractions: "
                      f"{', '.join(f'{v:.3f}' for v in gone.VAF)}")
            print(f"\n  Whether these were dropped before injection or simply")
            print(f"  never written to the collected table is the question;")
            print(f"  the addindel logs below will say which.")

            for sid in gone["sample"].unique()[:3]:
                log = f"{INDELS}/{sid}/addindel.log"
                if os.path.exists(log):
                    txt = open(log, errors="ignore").read()
                    hits = [l for l in txt.splitlines()
                            if any(w in l.lower() for w in
                                   ("skip", "fail", "warn", "error", "no reads"))]
                    if hits:
                        print(f"\n  {sid} addindel.log, relevant lines:")
                        for h in hits[:5]:
                            print(f"    {h.strip()[:100]}")

print()
print("=" * W)
print(" 3. WHAT A CALLING RUN NEEDS")
print("=" * W)
print(f"\n  tumour   {INDELS}/<sample>/<sample>_tumor_snv_indel.bam")
print(f"  normal   {COHORT}/<sample>/<sample>_normal.bam")
print(f"  panel    {WS}/simulation/panel/panel.bed")
print(f"  truth    {TRUTH}")
print(f"\n  The tumour BAM already carries the substitutions, so a run")
print(f"  against it recovers both kinds at once and the comparison has")
print(f"  to separate them by variant type rather than by file.")

if SAMTOOLS and len(ready):
    sid = ready.iloc[0]["sample"]
    bam = f"{INDELS}/{sid}/{sid}_tumor_snv_indel.bam"
    print(f"\n  read group in {sid}:")
    rg = subprocess.run(f"{SAMTOOLS} view -H {bam} | grep '^@RG' | head -1",
                        shell=True, capture_output=True, text=True).stdout.strip()
    print(f"    {rg[:140] if rg else 'none — Mutect2 will refuse without one'}")

out = os.path.expanduser("~/immune_escape_project/results/indel_file_survey.tsv")
f.to_csv(out, sep="\t", index=False)
print(f"\nwritten to {out}")
