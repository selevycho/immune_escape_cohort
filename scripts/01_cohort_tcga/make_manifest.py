#!/usr/bin/env python3
"""
Build the manifest that drives the whole simulated cohort.

Every downstream script reads one line of this file and knows everything it
needs: which TCGA patient supplies the mutations, which 1000 Genomes sample
supplies the reads, and where the CRAM lives.

Patient selection is stratified by panel mutation count, so the cohort spans
the whole range rather than clustering at the median. That matters because
the sensitivity curve is the main deliverable: a cohort of uniformly
hypermutated cases would say nothing about how the caller behaves when there
is little to find.

Backbone genomes are reused across patients. Thirty distinct CRAMs would mean
thirty slow downloads for no scientific gain - the reads are only a substrate
for the injected mutations. Sixteen backbones spread across five
superpopulations keep the HLA diversity that LOHHLA needs while cutting the
download cost roughly in half.

Usage:
  python make_manifest.py <liftover_dir> <panel_bed> <backbone_tsv> \
                          <sequence_index> <out_tsv> [n_brca] [n_ov] [min_mut]
"""
import sys, os
import numpy as np
import pandas as pd

LIFT_DIR = sys.argv[1]
PANEL = sys.argv[2]
BACKBONE_TSV = sys.argv[3]
SEQ_INDEX = sys.argv[4]
OUT_TSV = sys.argv[5]
N_BRCA = int(sys.argv[6]) if len(sys.argv) > 6 else 15
N_OV = int(sys.argv[7]) if len(sys.argv) > 7 else 15
MIN_MUT = int(sys.argv[8]) if len(sys.argv) > 8 else 15

MAX_MUT = 400          # beyond this BAMSurgeon runtime grows without benefit
MIN_VAF = 0.05
SEED = 42

NONSYN = {"Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
          "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Splice_Site"}

print("[1/5] loading the panel ...", flush=True)
panel = pd.read_csv(PANEL, sep="\t", header=None,
                    names=["chrom", "start", "end", "name"])
by_chrom = {c: g[["start", "end"]].to_numpy()
            for c, g in panel.groupby("chrom")}
total_mb = (panel.end - panel.start).sum() / 1e6
print("      %d intervals, %.1f Mb" % (len(panel), total_mb), flush=True)


def count_in_panel(sub):
    n = 0
    for c, pos in zip(sub.Chromosome_hg38.values,
                      sub.Start_Position_hg38.values):
        iv = by_chrom.get(c)
        if iv is None:
            continue
        p = int(pos) - 1
        if ((iv[:, 0] <= p) & (p < iv[:, 1])).any():
            n += 1
    return n


print("[2/5] counting panel mutations per patient ...", flush=True)
counts = {}
for cohort in ["brca", "ov"]:
    maf = pd.read_csv("%s/%s.hg38.maf.tsv" % (LIFT_DIR, cohort),
                      sep="\t", low_memory=False)
    maf = maf[maf.Variant_Classification.isin(NONSYN)]
    maf = maf[maf.Variant_Type == "SNP"]
    den = maf.t_ref_count.fillna(0) + maf.t_alt_count.fillna(0)
    maf = maf[(den > 0) & (maf.t_alt_count / den >= MIN_VAF)]

    rows = []
    for bc, sub in maf.groupby("Tumor_Sample_Barcode"):
        n = count_in_panel(sub)
        if n >= MIN_MUT:
            rows.append({"barcode": bc, "n_panel_mut": n,
                         "n_total_mut": len(sub)})
    c = pd.DataFrame(rows).sort_values("n_panel_mut", ascending=False)
    counts[cohort] = c
    print("      %s: %d eligible patients (>=%d panel mutations)"
          % (cohort.upper(), len(c), MIN_MUT), flush=True)
    if len(c):
        print("        range %d - %d, median %d"
              % (c.n_panel_mut.min(), c.n_panel_mut.max(),
                 c.n_panel_mut.median()), flush=True)

print("[3/5] stratified selection ...", flush=True)


def stratified(c, n):
    c = c[(c.n_panel_mut >= MIN_MUT) & (c.n_panel_mut <= MAX_MUT)].copy()
    if len(c) <= n:
        return c
    try:
        c["stratum"] = pd.qcut(c.n_panel_mut, 3,
                               labels=["low", "mid", "high"],
                               duplicates="drop")
    except ValueError:
        c["stratum"] = "all"
    per = max(1, n // c.stratum.nunique())
    out = [g.sample(min(per, len(g)), random_state=SEED)
           for _, g in c.groupby("stratum", observed=True)]
    sel = pd.concat(out)
    if len(sel) < n:
        rest = c[~c.barcode.isin(sel.barcode)]
        if len(rest):
            sel = pd.concat([sel, rest.sample(min(n - len(sel), len(rest)),
                                              random_state=SEED)])
    return sel.sort_values("n_panel_mut", ascending=False).head(n)


sel = {}
for cohort, n in [("brca", N_BRCA), ("ov", N_OV)]:
    s = stratified(counts[cohort], n)
    sel[cohort] = s
    print("      %s selected: %d  (%d - %d mutations, median %d)"
          % (cohort.upper(), len(s), s.n_panel_mut.min(),
             s.n_panel_mut.max(), s.n_panel_mut.median()), flush=True)

print("[4/5] resolving backbone CRAM URLs ...", flush=True)
bb = pd.read_csv(BACKBONE_TSV, sep="\t")
print("      backbones listed: %d" % len(bb), flush=True)

idx_rows = []
with open(SEQ_INDEX) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 11 or not f[0].endswith(".cram"):
            continue
        idx_rows.append({"cram_url": f[0].replace("ftp://", "https://", 1),
                         "err_id": f[2], "sample": f[9]})
idx = pd.DataFrame(idx_rows).drop_duplicates("sample")
print("      CRAM entries in the index: %d" % len(idx), flush=True)

bb = bb.merge(idx, on="sample", how="inner")
print("      backbones resolved: %d" % len(bb), flush=True)
for _, r in bb.iterrows():
    print("        %-9s %-4s %-4s %s"
          % (r["sample"], r.population, r.superpopulation, r.err_id),
          flush=True)

if len(bb) == 0:
    raise SystemExit("no backbone resolved - check the sample names")

# shuffle so cohorts do not end up systematically paired with
# particular ancestries - that would confound cohort with backbone
bb = bb.sample(frac=1, random_state=SEED).reset_index(drop=True)

print("[5/5] writing the manifest ...", flush=True)
rows = []
k = 0
for cohort, prefix in [("brca", "B"), ("ov", "O")]:
    for i, (_, p) in enumerate(sel[cohort].iterrows(), start=1):
        b = bb.iloc[k % len(bb)]
        k += 1
        rows.append({
            "sample_id": "%s%03d" % (prefix, i),
            "cohort": cohort,
            "tcga_barcode": p.barcode,
            "n_panel_mut": int(p.n_panel_mut),
            "n_total_mut": int(p.n_total_mut),
            "backbone": b["sample"],
            "population": b.population,
            "superpopulation": b.superpopulation,
            "err_id": b.err_id,
            "cram_url": b.cram_url,
            "status": "pending",
        })

man = pd.DataFrame(rows)
man.to_csv(OUT_TSV, sep="\t", index=False)

pd.set_option("display.width", 240)
print()
print("=================== MANIFEST ===================")
print("samples            : %d  (%d BRCA, %d OV)"
      % (len(man), (man.cohort == "brca").sum(), (man.cohort == "ov").sum()))
print("distinct backbones : %d" % man.backbone.nunique())
print("panel mutations    : %d - %d, median %d, total %d"
      % (man.n_panel_mut.min(), man.n_panel_mut.max(),
         man.n_panel_mut.median(), man.n_panel_mut.sum()))
print()
print("backbone reuse:")
print(man.backbone.value_counts().to_string())
print()
print(man[["sample_id", "cohort", "tcga_barcode", "n_panel_mut",
           "backbone", "superpopulation", "err_id"]].to_string(index=False))
print()
print("written to %s" % OUT_TSV)
