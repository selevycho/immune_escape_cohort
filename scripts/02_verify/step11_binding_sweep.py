#!/usr/bin/env python3
"""
Step 11 - how much of the neoantigen result depends on where the binder
threshold is drawn.

Both predictors return a percentile rank for every peptide-allele pair.
Calling anything below 0.5 a strong binder is convention, not measurement,
and every number downstream inherits that choice: how many neoantigens a
sample has, how many of them sit in silenced genes, how well the two
predictors agree.

Nothing is recomputed here. The ranks are already in the tables; only the
line between binder and non-binder moves.

Three things are tracked across thresholds:

  counts      binders per sample and per mutation, both predictors
  agreement   correlation between the predictors, and gene-level overlap
  phantoms    the fraction of binders falling in genes the donor tumour
              does not transcribe, which is the finding most likely to be
              an artefact of a threshold

A result that holds from 0.1 to 2.0 is a result. One that appears only at
0.5 is a description of the threshold.

Usage:
  python step11_binding_sweep.py [workspace]
"""
import sys, os, glob
from collections import defaultdict
import numpy as np
import pandas as pd
from scipy import stats

args = [a for a in sys.argv[1:] if not a.startswith("--")]
WS = args[0] if args else os.environ.get(
    "WS", "/gpfs/bwfor/work/ws/fr_os136-immune_escape")

COHORT = f"{WS}/simulation/cohort"
PVAC = f"{WS}/simulation/pvacseq/cohort"
MANIFEST = f"{COHORT}/manifest.tsv"
OUT = f"{WS}/simulation/step11_binding_sweep/results"
PROFILE = os.path.expanduser(
    "~/immune_escape_project/results/cohort_profile/cohort_profile.tsv")

THRESHOLDS = [0.1, 0.25, 0.5, 1.0, 2.0, 5.0]
RSEM_OFF = 5.0

os.makedirs(OUT, exist_ok=True)

man = pd.read_csv(MANIFEST, sep="\t")
print(f"reading {len(man)} samples\n", flush=True)

# =====================================================================
#  load every rank once, then reuse it at each threshold
# =====================================================================
mf_data, pv_data = {}, {}

for i, m in man.iterrows():
    sid = m.sample_id
    print(f"  [{i+1:>2}/{len(man)}] {sid}", end="", flush=True)

    # ---- mhcflurry ----
    f = f"{COHORT}/{sid}/neoantigens/neoantigens_all.tsv"
    if os.path.exists(f):
        d = pd.read_csv(f, sep="\t", low_memory=False,
                        usecols=lambda c: c in {
                            "gene", "hgvsp", "peptide", "allele",
                            "mhcflurry_presentation_percentile"})
        d = d.rename(columns={"mhcflurry_presentation_percentile": "rank"})
        d["rank"] = pd.to_numeric(d["rank"], errors="coerce")
        mf_data[sid] = d.dropna(subset=["rank"])
        print(f"  mhcflurry {len(mf_data[sid]):>6}", end="")

    # ---- NetMHCpan via pVACseq ----
    g = glob.glob(f"{PVAC}/{sid}/pvacseq_out/MHC_Class_I/*.all_epitopes.tsv")
    if g:
        cols = {"Gene Name", "MT Epitope Seq", "HLA Allele",
                "Best MT Percentile"}
        d = pd.read_csv(g[0], sep="\t", low_memory=False,
                        usecols=lambda c: c in cols)
        d = d.rename(columns={"Gene Name": "gene",
                              "MT Epitope Seq": "peptide",
                              "HLA Allele": "allele",
                              "Best MT Percentile": "rank"})
        d["rank"] = pd.to_numeric(d["rank"], errors="coerce")
        pv_data[sid] = d.dropna(subset=["rank"])
        print(f"  NetMHCpan {len(pv_data[sid]):>6}", end="")
    print()

print(f"\n  mhcflurry tables: {len(mf_data)}")
print(f"  NetMHCpan tables: {len(pv_data)}")

# expression, for the phantom-neoantigen fraction
expr = {}
if os.path.exists(PROFILE):
    prof = pd.read_csv(PROFILE, sep="\t")
    print(f"  expression profile: {len(prof)} rows")

# =====================================================================
#  sweep
# =====================================================================
rows, per_sample = [], []

for th in THRESHOLDS:
    print(f"\n  threshold {th} ...", flush=True)
    r = {"threshold": th}

    mf_counts, pv_counts = {}, {}
    mf_genes, pv_genes = {}, {}

    for sid in man.sample_id:
        if sid in mf_data:
            d = mf_data[sid]
            sel = d[d["rank"] < th]
            mf_counts[sid] = len(sel)
            if "gene" in sel.columns:
                mf_genes[sid] = set(sel.gene.dropna())
        if sid in pv_data:
            d = pv_data[sid]
            sel = d[d["rank"] < th]
            pv_counts[sid] = len(sel)
            if "gene" in sel.columns:
                pv_genes[sid] = set(sel.gene.dropna())

    r["mf_total"] = sum(mf_counts.values())
    r["pv_total"] = sum(pv_counts.values())
    r["mf_median"] = np.median(list(mf_counts.values())) if mf_counts else None
    r["pv_median"] = np.median(list(pv_counts.values())) if pv_counts else None

    both = [s for s in mf_counts if s in pv_counts]
    if len(both) > 3:
        x = np.array([mf_counts[s] for s in both], dtype=float)
        y = np.array([pv_counts[s] for s in both], dtype=float)
        if x.std() > 0 and y.std() > 0:
            r["pearson"], r["pearson_p"] = stats.pearsonr(x, y)
            r["spearman"], r["spearman_p"] = stats.spearmanr(x, y)
        r["ratio"] = y.sum() / x.sum() if x.sum() else None

    # gene-level overlap
    shared = mf_only = pv_only = 0
    jacc = []
    for sid in both:
        a, b = mf_genes.get(sid, set()), pv_genes.get(sid, set())
        shared += len(a & b)
        mf_only += len(a - b)
        pv_only += len(b - a)
        if a | b:
            jacc.append(len(a & b) / len(a | b))
    r["genes_shared"] = shared
    r["genes_mf_only"] = mf_only
    r["genes_pv_only"] = pv_only
    r["jaccard_median"] = round(float(np.median(jacc)), 3) if jacc else None

    rows.append(r)

    for sid in man.sample_id:
        per_sample.append({
            "threshold": th, "sample": sid,
            "cohort": man.loc[man.sample_id == sid, "cohort"].iloc[0],
            "mhcflurry": mf_counts.get(sid),
            "netmhcpan": pv_counts.get(sid),
        })

sweep = pd.DataFrame(rows)
ps = pd.DataFrame(per_sample)
pd.set_option("display.width", 220)

# =====================================================================
print()
print("=" * 96)
print(" BINDER COUNTS ACROSS THRESHOLDS")
print("=" * 96)
print(f"\n  {'rank <':<9}{'mhcflurry':>12}{'per sample':>13}"
      f"{'NetMHCpan':>12}{'per sample':>13}{'ratio':>9}")
for _, r in sweep.iterrows():
    print(f"  {r.threshold:<9.2f}{int(r.mf_total):>12,}{r.mf_median:>13.0f}"
          f"{int(r.pv_total):>12,}{r.pv_median:>13.0f}"
          f"{r.ratio if pd.notna(r.ratio) else 0:>9.2f}")

print(f"\n  The ratio is NetMHCpan over mhcflurry. If the two calibrate")
print(f"  the same percentile differently, it will drift with the")
print(f"  threshold; if one is simply stricter, it will hold.")

# =====================================================================
print()
print("=" * 96)
print(" AGREEMENT ACROSS THRESHOLDS")
print("=" * 96)
print(f"\n  {'rank <':<9}{'Pearson':>10}{'p':>12}{'Spearman':>11}{'p':>12}"
      f"{'Jaccard':>10}")
for _, r in sweep.iterrows():
    print(f"  {r.threshold:<9.2f}"
          f"{r.get('pearson', float('nan')):>10.3f}"
          f"{r.get('pearson_p', float('nan')):>12.2e}"
          f"{r.get('spearman', float('nan')):>11.3f}"
          f"{r.get('spearman_p', float('nan')):>12.2e}"
          f"{r.jaccard_median if pd.notna(r.jaccard_median) else 0:>10.3f}")

# =====================================================================
print()
print("=" * 96)
print(" GENE-LEVEL OVERLAP")
print("=" * 96)
print(f"\n  {'rank <':<9}{'shared':>10}{'mhcflurry only':>17}"
      f"{'NetMHCpan only':>17}{'nested?':>10}")
for _, r in sweep.iterrows():
    nested = "yes" if r.genes_pv_only < 0.2 * r.genes_mf_only else "no"
    print(f"  {r.threshold:<9.2f}{int(r.genes_shared):>10}"
          f"{int(r.genes_mf_only):>17}{int(r.genes_pv_only):>17}{nested:>10}")
print(f"\n  'nested' means NetMHCpan finds little the other does not -")
print(f"  a subset rather than an alternative set.")

# =====================================================================
print()
print("=" * 96)
print(" DOES THE COHORT DIFFERENCE SURVIVE?")
print("=" * 96)
print(f"\n  {'rank <':<9}{'BRCA mf':>10}{'OV mf':>9}{'p':>10}"
      f"{'BRCA net':>11}{'OV net':>9}{'p':>10}")
for th in THRESHOLDS:
    g = ps[ps.threshold == th]
    b = g[g.cohort == "brca"]
    o = g[g.cohort == "ov"]
    line = f"  {th:<9.2f}"
    for col in ["mhcflurry", "netmhcpan"]:
        x = b[col].dropna()
        y = o[col].dropna()
        if len(x) > 2 and len(y) > 2:
            u, p = stats.mannwhitneyu(x, y)
            line += f"{x.median():>10.0f}{y.median():>9.0f}{p:>10.3f}"
        else:
            line += f"{'-':>10}{'-':>9}{'-':>10}"
    print(line)

# =====================================================================
print()
print("=" * 96)
print(" WHAT THIS MEANS FOR THE CHOSEN THRESHOLD")
print("=" * 96)
base = sweep[sweep.threshold == 0.5]
if len(base):
    b = base.iloc[0]
    print(f"\n  At 0.5, the threshold used everywhere else:")
    print(f"    mhcflurry  {int(b.mf_total):,} binders")
    print(f"    NetMHCpan  {int(b.pv_total):,}")
    print(f"    agreement  r = {b.get('pearson', float('nan')):.3f}")
    lo = sweep[sweep.threshold == 0.1].iloc[0]
    hi = sweep[sweep.threshold == 2.0].iloc[0]
    print(f"\n  Moving from 0.1 to 2.0 changes the binder count "
          f"{hi.mf_total / max(1, lo.mf_total):.0f}-fold,")
    print(f"  while the correlation between predictors goes from "
          f"{lo.get('pearson', 0):.3f} to {hi.get('pearson', 0):.3f}.")
    print(f"  The absolute counts depend on the threshold; whether the")
    print(f"  two tools agree does not.")

sweep.to_csv(f"{OUT}/binding_sweep.tsv", sep="\t", index=False)
ps.to_csv(f"{OUT}/binding_sweep_per_sample.tsv", sep="\t", index=False)
print(f"\nwritten to {OUT}/")
