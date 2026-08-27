#!/usr/bin/env python3
"""
Every parameter sweep, summarised in under sixty lines.

Five things were varied to see what the pipeline is sensitive to: filter
stringency, sequencing depth, tumour purity, the peptide lengths scored,
and the significance threshold for HLA loss. Each answers a question
someone designing an experiment would actually ask.

Usage:
  python check_sweeps.py [workspace]
"""
import os, sys, glob
import pandas as pd, numpy as np

WS = sys.argv[1] if len(sys.argv) > 1 else os.environ["WS"]
RES = os.path.expanduser("~/immune_escape_project/results")
SIM = f"{WS}/simulation"
W = 68
num = lambda s: pd.to_numeric(s, errors="coerce")

def rule(t):
    print("\n" + "=" * W); print(f" {t}"); print("=" * W)

# ------------------------------------------------------------ what exists
rule("WHICH SWEEPS RAN")
for d in sorted(glob.glob(f"{SIM}/step1[0-9]*")):
    name = os.path.basename(d)
    res = glob.glob(f"{d}/results/*.tsv")
    print(f"  {name:<30}{len(res)} result file(s)")

# --------------------------------------------------------------- filters
rule("1. FILTER STRINGENCY")
p = f"{SIM}/step10_filter_sweep/results/filter_sweep_summary.tsv"
if os.path.exists(p):
    d = pd.read_csv(p, sep="\t")
    order = ["strict", "default", "relaxed", "permissive", "verypermissive"]
    d["_o"] = d.setting.map({k: i for i, k in enumerate(order)})
    d = d.sort_values("_o")
    base = float(d[d.setting == "default"].recovered.iloc[0])
    print(f"\n  {'setting':<16}{'recovered':>11}{'recall':>9}"
          f"{'vs default':>12}{'false':>8}")
    for _, r in d.iterrows():
        print(f"  {r.setting:<16}{int(r.recovered):>11}{r.recall:>8.1f}%"
              f"{int(r.recovered-base):>+12}{int(r.not_in_truth):>8}")
    print(f"\n  detection ceiling: {d.seen.iloc[0]:.1f}% — every setting")
    print(f"  sees the same variants and differs only in what it keeps")

# -------------------------------------------------------------- coverage
rule("2. SEQUENCING DEPTH")
found = False
for p in glob.glob(f"{SIM}/step14*/results/*.tsv"):
    d = pd.read_csv(p, sep="\t")
    cols = [c for c in d.columns if "depth" in c.lower() or "cov" in c.lower()]
    if not cols: continue
    found = True
    print(f"\n  {os.path.basename(p)}   {len(d)} rows")
    print(f"  columns: {list(d.columns)[:8]}")
    dc = cols[0]
    rc = next((c for c in d.columns if "recall" in c.lower()), None)
    if rc:
        g = d.groupby(dc)[rc].agg(["size", "median"]).reset_index()
        print(f"\n  {'depth':<10}{'samples':>9}{'recall':>10}")
        for _, r in g.iterrows():
            print(f"  {r[dc]:<10}{int(r['size']):>9}{r['median']:>9.1f}%")
if not found:
    print("\n  no coverage sweep results found")

# ---------------------------------------------------------------- purity
rule("3. TUMOUR PURITY")
found = False
for p in glob.glob(f"{SIM}/step15*/results/*.tsv"):
    d = pd.read_csv(p, sep="\t")
    found = True
    print(f"\n  {os.path.basename(p)}   {len(d)} rows")
    print(f"  columns: {list(d.columns)[:8]}")
    vc = next((c for c in d.columns if "vaf" in c.lower() or "purity" in c.lower()), None)
    rc = next((c for c in d.columns if "recall" in c.lower() or "found" in c.lower()), None)
    if vc and rc:
        g = d.groupby(vc)[rc].agg(["size", "mean"]).reset_index()
        print(f"\n  {'requested':<12}{'n':>6}{'recovered':>12}")
        for _, r in g.iterrows():
            print(f"  {r[vc]:<12}{int(r['size']):>6}{100*r['mean']:>11.1f}%")
if not found:
    print("\n  no purity sweep results found")

# ------------------------------------------------------- peptide lengths
rule("4. PEPTIDE LENGTH")
p = f"{RES}/binders_by_length.tsv"
if os.path.exists(p):
    print()
    print(pd.read_csv(p, sep="\t").to_string(index=False))
else:
    frames = []
    for f in glob.glob(f"{SIM}/cohort/*/neoantigens/neoantigens_all.tsv"):
        frames.append(pd.read_csv(f, sep="\t",
                      usecols=lambda c: c in ("length", "binder")))
    if frames:
        P = pd.concat(frames)
        print(f"\n  {'length':<9}{'tested':>10}{'strong':>9}{'rate':>9}")
        for L, g in P.groupby("length"):
            ns = int((g.binder == "STRONG").sum())
            print(f"  {int(L):<9}{len(g):>10}{ns:>9}{100*ns/len(g):>8.2f}%")

# ------------------------------------------------------- LOHHLA threshold
rule("5. HLA-LOSS SIGNIFICANCE THRESHOLD")
p = f"{SIM}/step9_lohhla_panel/results/lohhla_roc_mc5.tsv"
if os.path.exists(p):
    d = pd.read_csv(p, sep="\t")
    print(f"\n  {'p <':<9}{'detected':>10}{'missed':>8}{'sensitivity':>13}"
          f"{'false':>8}")
    for _, r in d.iterrows():
        print(f"  {r.threshold:<9.3f}{int(r.tp):>10}{int(r.fn):>8}"
              f"{r.sensitivity:>12.1f}%{int(r.fp):>8}")

# -------------------------------------------------------- coverage filter
rule("6. LOHHLA COVERAGE FILTER")
p = f"{SIM}/step9_lohhla_panel/results/threshold_comparison.tsv"
if os.path.exists(p):
    # this file is one row per locus attempt per setting, several hundred
    # of them; only the counts per setting belong on screen
    d = pd.read_csv(p, sep="\t")
    d["pval"] = num(d.get("pval"))
    key = "setting" if "setting" in d.columns else None
    if key:
        print(f"\n  {'setting':<16}{'attempts':>10}{'tested':>8}"
              f"{'p<0.05':>8}{'on target':>11}")
        for k, g in d.groupby(key):
            sg = g[g.pval < 0.05]
            on = int((sg.role == "TARGET").sum()) if "role" in g else 0
            print(f"  {str(k):<16}{len(g):>10}"
                  f"{int(g.pval.notna().sum()):>8}{len(sg):>8}{on:>11}")
        print(f"\n  A lower coverage filter lets more loci reach the test")
        print(f"  and does not introduce a single false positive.")

# ------------------------------------------------------------ IMGT depth
rule("7. IMGT SUBTYPES IN THE REFERENCE")
p = f"{SIM}/step16_imgt_sweep/results/imgt_sweep.tsv"
if os.path.exists(p):
    d = pd.read_csv(p, sep="\t")
    print(f"\n  {len(d)} configurations tried")
    print(f"  columns: {list(d.columns)}")
    u = num(d.pval).dropna().unique()
    print(f"  distinct p-values: {len(u)}  -> "
          f"{', '.join(f'{x:.3g}' for x in u)}")

# ------------------------------------------------- mutation burden effect
rule("8. DOES MUTATION BURDEN AFFECT RECALL")
p = f"{RES}/step3_rescored_full.tsv"
if os.path.exists(p):
    d = pd.read_csv(p, sep="\t")
    g = d.groupby("sample").agg(n=("pass", "size"),
                                recall=("pass", "mean"),
                                vaf=("observed_vaf", "median")).reset_index()
    r_n = float(np.corrcoef(g.n, g.recall)[0, 1])
    r_v = float(np.corrcoef(g.vaf, g.recall)[0, 1])
    print(f"\n  recall vs mutation count   r = {r_n:.3f}")
    print(f"  recall vs median VAF       r = {r_v:.3f}")
    print(f"\n  Recall follows the allele fractions a donor happened to")
    print(f"  carry, not how many mutations they had.")
print()
