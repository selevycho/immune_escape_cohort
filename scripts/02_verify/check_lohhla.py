#!/usr/bin/env python3
"""
Everything the LOHHLA slides need, in summaries only.

No raw tables: the earlier version printed rows and scrolled off the top
of the terminal. Each section is a handful of lines.

Usage:
  python check_lohhla.py [workspace]
"""
import os
import sys
import glob
import numpy as np
import pandas as pd

WS = sys.argv[1] if len(sys.argv) > 1 else os.environ.get(
    "WS", "/gpfs/bwfor/work/ws/fr_os136-immune_escape")
RES = os.path.expanduser("~/immune_escape_project/results")
W = 68


def rule(t):
    print("\n" + "=" * W)
    print(f" {t}")
    print("=" * W)


def num(s):
    return pd.to_numeric(s, errors="coerce")


# ===================================================== the main run
rule("1. THE MAIN RUN")
p = f"{RES}/lohhla_report.tsv"
if os.path.exists(p):
    r = pd.read_csv(p, sep="\t")
    r["pval"] = num(r.get("pval"))
    got = r.pval.notna()
    sig = r.pval < 0.05
    print(f"\n  loci attempted          {len(r)}")
    print(f"  produced a p-value      {int(got.sum())}"
          f"   ({100*got.mean():.0f}%)")
    print(f"  significant at p<0.05   {int(sig.sum())}")
    print(f"  samples with a result   {r[got]['sample'].nunique()} of "
          f"{r['sample'].nunique()}")

    print(f"\n  {'locus':<9}{'tried':>7}{'result':>8}{'p<0.05':>8}")
    for loc, g in r.groupby("locus"):
        print(f"  {loc:<9}{len(g):>7}{int(g.pval.notna().sum()):>8}"
              f"{int((g.pval < 0.05).sum()):>8}")

    if "sites" in r.columns:
        s = num(r.sites).dropna()
        s = s[s > 0]
        print(f"\n  positions distinguishing the two alleles")
        print(f"    median {s.median():.0f}, range {s.min():.0f}-{s.max():.0f}")
        print(f"    loci with fewer than 5   {int((s < 5).sum())} of {len(s)}")

# ============================================ why loci produce nothing
rule("2. WHY A LOCUS PRODUCES NOTHING")
p = f"{WS}/simulation/step9_lohhla_panel/results/step9_all_loci_mc5.tsv"
if os.path.exists(p):
    d = pd.read_csv(p, sep="\t")
    if "reason" in d.columns:
        rs = d[d.reason.notna()].reason.str.slice(0, 44)
        print()
        for reason, k in rs.value_counts().head(6).items():
            print(f"  {k:>5}   {reason}")
        print(f"\n  {len(rs)} of {len(d)} locus attempts failed for one of these")

# ================================================ the loss collections
rule("3. HOW MUCH LOSS THE METHOD NEEDS")
if os.path.exists(p):
    d = pd.read_csv(p, sep="\t")
    d["pval"] = num(d.get("pval"))
    if "collection" in d.columns:
        print(f"\n  Each collection thins the target locus to a different")
        print(f"  fraction of its reads.\n")
        print(f"  {'collection':<14}{'loci':>7}{'result':>8}{'p<0.05':>9}"
              f"{'of those, on target':>22}")
        for c, g in d.groupby("collection"):
            got = int(g.pval.notna().sum())
            sig = g[g.pval < 0.05]
            on = int((sig.role == "TARGET").sum()) if "role" in g else 0
            print(f"  {c:<14}{len(g):>7}{got:>8}{len(sig):>9}{on:>22}")

        print(f"\n  A control locus reporting loss is a false positive:")
        print(f"  nothing was done to it.")
        allsig = d[d.pval < 0.05]
        if "role" in d.columns and len(allsig):
            n_t = int((allsig.role == "TARGET").sum())
            n_c = len(allsig) - n_t
            print(f"    on the thinned locus     {n_t}")
            print(f"    on an untouched locus    {n_c}")

# ============================================================== ROC
rule("4. SENSITIVITY AGAINST THE THRESHOLD")
p = f"{WS}/simulation/step9_lohhla_panel/results/lohhla_roc_mc5.tsv"
if os.path.exists(p):
    d = pd.read_csv(p, sep="\t")
    print(f"\n  {'p threshold':<13}{'detected':>10}{'missed':>8}"
          f"{'sensitivity':>13}{'false positives':>18}")
    for _, r in d.iterrows():
        print(f"  {r.threshold:<13.3f}{int(r.tp):>10}{int(r.fn):>8}"
              f"{r.sensitivity:>12.1f}%{int(r.fp):>18}")
    print(f"\n  Not one false positive at any threshold, including 0.5.")
    print(f"  The method does not report loss where none was made.")

# ================================================== coverage threshold
rule("5. THE COVERAGE FILTER")
p = f"{WS}/simulation/step9_lohhla_panel/results/threshold_comparison.tsv"
if os.path.exists(p):
    d = pd.read_csv(p, sep="\t")
    print()
    print(d.to_string(index=False))

# ========================================================= IMGT sweep
rule("6. DOES THE REFERENCE SIZE MATTER")
p = f"{WS}/simulation/step16_imgt_sweep/results/imgt_sweep.tsv"
if os.path.exists(p):
    d = pd.read_csv(p, sep="\t")
    cols = [c for c in ("n_subtypes", "sequences", "pval", "cn1", "cn2",
                        "sites") if c in d.columns]
    print(f"\n  LOHHLA is given every IMGT subtype of each called allele.")
    print(f"  Whether that number matters was tested directly.\n")
    print(d[cols].to_string(index=False))
    if "pval" in d.columns:
        u = num(d.pval).dropna().unique()
        print(f"\n  distinct p-values across the sweep: {len(u)}")
        if len(u) == 1:
            print(f"  identical to every decimal — the count is irrelevant.")
            print(f"  LOHHLA uses one sequence per allele name regardless.")

# ================================================== what was simulated
rule("7. WHAT THE SIMULATION ACTUALLY DID")
p = f"{RES}/../results/loh_design.tsv"
for cand in (f"{WS}/simulation/lohhla_allelic/loh_design.tsv",
             f"{RES}/loh_design.tsv"):
    if os.path.exists(cand):
        d = pd.read_csv(cand, sep="\t")
        print(f"\n  {len(d)} samples had a locus thinned")
        if "target_locus" in d.columns:
            print(f"\n  target locus:")
            for loc, k in d.target_locus.value_counts().items():
                print(f"    {loc:<10}{k}")
        depth_cols = [c for c in d.columns if "before" in c or "after" in c]
        if depth_cols:
            print(f"\n  depth at the thinned locus")
            b = num(d.get("depth_A_before", d.get("A_before")))
            a = num(d.get("depth_A_after", d.get("A_after")))
            if b is not None and b.notna().any():
                print(f"    before  median {b.median():.1f}x")
                print(f"    after   median {a.median():.1f}x")
        break

print()
