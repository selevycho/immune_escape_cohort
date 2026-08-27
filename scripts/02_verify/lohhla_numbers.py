#!/usr/bin/env python3
"""Every LOHHLA number the slides need, in thirty lines."""
import os, sys, glob
import pandas as pd, numpy as np

WS = sys.argv[1] if len(sys.argv) > 1 else os.environ["WS"]
RES = os.path.expanduser("~/immune_escape_project/results")
S9 = f"{WS}/simulation/step9_lohhla_panel/results"
n = lambda s: pd.to_numeric(s, errors="coerce")

r = pd.read_csv(f"{RES}/lohhla_report.tsv", sep="\t")
r["pval"] = n(r.pval)
got, sig = r.pval.notna(), r.pval < 0.05

print(f"MAIN RUN")
print(f"  loci attempted        {len(r)}")
print(f"  produced a p-value    {int(got.sum())}  ({100*got.mean():.0f}%)")
print(f"  significant p<0.05    {int(sig.sum())}")
print(f"  samples with result   {r[got]['sample'].nunique()} of {r['sample'].nunique()}")
if "sites" in r:
    s = n(r.sites).dropna(); s = s[s > 0]
    print(f"  distinguishing sites  median {s.median():.0f}, "
          f"range {s.min():.0f}-{s.max():.0f}")

d = pd.read_csv(f"{S9}/step9_all_loci_mc5.tsv", sep="\t")
d["pval"] = n(d.pval)
print(f"\nWHY LOCI FAIL   ({int(d.pval.isna().sum())} of {len(d)} attempts)")
for k, v in d[d.reason.notna()].reason.str.slice(0,40).value_counts().head(4).items():
    print(f"  {v:>5}  {k}")

print(f"\nBY COLLECTION")
print(f"  {'':<12}{'loci':>6}{'result':>8}{'p<0.05':>8}{'on target':>11}")
for c, g in d.groupby("collection"):
    sg = g[g.pval < 0.05]
    on = int((sg.role == "TARGET").sum()) if "role" in g else 0
    print(f"  {c:<12}{len(g):>6}{int(g.pval.notna().sum()):>8}{len(sg):>8}{on:>11}")

a = d[d.pval < 0.05]
print(f"\n  all significant: {len(a)}, on target {int((a.role=='TARGET').sum())}, "
      f"on controls {int((a.role!='TARGET').sum())}")

roc = pd.read_csv(f"{S9}/lohhla_roc_mc5.tsv", sep="\t")
print(f"\nSENSITIVITY   (false positives at every threshold: "
      f"{int(roc.fp.sum())})")
for _, x in roc.iterrows():
    print(f"  p<{x.threshold:<7.3f}{int(x.tp):>4} of {int(x.tp+x.fn):<4}"
          f"{x.sensitivity:>7.1f}%")

sw = pd.read_csv(f"{S9.replace('step9_lohhla_panel','step16_imgt_sweep')}/imgt_sweep.tsv", sep="\t")
u = n(sw.pval).dropna().unique()
print(f"\nIMGT SUBTYPE SWEEP")
print(f"  configurations tried  {len(sw)}")
print(f"  distinct p-values     {len(u)}  -> {', '.join(f'{x:.2g}' for x in u)}")
print(f"  the subtype count changes nothing")

for c in (f"{WS}/simulation/lohhla_allelic/loh_design.tsv",
          f"{RES}/loh_design.tsv"):
    if os.path.exists(c):
        g = pd.read_csv(c, sep="\t")
        b = n(g.get("depth_A_before", g.get("A_before")))
        af = n(g.get("depth_A_after", g.get("A_after")))
        print(f"\nWHAT WAS SIMULATED")
        print(f"  samples thinned       {len(g)}")
        print(f"  target locus          " +
              ", ".join(f"{k} {v}" for k, v in g.target_locus.value_counts().items()))
        print(f"  depth before / after  {b.median():.0f}x -> {af.median():.0f}x")
        break
