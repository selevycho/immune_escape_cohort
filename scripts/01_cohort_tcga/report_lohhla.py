#!/usr/bin/env python3
"""
Report the LOHHLA results, sample by sample.

Reads every prediction that came back with content and lays it beside the
design: which locus was thinned, what LOHHLA said about it, and what it
said about the loci that were left alone.

Copy numbers are the quantity of interest. Under the simulation the
thinned allele should carry roughly a third of the reads of its partner,
so a wide gap between the two copy number estimates is the signal, and
PVal is the test of whether that gap is larger than the noise.

Samples where LOHHLA could not finish are listed separately with the
reason, since the reason - too few positions distinguishing the alleles -
is itself a result about panel data.

Usage:
  python report_lohhla.py [workspace]
"""
import sys, os, glob
import pandas as pd

WS = sys.argv[1] if len(sys.argv) > 1 else os.environ.get(
    "WS", "/gpfs/bwfor/work/ws/fr_os136-immune_escape")

LOH = f"{WS}/simulation/lohhla_allelic"
MANIFEST = f"{WS}/simulation/cohort/manifest.tsv"
DESIGN = f"{LOH}/loh_design.tsv"
ALPHA = 0.05

man = pd.read_csv(MANIFEST, sep="\t")
design = pd.read_csv(DESIGN, sep="\t") if os.path.exists(DESIGN) else pd.DataFrame()
target_of = dict(zip(design["sample"], design.target_locus)) if len(design) else {}
alleles_of = dict(zip(design["sample"], design.target_alleles)) if len(design) else {}

rows, failed = [], []

for _, m in man.iterrows():
    sid = m.sample_id
    got_any = False
    for locus in ["a", "b", "c"]:
        files = [f for f in glob.glob(f"{LOH}/{sid}/out_{locus}/*HLAlossPrediction*.txt")
                 if os.path.getsize(f) > 200
                 and "homozygous" not in f and "No_Suitable" not in f]
        if not files:
            continue
        try:
            t = pd.read_csv(files[0], sep="\t")
        except Exception:
            continue
        if not len(t):
            continue
        got_any = True
        r = t.iloc[0]
        loc = f"HLA-{locus.upper()}"
        cn1 = pd.to_numeric(r.get("HLA_type1copyNum_withBAFBin"), errors="coerce")
        cn2 = pd.to_numeric(r.get("HLA_type2copyNum_withBAFBin"), errors="coerce")
        pv = pd.to_numeric(r.get("PVal"), errors="coerce")
        pu = pd.to_numeric(r.get("PVal_unique"), errors="coerce")
        sites = r.get("numMisMatchSitesCov")
        rows.append({
            "sample": sid, "cohort": m.cohort, "locus": loc,
            "role": "TARGET" if loc == target_of.get(sid) else "control",
            "allele1": str(r.get("HLA_A_type1", "")).replace("hla_", ""),
            "allele2": str(r.get("HLA_A_type2", "")).replace("hla_", ""),
            "cn1": cn1, "cn2": cn2,
            "cn_ratio": (min(cn1, cn2) / max(cn1, cn2)) if pd.notna(cn1) and pd.notna(cn2) and max(cn1, cn2) > 0 else None,
            "pval": pv, "pval_unique": pu, "sites": sites,
            "detected": bool(pd.notna(pv) and pv < ALPHA),
        })

    if not got_any:
        reasons = set()
        for locus in ["a", "b", "c"]:
            e = f"{LOH}/{sid}/lohhla_{locus}.err"
            if not os.path.exists(e):
                continue
            for line in open(e, errors="ignore"):
                if line.startswith("Error"):
                    if "t.test" in line or "constant" in line:
                        reasons.add("no positions distinguishing the alleles")
                    elif "indelTotals" in line:
                        reasons.add("no indels for the edit-distance step")
                    elif "combinedTable" in line:
                        reasons.add("no locus produced a table")
                    else:
                        reasons.add(line.strip()[:50])
                    break
        failed.append((sid, m.cohort, target_of.get(sid, "?"),
                       "; ".join(sorted(reasons)) or "no error recorded"))

a = pd.DataFrame(rows)
if a.empty:
    print("no LOHHLA results")
    sys.exit(0)

print("=" * 104)
print(f" LOHHLA RESULTS - {a['sample'].nunique()} of {len(man)} samples completed")
print("=" * 104)
print()

# ---------------- per sample ----------------
for sid in sorted(a["sample"].unique()):
    g = a[a["sample"] == sid]
    coh = g.cohort.iloc[0]
    tgt = target_of.get(sid, "?")
    al = alleles_of.get(sid, "")
    print(f"  {sid}  ({coh})   thinned: {tgt} {al}")
    for _, r in g.iterrows():
        mark = "  <-- target" if r.role == "TARGET" else ""
        if r.role == "TARGET":
            verdict = "DETECTED" if r.detected else "missed"
        else:
            verdict = "FALSE POSITIVE" if r.detected else "clean"
        pv = f"{r.pval:.4g}" if pd.notna(r.pval) else "NA"
        ratio = f"{r.cn_ratio:.2f}" if r.cn_ratio is not None else "  - "
        print(f"      {r.locus:<7} CN {r.cn1:>5.2f} / {r.cn2:<5.2f}  ratio {ratio}"
              f"  P {pv:>10}  sites {str(r.sites):>3}   {verdict}{mark}")
        print(f"              {r.allele1}  vs  {r.allele2}")
    print()

# ---------------- summary ----------------
print("=" * 104)
print(" SUMMARY")
print("=" * 104)

t = a[a.role == "TARGET"]
c = a[a.role == "control"]

print(f"\n  target loci   : {len(t)} tested, {int(t.detected.sum())} detected"
      f"  ({100*t.detected.mean():.0f}% sensitivity)")
print(f"  control loci  : {len(c)} tested, {int(c.detected.sum())} called"
      f"  ({100*(1-c.detected.mean()):.0f}% specificity)" if len(c) else "")

print(f"\n  detected targets, by number of distinguishing sites:")
det = t[t.detected].copy()
mis = t[~t.detected].copy()
for lbl, g in [("detected", det), ("missed", mis)]:
    s = pd.to_numeric(g.sites, errors="coerce").dropna()
    if len(s):
        print(f"    {lbl:<10} n={len(g):<3} sites: median {s.median():.0f}, "
              f"range {s.min():.0f}-{s.max():.0f}")

print(f"\n  by cohort:")
for coh, g in a.groupby("cohort"):
    gt = g[g.role == "TARGET"]
    gc = g[g.role == "control"]
    print(f"    {coh.upper():<5} targets {int(gt.detected.sum())}/{len(gt)},"
          f" controls {int(gc.detected.sum())}/{len(gc)} false")

print(f"\n  by locus:")
for loc, g in a.groupby("locus"):
    gt = g[g.role == "TARGET"]
    if len(gt):
        print(f"    {loc:<7} {int(gt.detected.sum())}/{len(gt)} detected as target,"
              f" {len(g[g.role=='control'])} used as control")

# ---------------- failures ----------------
if failed:
    print()
    print("=" * 104)
    print(f" NOT COMPLETED - {len(failed)} samples")
    print("=" * 104)
    print(f"\n  {'sample':<8}{'cohort':<8}{'target':<9}reason")
    for sid, coh, tgt, why in failed:
        print(f"  {sid:<8}{coh:<8}{tgt:<9}{why}")

    print(f"\n  reasons:")
    from collections import Counter
    for why, n in Counter(w for _, _, _, w in failed).most_common():
        print(f"    {n:>3}  {why}")

out = os.path.expanduser("~/immune_escape_project/results/lohhla_report.tsv")
os.makedirs(os.path.dirname(out), exist_ok=True)
a.to_csv(out, sep="\t", index=False)
print(f"\nwritten to {out}")
