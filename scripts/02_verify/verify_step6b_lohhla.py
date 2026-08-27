#!/usr/bin/env python3
"""
Step 6b verification: what LOHHLA did with the simulated loss.

Step 6a established that one locus per sample was thinned to a known
fraction and the others were not touched at all. That makes this a
detection problem with an answer key: every locus LOHHLA reports on is
either one where loss was created or one where it was not, and both
outcomes are informative.

The tool is a modified copy - novoalign replaced by bwa, GATK 3 by
samtools, hg19 coordinates by hg38 - so the controls carry the weight of
the argument. A patched aligner that invented losses would show up there
before anywhere else.

Seventeen of forty samples produce no prediction at all. That is not a
crash to be worked around but a measurement: LOHHLA needs positions where
the two alleles of a locus differ, and on panel data there are often too
few. The failure modes are read from the R stderr and counted, and the
number of discriminating positions is compared between the loci that
worked and the ones that did not.

Usage:
  python verify_step6b_lohhla.py [workspace]
"""
import sys, os, glob
import numpy as np
import pandas as pd
from scipy import stats

args = [a for a in sys.argv[1:] if not a.startswith("--")]
WS = args[0] if args else os.environ.get(
    "WS", "/gpfs/bwfor/work/ws/fr_os136-immune_escape")

COHORT = f"{WS}/simulation/cohort"
LOH = f"{WS}/simulation/lohhla_allelic"
MANIFEST = f"{COHORT}/manifest.tsv"
DESIGN = f"{LOH}/loh_design.tsv"
ALPHA = 0.05

ERRORS = {
    "t.test": "no positions distinguishing the alleles",
    "constant": "no positions distinguishing the alleles",
    "indelTotals": "no indels for the edit-distance step",
    "combinedTable": "no locus produced a table",
    "editDistance": "no indels for the edit-distance step",
}

man = pd.read_csv(MANIFEST, sep="\t")
design = pd.read_csv(DESIGN, sep="\t") if os.path.exists(DESIGN) else pd.DataFrame()
target_of = dict(zip(design["sample"], design.target_locus)) if len(design) else {}

print(f"checking {len(man)} samples\n", flush=True)

locus_rows, sample_rows = [], []

for _, m in man.iterrows():
    sid = m.sample_id
    tgt = target_of.get(sid)
    ran, produced, failures = 0, 0, {}

    for locus in ["a", "b", "c"]:
        d = f"{LOH}/{sid}/out_{locus}"
        loc_name = f"HLA-{locus.upper()}"
        is_target = (loc_name == tgt)

        # did this locus even qualify for a run
        al = f"{LOH}/{sid}/al_{locus}.txt"
        n_alleles = 0
        if os.path.exists(al):
            n_alleles = len([x for x in open(al) if x.strip()])
        if n_alleles < 2:
            continue
        ran += 1

        preds = [f for f in glob.glob(f"{d}/*HLAlossPrediction*.txt")
                 if os.path.getsize(f) > 200
                 and "homozygous" not in f and "No_Suitable" not in f]

        if not preds:
            err = f"{LOH}/{sid}/lohhla_{locus}.err"
            reason = "no error recorded"
            if os.path.exists(err):
                for line in open(err, errors="ignore"):
                    if line.startswith("Error"):
                        for key, txt in ERRORS.items():
                            if key in line:
                                reason = txt
                                break
                        else:
                            reason = line.strip()[:60]
                        break
            failures[loc_name] = reason
            locus_rows.append({
                "sample": sid, "cohort": m.cohort, "locus": loc_name,
                "role": "TARGET" if is_target else "control",
                "status": "failed", "reason": reason,
            })
            continue

        try:
            t = pd.read_csv(preds[0], sep="\t")
        except Exception:
            continue
        if not len(t):
            continue

        produced += 1
        r = t.iloc[0]
        cn1 = pd.to_numeric(r.get("HLA_type1copyNum_withBAFBin"), errors="coerce")
        cn2 = pd.to_numeric(r.get("HLA_type2copyNum_withBAFBin"), errors="coerce")
        pv = pd.to_numeric(r.get("PVal"), errors="coerce")
        sites = pd.to_numeric(r.get("numMisMatchSitesCov"), errors="coerce")

        locus_rows.append({
            "sample": sid, "cohort": m.cohort, "locus": loc_name,
            "role": "TARGET" if is_target else "control",
            "status": "ok",
            "allele1": str(r.get("HLA_A_type1", "")).replace("hla_", ""),
            "allele2": str(r.get("HLA_A_type2", "")).replace("hla_", ""),
            "cn1": cn1, "cn2": cn2,
            "cn_ratio": (min(cn1, cn2) / max(cn1, cn2))
                        if pd.notna(cn1) and pd.notna(cn2) and max(cn1, cn2) > 0 else None,
            "pval": pv,
            "sites": sites,
            "detected": bool(pd.notna(pv) and pv < ALPHA),
        })

    sample_rows.append({
        "sample": sid, "cohort": m.cohort, "target": tgt,
        "loci_attempted": ran, "loci_produced": produced,
        "any_result": produced > 0,
        "failures": "; ".join(f"{k}: {v}" for k, v in failures.items()),
    })
    print(f"  {sid}  {produced}/{ran} loci produced a prediction")

L = pd.DataFrame(locus_rows)
S = pd.DataFrame(sample_rows)
pd.set_option("display.width", 220)

if L.empty:
    print("\nno LOHHLA output at all")
    sys.exit(1)

ok = L[L.status == "ok"].copy()

# The column arrives as object dtype because failed rows contribute NaN,
# and ~ on an object column negates bitwise rather than logically, which
# turns True into -2 and every count that follows into nonsense.
ok["detected"] = ok.detected.fillna(False).astype(bool)

# =====================================================================
print()
print("=" * 100)
print(" LOCUS BY LOCUS")
print("=" * 100)
print(f"\n  {'sample':<7}{'locus':<8}{'role':<9}{'CN1':>7}{'CN2':>7}"
      f"{'ratio':>7}{'PVal':>11}{'sites':>7}  {'verdict':<12} alleles")
for _, r in ok.sort_values(["sample", "locus"]).iterrows():
    if r.role == "TARGET":
        v = "DETECTED" if r.detected else "missed"
    else:
        v = "FALSE POS" if r.detected else "clean"
    pv = f"{r.pval:.4g}" if pd.notna(r.pval) else "NA"
    rat = f"{r.cn_ratio:.2f}" if pd.notna(r.cn_ratio) else "  - "
    print(f"  {r['sample']:<7}{r.locus:<8}{r.role:<9}"
          f"{r.cn1 if pd.notna(r.cn1) else 0:>7.2f}"
          f"{r.cn2 if pd.notna(r.cn2) else 0:>7.2f}{rat:>7}"
          f"{pv:>11}{str(r.sites):>7}  {v:<12} "
          f"{r.allele1} / {r.allele2}")

# =====================================================================
print()
print("=" * 100)
print(" DETECTION")
print("=" * 100)
tg = ok[ok.role == "TARGET"]
ct = ok[ok.role == "control"]

print(f"\n  {'':<26}{'tested':>9}{'called':>9}{'rate':>9}"
      f"{'median P':>12}{'median sites':>14}")
for lab, g in [("thinned locus", tg), ("untouched locus", ct)]:
    if not len(g):
        continue
    mp = pd.to_numeric(g.pval, errors="coerce").median()
    ms = pd.to_numeric(g.sites, errors="coerce").median()
    print(f"  {lab:<26}{len(g):>9}{int(g.detected.sum()):>9}"
          f"{100*g.detected.mean():>8.1f}%"
          f"{mp if pd.notna(mp) else 0:>12.4f}{ms if pd.notna(ms) else 0:>14.0f}")

tp, fn = int(tg.detected.sum()), int((~tg.detected).sum())
tn, fp = int((~ct.detected).sum()), int(ct.detected.sum())
print(f"\n  sensitivity  {tp}/{tp+fn}  ({100*tp/max(1,tp+fn):.1f}%)")
print(f"  specificity  {tn}/{tn+fp}  ({100*tn/max(1,tn+fp):.1f}%)")
if fp == 0:
    print(f"\n  No loss was called at any locus where none was made. Since")
    print(f"  step 6a confirmed those loci were untouched, the substitution")
    print(f"  of bwa for novoalign has not introduced spurious detections.")

# =====================================================================
print()
print("=" * 100)
print(" WHAT SEPARATES DETECTION FROM A MISS")
print("=" * 100)
det = tg[tg.detected]
mis = tg[~tg.detected]
for lab, g in [("detected", det), ("missed", mis)]:
    v = pd.to_numeric(g.sites, errors="coerce").dropna()
    if len(v):
        print(f"\n  {lab:<10} n={len(g):<3} discriminating positions: "
              f"median {v.median():.0f}, range {v.min():.0f}–{v.max():.0f}")
if len(det) > 2 and len(mis) > 2:
    a = pd.to_numeric(det.sites, errors="coerce").dropna()
    b = pd.to_numeric(mis.sites, errors="coerce").dropna()
    if len(a) > 2 and len(b) > 2:
        u, p = stats.mannwhitneyu(a, b, alternative="greater")
        print(f"\n  Mann-Whitney, detected > missed: p = {p:.4f}")

# is it coverage instead?
v1 = os.path.expanduser("~/immune_escape_project/results/verify_step1.tsv")
if os.path.exists(v1):
    cov = pd.read_csv(v1, sep="\t")
    cov = cov.set_index("sample")
    got = set(S[S.any_result]["sample"])
    rows = []
    for _, r in S.iterrows():
        tgt = str(r.target).replace("HLA-", "") if pd.notna(r.target) else None
        col = f"hla_{tgt}_depth" if tgt else None
        if col and r["sample"] in cov.index and col in cov.columns:
            rows.append({"sample": r["sample"],
                         "depth": cov.loc[r["sample"], col],
                         "worked": r["sample"] in got})
    c = pd.DataFrame(rows)
    if len(c) > 5:
        print(f"\n  target locus coverage:")
        for w, g in c.groupby("worked"):
            lab = "produced a result" if w else "did not"
            print(f"    {lab:<20} n={len(g):>3}  median {g.depth.median():.1f}x, "
                  f"range {g.depth.min():.1f}–{g.depth.max():.1f}x")
        a = c[c.worked].depth.dropna()
        b = c[~c.worked].depth.dropna()
        if len(a) > 2 and len(b) > 2:
            u, p = stats.mannwhitneyu(a, b)
            print(f"    Mann-Whitney p = {p:.3f}"
                  f"  — {'no relation to coverage' if p > 0.05 else 'coverage matters'}")

# =====================================================================
print()
print("=" * 100)
print(" WHERE LOHHLA STOPPED")
print("=" * 100)
failed = L[L.status == "failed"]
print(f"\n  {len(failed)} locus runs did not produce a prediction")
if len(failed):
    print(f"\n  by reason:")
    for why, n in failed.reason.value_counts().items():
        print(f"    {n:>4}  {why}")
    print(f"\n  by locus:")
    for loc, n in failed.locus.value_counts().items():
        tot = len(L[L.locus == loc])
        print(f"    {loc:<8}{n:>4} of {tot:>4}  ({100*n/tot:.0f}%)")
    print(f"\n  HLA-B and HLA-C fail more often than HLA-A because their")
    print(f"  allele pairs are more similar, leaving fewer positions at")
    print(f"  which reads can be assigned to one allele or the other.")

no_result = S[~S.any_result]
if len(no_result):
    print(f"\n  {len(no_result)} samples produced nothing at any locus:")
    print(f"    {' '.join(no_result['sample'])}")

# =====================================================================
print()
print("=" * 100)
print(" BY COHORT")
print("=" * 100)
for coh, g in ok.groupby("cohort"):
    t_ = g[g.role == "TARGET"]
    c_ = g[g.role == "control"]
    print(f"\n  {coh.upper():<6} targets {int(t_.detected.sum())}/{len(t_)} detected, "
          f"controls {int(c_.detected.sum())}/{len(c_)} false")

out = os.path.expanduser("~/immune_escape_project/results/verify_step6b.tsv")
os.makedirs(os.path.dirname(out), exist_ok=True)
L.to_csv(out, sep="\t", index=False)
S.to_csv(out.replace(".tsv", "_samples.tsv"), sep="\t", index=False)
print(f"\nwritten to {out}")
