#!/usr/bin/env python3
"""
Did the LOH simulation build correctly, and did LOHHLA find what was put there?

Two questions, answered separately.

Step 6a is checked by construction: one locus per sample was thinned to a
known fraction and the other two were left alone. The check is arithmetic -
does the observed depth ratio match what was asked for, and did the control
loci stay where they were.

Step 6b is checked as a detection problem. Each sample gives up to three
independent tests, one per heterozygous locus. The locus that was thinned
should come back significant; the loci that were not should not. Counting
both gives sensitivity and the false positive rate on the same data, which
is the point of simulating in the first place.

Usage:
  python check_lohhla.py [workspace]
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
print(f"cohort: {len(man)} samples\n")

# =====================================================================
# STEP 6a
# =====================================================================
print("=" * 92)
print(" STEP 6a - was the loss simulated as designed?")
print("=" * 92)

if not os.path.exists(DESIGN):
    print("  loh_design.tsv missing - step 6a has not run")
    sys.exit(1)

d = pd.read_csv(DESIGN, sep="\t")
print(f"  samples in design : {len(d)} of {len(man)}")

missing = sorted(set(man.sample_id) - set(d["sample"]))
if missing:
    print(f"  NOT SIMULATED     : {' '.join(missing)}")

print(f"\n  target locus distribution:")
for loc, n in d.target_locus.value_counts().items():
    print(f"    {loc:<8} {n}")

# depth ratios
print(f"\n  {'sample':<7}{'target':<8}{'before':>8}{'after':>8}{'ratio':>8}   controls")
bad_target, bad_control = [], []
for _, r in d.iterrows():
    loc = r.target_locus.replace("HLA-", "")
    b = r.get(f"depth_{loc}_before")
    a = r.get(f"depth_{loc}_after")
    if pd.isna(b) or b == 0:
        continue
    ratio = a / b
    want = r.keep_fraction

    ctrl_txt = []
    for c in "ABC":
        if c == loc:
            continue
        cb = r.get(f"depth_{c}_before")
        ca = r.get(f"depth_{c}_after")
        if pd.isna(cb):
            continue
        if abs(ca - cb) > 0.5:
            ctrl_txt.append(f"{c} MOVED {cb:.1f}->{ca:.1f}")
            bad_control.append((r["sample"], c, cb, ca))
        else:
            ctrl_txt.append(f"{c} {ca:.1f}")

    flag = ""
    if abs(ratio - want) > 0.05:
        flag = "  <-- off target"
        bad_target.append((r["sample"], ratio, want))

    print(f"  {r['sample']:<7}{r.target_locus:<8}{b:>8.1f}{a:>8.1f}{ratio:>8.2f}   "
          f"{', '.join(ctrl_txt)}{flag}")

print()
if not bad_target and not bad_control and not missing:
    print("  6a: all samples thinned to the requested fraction, controls untouched")
else:
    if missing:
        print(f"  6a: {len(missing)} sample(s) never simulated")
    if bad_target:
        print(f"  6a: {len(bad_target)} sample(s) missed the target ratio")
    if bad_control:
        print(f"  6a: {len(bad_control)} control locus/loci moved when they should not have")

# =====================================================================
# STEP 6b
# =====================================================================
print()
print("=" * 92)
print(" STEP 6b - did LOHHLA detect it?")
print("=" * 92)

target_of = dict(zip(d["sample"], d.target_locus))

rows = []
for sid in man.sample_id:
    for locus in ["a", "b", "c"]:
        preds = glob.glob(f"{LOH}/{sid}/out_{locus}/*HLAlossPrediction*.txt")
        preds = [p for p in preds
                 if os.path.getsize(p) > 200
                 and "homozygous" not in p and "No_Suitable" not in p]
        if not preds:
            continue
        try:
            t = pd.read_csv(preds[0], sep="\t")
        except Exception:
            continue
        if not len(t):
            continue
        r = t.iloc[0]
        loc_name = f"HLA-{locus.upper()}"
        rows.append({
            "sample": sid,
            "cohort": "brca" if sid[0] == "B" else "ov",
            "locus": loc_name,
            "is_target": loc_name == target_of.get(sid),
            "allele1": r.get("HLA_A_type1"),
            "allele2": r.get("HLA_A_type2"),
            "cn1": r.get("HLA_type1copyNum_withBAFBin"),
            "cn2": r.get("HLA_type2copyNum_withBAFBin"),
            "pval": r.get("PVal"),
            "pval_unique": r.get("PVal_unique"),
            "sites": r.get("numMisMatchSitesCov"),
            "loss_allele": r.get("LossAllele"),
        })

if not rows:
    print("  no LOHHLA predictions at all")
    sys.exit(0)

a = pd.DataFrame(rows)
a["detected"] = pd.to_numeric(a.pval, errors="coerce") < ALPHA

print(f"  locus-level tests completed: {len(a)}")
print(f"  samples with at least one  : {a['sample'].nunique()} of {len(man)}")

no_run = sorted(set(man.sample_id) - set(a["sample"]))
if no_run:
    print(f"  no result for              : {' '.join(no_run)}")

# --- the core table ---
print()
print("  " + "-" * 88)
print(f"  {'':<24}{'tests':>8}{'detected':>10}{'rate':>8}{'median P':>12}{'median sites':>14}")
print("  " + "-" * 88)
for is_t, label in [(True, "thinned locus (target)"), (False, "untouched locus (control)")]:
    g = a[a.is_target == is_t]
    if not len(g):
        continue
    det = int(g.detected.sum())
    mp = pd.to_numeric(g.pval, errors="coerce").median()
    ms = pd.to_numeric(g.sites, errors="coerce").median()
    print(f"  {label:<24}{len(g):>8}{det:>10}{100*det/len(g):>7.1f}%{mp:>12.4f}{ms:>14.0f}")
print("  " + "-" * 88)

tp = a[(a.is_target) & (a.detected)]
fn = a[(a.is_target) & (~a.detected)]
fp = a[(~a.is_target) & (a.detected)]
tn = a[(~a.is_target) & (~a.detected)]

print()
print(f"  sensitivity : {len(tp)}/{len(tp)+len(fn)}"
      f"  ({100*len(tp)/max(1,len(tp)+len(fn)):.1f}%)   loss found where it was created")
print(f"  specificity : {len(tn)}/{len(tn)+len(fp)}"
      f"  ({100*len(tn)/max(1,len(tn)+len(fp)):.1f}%)   no loss called where none was made")

# --- per sample ---
print()
print("=" * 92)
print(" PER LOCUS")
print("=" * 92)
print(f"  {'sample':<7}{'locus':<8}{'role':<9}{'CN1':>7}{'CN2':>7}{'PVal':>10}"
      f"{'sites':>7}  {'result':<10} allele pair")
for _, r in a.sort_values(["sample", "locus"]).iterrows():
    role = "TARGET" if r.is_target else "control"
    if r.is_target:
        res = "detected" if r.detected else "MISSED"
    else:
        res = "FALSE POS" if r.detected else "clean"
    cn1 = r.cn1 if pd.notna(r.cn1) else 0
    cn2 = r.cn2 if pd.notna(r.cn2) else 0
    pv = pd.to_numeric(r.pval, errors="coerce")
    print(f"  {r['sample']:<7}{r.locus:<8}{role:<9}{cn1:>7.2f}{cn2:>7.2f}"
          f"{pv:>10.4f}{str(r.sites):>7}  {res:<10} "
          f"{str(r.allele1).replace('hla_','')} / {str(r.allele2).replace('hla_','')}")

# --- by cohort ---
print()
print("=" * 92)
print(" BY COHORT")
print("=" * 92)
for coh, g in a.groupby("cohort"):
    t = g[g.is_target]
    c = g[~g.is_target]
    print(f"  {coh.upper():<6} targets {int(t.detected.sum())}/{len(t)} detected, "
          f"controls {int(c.detected.sum())}/{len(c)} false positives")

out = os.path.expanduser("~/immune_escape_project/results/lohhla_summary.tsv")
os.makedirs(os.path.dirname(out), exist_ok=True)
a.to_csv(out, sep="\t", index=False)
print(f"\nwritten to {out}")
