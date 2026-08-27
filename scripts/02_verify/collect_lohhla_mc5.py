#!/usr/bin/env python3
"""
Collect the LOHHLA results produced at minCoverageFilter 5.

Step 9b wrote to out_<locus>; step 9d wrote to out_<locus>_mc5 with the
same BAMs and the same allele references. This gathers the second set into
the same shape as the first so the two can be compared directly.

Usage:
  python collect_lohhla_mc5.py [workspace] [suffix]

  python collect_lohhla_mc5.py $WS _mc5
"""
import sys, os, glob
import numpy as np
import pandas as pd

args = [a for a in sys.argv[1:] if not a.startswith("--")]
WS = args[0] if args else os.environ.get(
    "WS", "/gpfs/bwfor/work/ws/fr_os136-immune_escape")
SUFFIX = args[1] if len(args) > 1 else "_mc5"

PANEL = f"{WS}/simulation/step9_lohhla_panel"
DESIGN = f"{PANEL}/collections.tsv"
MANIFEST = f"{WS}/simulation/cohort/manifest.tsv"
OUT = f"{PANEL}/results"
ALPHA = 0.05

COLLECTIONS = ["A_untouched", "B_loss20", "C_loss35", "D_loss50"]

ERRORS = {
    "t.test": "no positions distinguishing the alleles",
    "constant": "no positions distinguishing the alleles",
    "indelTotals": "no indels for the edit-distance step",
    "editDistance": "no indels for the edit-distance step",
    "combinedTable": "no locus produced a table",
}

os.makedirs(OUT, exist_ok=True)
design = pd.read_csv(DESIGN, sep="\t")
man = pd.read_csv(MANIFEST, sep="\t")

target_of = {(r.collection, r["sample"]): r.target_locus
             for _, r in design.iterrows()}
het_of = {(r.collection, r["sample"]): str(r.het_loci)
          for _, r in design.iterrows()}

print(f"reading out_<locus>{SUFFIX}\n")
rows = []

for col in COLLECTIONS:
    for _, m in man.iterrows():
        sid = m.sample_id
        key = (col, sid)
        if key not in target_of:
            continue
        tgt, het = target_of[key], het_of[key]

        for locus in ["a", "b", "c"]:
            loc_name = f"HLA-{locus.upper()}"
            if locus.upper() not in het:
                continue

            d = f"{PANEL}/{col}/{sid}/out_{locus}{SUFFIX}"
            preds = [f for f in glob.glob(f"{d}/*HLAlossPrediction*.txt")
                     if os.path.getsize(f) > 200
                     and "homozygous" not in f and "No_Suitable" not in f]
            role = "TARGET" if loc_name == tgt else "control"

            if not preds:
                reason = "no error recorded"
                for cand in [f"{PANEL}/{col}/{sid}/lohhla_{locus}{SUFFIX}.err",
                             f"{PANEL}/{col}/{sid}/lohhla_{locus}.err"]:
                    if not os.path.exists(cand):
                        continue
                    for line in open(cand, errors="ignore"):
                        if line.startswith("Error"):
                            for k, txt in ERRORS.items():
                                if k in line:
                                    reason = txt
                                    break
                            else:
                                reason = line.strip()[:50]
                            break
                    break
                rows.append({"collection": col, "sample": sid,
                             "cohort": m.cohort, "locus": loc_name,
                             "role": role, "status": "no result",
                             "reason": reason})
                continue

            try:
                t = pd.read_csv(preds[0], sep="\t")
            except Exception:
                continue
            if not len(t):
                continue

            r = t.iloc[0]
            cn1 = pd.to_numeric(r.get("HLA_type1copyNum_withBAFBin"), errors="coerce")
            cn2 = pd.to_numeric(r.get("HLA_type2copyNum_withBAFBin"), errors="coerce")
            pv = pd.to_numeric(r.get("PVal"), errors="coerce")
            sites = pd.to_numeric(r.get("numMisMatchSitesCov"), errors="coerce")

            rows.append({
                "collection": col, "sample": sid, "cohort": m.cohort,
                "locus": loc_name, "role": role, "status": "ok",
                "allele1": str(r.get("HLA_A_type1", "")).replace("hla_", ""),
                "allele2": str(r.get("HLA_A_type2", "")).replace("hla_", ""),
                "cn1": cn1, "cn2": cn2,
                "cn_ratio": (min(cn1, cn2) / max(cn1, cn2))
                            if pd.notna(cn1) and pd.notna(cn2) and max(cn1, cn2) > 0
                            else np.nan,
                "pval": pv, "sites": sites,
                "detected": bool(pd.notna(pv) and pv < ALPHA),
            })

L = pd.DataFrame(rows)
if L.empty:
    print("nothing found")
    sys.exit(1)

L["detected"] = L.detected.fillna(False).astype(bool)
ok = L[L.status == "ok"]

print(f"  locus runs      {len(L)}")
print(f"  produced        {len(ok)}")
print(f"  targets         {len(ok[ok.role=='TARGET'])}")
print(f"  controls        {len(ok[ok.role=='control'])}")

t = ok[ok.role == "TARGET"]
c = ok[ok.role == "control"]
if len(t):
    print(f"\n  sensitivity  {int(t.detected.sum())}/{len(t)}"
          f"  ({100*t.detected.mean():.1f}%)")
if len(c):
    print(f"  specificity  {int((~c.detected).sum())}/{len(c)}"
          f"  ({100*(1-c.detected.mean()):.1f}%)")

L.to_csv(f"{OUT}/step9_all_loci{SUFFIX}.tsv", sep="\t", index=False)
ok.to_csv(f"{OUT}/step9_predictions{SUFFIX}.tsv", sep="\t", index=False)
print(f"\nwritten to {OUT}/step9_all_loci{SUFFIX}.tsv")
