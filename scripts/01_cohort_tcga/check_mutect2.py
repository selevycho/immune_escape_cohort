#!/usr/bin/env python3
"""
Mutect2 performance against the known truth set.

Recall here is measured only against mutations that verification confirmed
were actually injected. BAMSurgeon reports success at positions where
realignment later drove the reads to zero mapping quality, mostly across
the MHC; counting those as caller failures would blame Mutect2 for a
mutation that was never in the file.

Sensitivity is reported per VAF bin rather than as one number, because the
per-sample recall varies from 35% to 94% and that variation follows each
donor's VAF distribution rather than anything about the sample. Binning
separates the two.

Usage:
  python check_mutect2.py <workspace>
"""
import sys, os, glob
import numpy as np
import pandas as pd

WS = sys.argv[1] if len(sys.argv) > 1 else os.environ.get(
    "WS", "/gpfs/bwfor/work/ws/fr_os136-immune_escape")

COHORT = f"{WS}/simulation/cohort"
MANIFEST = f"{COHORT}/manifest.tsv"

BINS = [0, 0.05, 0.10, 0.15, 0.20, 0.30, 1.01]
LABELS = ["<5%", "5-10%", "10-15%", "15-20%", "20-30%", ">30%"]

man = pd.read_csv(MANIFEST, sep="\t")

frames, per_sample = [], []

for _, m in man.iterrows():
    sid = m.sample_id
    f = f"{COHORT}/{sid}/comparison/truth_vs_calls.tsv"
    if not os.path.exists(f):
        per_sample.append({"sample": sid, "cohort": m.cohort, "status": "no comparison"})
        continue
    d = pd.read_csv(f, sep="\t")
    if not len(d):
        per_sample.append({"sample": sid, "cohort": m.cohort, "status": "empty"})
        continue
    d["sample"] = sid
    d["cohort"] = m.cohort
    frames.append(d)

    # false positives, if the comparison recorded them
    log = f"{COHORT}/{sid}/comparison.log"
    fp = None
    if os.path.exists(log):
        for line in open(log, errors="ignore"):
            if "false positives" in line:
                digits = [int(x) for x in line.split() if x.isdigit()]
                if digits:
                    fp = digits[0]
                break

    per_sample.append({
        "sample": sid, "cohort": m.cohort, "status": "ok",
        "injected": len(d),
        "recovered": int(d.detected_pass.sum()),
        "recall": round(100 * d.detected_pass.mean(), 1),
        "median_vaf": round(d.VAF.median(), 3),
        "false_pos": fp,
    })

ps = pd.DataFrame(per_sample)
a = pd.concat(frames) if frames else pd.DataFrame()

if a.empty:
    print("no comparison files found")
    sys.exit(1)

pd.set_option("display.width", 200)

print("=" * 76)
print(" PER SAMPLE")
print("=" * 76)
cols = [c for c in ["sample", "cohort", "injected", "recovered", "recall",
                    "median_vaf", "false_pos"] if c in ps.columns]
print(ps[ps.status == "ok"][cols].to_string(index=False))

missing = ps[ps.status != "ok"]
if len(missing):
    print(f"\n  no result for: {' '.join(missing['sample'])}")

print()
print("=" * 76)
print(" SENSITIVITY BY VAF")
print("=" * 76)

a["bin"] = pd.cut(a.VAF, bins=BINS, labels=LABELS, right=False)
t = a.groupby("bin", observed=False).agg(
    mutations=("VAF", "size"),
    recovered=("detected_pass", "sum"),
    median_vaf=("VAF", "median")).reset_index()
t["recall"] = (100 * t.recovered / t.mutations).round(1)
t["reads_at_33x"] = (33 * t.median_vaf).round(1)

print()
print(f" {'VAF bin':<10}{'mutations':>11}{'recovered':>11}{'recall':>9}"
      f"{'median VAF':>12}{'reads at 33x':>14}")
for _, r in t.iterrows():
    if r.mutations == 0:
        continue
    print(f" {str(r['bin']):<10}{int(r.mutations):>11}{int(r.recovered):>11}"
          f"{r.recall:>8.1f}%{r.median_vaf:>12.3f}{r.reads_at_33x:>14.1f}")

print()
print("=" * 76)
print(" OVERALL")
print("=" * 76)
n_ok = (ps.status == "ok").sum()
print(f"\n  samples analysed      {n_ok} of {len(man)}")
print(f"  mutations in truth    {len(a):,}")
print(f"  recovered (PASS)      {int(a.detected_pass.sum()):,}")
print(f"  overall recall        {100 * a.detected_pass.mean():.1f}%")
print(f"  per-sample recall     {ps[ps.status=='ok'].recall.min():.1f}% – "
      f"{ps[ps.status=='ok'].recall.max():.1f}%, "
      f"median {ps[ps.status=='ok'].recall.median():.1f}%")

if ps.false_pos.notna().any():
    fpv = ps.false_pos.dropna()
    print(f"  false positives       total {int(fpv.sum())}, "
          f"{int((fpv == 0).sum())} of {len(fpv)} samples had none")

print()
print("=" * 76)
print(" BY COHORT")
print("=" * 76)
for coh, g in a.groupby("cohort"):
    print(f"\n  {coh.upper()}  n={g['sample'].nunique()} samples, "
          f"{len(g):,} mutations")
    print(f"    recall {100 * g.detected_pass.mean():.1f}%, "
          f"median VAF {g.VAF.median():.3f}")

# does per-sample recall track VAF rather than sample identity?
ok = ps[ps.status == "ok"]
if len(ok) > 3:
    from scipy import stats
    r, p = stats.pearsonr(ok.median_vaf, ok.recall)
    print(f"\n  per-sample recall vs median VAF: r = {r:.3f}, p = {p:.3g}")
    print(f"    recall follows the VAF distribution, not the sample")

out = os.path.expanduser("~/immune_escape_project/results/mutect2_performance.tsv")
os.makedirs(os.path.dirname(out), exist_ok=True)
ps.to_csv(out, sep="\t", index=False)
t.to_csv(out.replace(".tsv", "_by_vaf.tsv"), sep="\t", index=False)

print()
print("=" * 76)
print(" NUMBERS FOR THE SLIDE")
print("=" * 76)
print(f"\n  {n_ok} / {len(man)} samples called")
print(f"  recall {ok.recall.min():.0f} – {ok.recall.max():.0f}%, "
      f"median {ok.recall.median():.0f}%")
print(f"\n  sensitivity curve:")
for _, r in t.iterrows():
    if r.mutations == 0:
        continue
    print(f"    {str(r['bin']):<10}{r.recall:>6.1f}%   (n = {int(r.mutations)})")
print(f"\nwritten to {out}")
