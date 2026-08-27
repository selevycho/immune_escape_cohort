#!/usr/bin/env python3
"""
Compare LOHHLA at two coverage thresholds on identical data.

The same collections were run at minCoverageFilter 3 and 5. Nothing else
differs - same BAMs, same allele references, same script - so the
comparison isolates what the threshold does.

Three things are worth knowing. Whether the looser setting reaches the
statistics on more loci, which it should, since it accepts positions
covered less deeply. Whether the extra positions carry signal or noise,
which shows in whether the number of discriminating positions still
predicts detection. And whether the looser setting calls loss where none
was made.

Usage:
  python compare_lohhla_thresholds.py [workspace]
"""
import sys, os
import numpy as np
import pandas as pd
from scipy import stats

args = [a for a in sys.argv[1:] if not a.startswith("--")]
WS = args[0] if args else os.environ.get(
    "WS", "/gpfs/bwfor/work/ws/fr_os136-immune_escape")

RES = f"{WS}/simulation/step9_lohhla_panel/results"

frames = {}
for tag, label in [("_mc3", "threshold 3"), ("_mc5", "threshold 5")]:
    p = f"{RES}/step9_all_loci{tag}.tsv"
    if os.path.exists(p):
        d = pd.read_csv(p, sep="\t")
        d["setting"] = label
        frames[label] = d
    else:
        print(f"  missing: {p}")

if len(frames) < 2:
    print("\nrun step9c_analyse_panel.py at both thresholds first:")
    print("  python step9c_analyse_panel.py $WS")
    print("  python step9c_analyse_panel.py $WS --mincov=5")
    sys.exit(1)

A = pd.concat(frames.values())
A["detected"] = A.detected.fillna(False).astype(bool)
ok = A[A.status == "ok"].copy()

print("=" * 92)
print(" HOW MANY LOCI REACHED THE STATISTICS")
print("=" * 92)
print(f"\n  {'setting':<16}{'attempted':>11}{'produced':>11}{'rate':>8}")
for lab, d in frames.items():
    o = d[d.status == "ok"]
    print(f"  {lab:<16}{len(d):>11}{len(o):>11}{100*len(o)/max(1,len(d)):>7.0f}%")

print()
print("=" * 92)
print(" SPECIFICITY")
print("=" * 92)
print(f"\n  {'setting':<16}{'controls':>10}{'false':>8}{'rate':>9}")
for lab in frames:
    c = ok[(ok.setting == lab) & (ok.role == "control")]
    if len(c):
        print(f"  {lab:<16}{len(c):>10}{int(c.detected.sum()):>8}"
              f"{100*c.detected.mean():>8.1f}%")

print(f"\n  collection A alone, where nothing was altered anywhere:")
for lab in frames:
    c = ok[(ok.setting == lab) & (ok.collection == "A_untouched")]
    if len(c):
        print(f"    {lab:<16}{int(c.detected.sum())} of {len(c)}")

print()
print("=" * 92)
print(" SENSITIVITY")
print("=" * 92)
print(f"\n  {'setting':<16}{'kept':<8}{'targets':>9}{'detected':>10}{'rate':>8}")
for lab in frames:
    for col, keep in [("B_loss20", 0.20), ("C_loss35", 0.35), ("D_loss50", 0.50)]:
        t = ok[(ok.setting == lab) & (ok.collection == col) &
               (ok.role == "TARGET")]
        if len(t):
            print(f"  {lab:<16}{keep:<8.2f}{len(t):>9}"
                  f"{int(t.detected.sum()):>10}{100*t.detected.mean():>7.1f}%")
    t = ok[(ok.setting == lab) & (ok.role == "TARGET")]
    print(f"  {'':<16}{'all':<8}{len(t):>9}{int(t.detected.sum()):>10}"
          f"{100*t.detected.mean():>7.1f}%")
    print()

print("=" * 92)
print(" DO THE EXTRA POSITIONS CARRY SIGNAL?")
print("=" * 92)
for lab in frames:
    t = ok[(ok.setting == lab) & (ok.role == "TARGET")]
    det = pd.to_numeric(t[t.detected].sites, errors="coerce").dropna()
    mis = pd.to_numeric(t[~t.detected].sites, errors="coerce").dropna()
    print(f"\n  {lab}")
    if len(det):
        print(f"    detected  n={len(det):<4} positions: median {det.median():.0f}, "
              f"range {det.min():.0f}-{det.max():.0f}")
    if len(mis):
        print(f"    missed    n={len(mis):<4} positions: median {mis.median():.0f}, "
              f"range {mis.min():.0f}-{mis.max():.0f}")
    if len(det) > 2 and len(mis) > 2:
        u, p = stats.mannwhitneyu(det, mis, alternative="greater")
        verdict = ("positions still predict detection" if p < 0.05
                   else "the association has gone")
        print(f"    Mann-Whitney p = {p:.4f}  -  {verdict}")

print()
print("=" * 92)
print(" DETECTION RATE BY NUMBER OF POSITIONS")
print("=" * 92)
bins = [0, 3, 5, 10, 20, 10000]
labs = ["0-2", "3-4", "5-9", "10-19", "20+"]
print(f"\n  {'positions':<12}", end="")
for lab in frames:
    print(f"{lab:>18}", end="")
print()
for i, lb in enumerate(labs):
    line = f"  {lb:<12}"
    for lab in frames:
        t = ok[(ok.setting == lab) & (ok.role == "TARGET")].copy()
        t["bin"] = pd.cut(pd.to_numeric(t.sites, errors="coerce"),
                          bins=bins, labels=labs, right=False)
        g = t[t.bin == lb]
        if len(g):
            line += f"{int(g.detected.sum()):>8} of {len(g):<3}"
        else:
            line += f"{'-':>18}"
    print(line)

print(f"\n  If more positions mean more information, detection should rise")
print(f"  down the column. Where it does not, the extra positions admitted")
print(f"  by the looser threshold are adding variance rather than signal.")

print()
print("=" * 92)
print(" VERDICT")
print("=" * 92)
best = None
for lab in frames:
    t = ok[(ok.setting == lab) & (ok.role == "TARGET")]
    c = ok[(ok.setting == lab) & (ok.role == "control")]
    if not len(t) or not len(c):
        continue
    sens = 100 * t.detected.mean()
    spec = 100 * (1 - c.detected.mean())
    print(f"\n  {lab}")
    print(f"    sensitivity {int(t.detected.sum())}/{len(t)}  ({sens:.1f}%)")
    print(f"    specificity {int((~c.detected).sum())}/{len(c)}  ({spec:.1f}%)")

out = f"{RES}/threshold_comparison.tsv"
A.to_csv(out, sep="\t", index=False)
print(f"\nwritten to {out}")
