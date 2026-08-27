#!/usr/bin/env python3
"""
Where to draw the line on LOHHLA's p-value.

Every result so far calls loss at p < 0.05 because that is the convention,
not because anything about this data recommends it. The four collections
make the choice testable: collection A has no loss anywhere, so every
detection in it is false, and the loss collections have exactly one
altered locus each.

The trade is the usual one. A looser threshold finds more of the losses
that were made and more that were not. What matters is the rate at which
each rises, and whether some threshold sits at a corner where one moves
faster than the other.

Two cautions about what follows.

The sample is small - 59 targets and 31 controls at threshold 5 - so the
curve is coarse and its confidence intervals are wide. Bootstrap intervals
are computed rather than left implicit, because a specificity of 100% on
31 observations is compatible with a true rate near 90%.

The three loss collections share loci: the same sample contributes at
0.20, 0.35 and 0.50. Those observations are not independent, and pooling
them inflates the effective sample size. Sensitivity is therefore also
reported per collection, where each locus appears once.

Usage:
  python analyse_lohhla_roc.py [workspace] [--mincov=5]
"""
import sys, os
import numpy as np
import pandas as pd

args = [a for a in sys.argv[1:] if not a.startswith("--")]
WS = args[0] if args else os.environ.get(
    "WS", "/gpfs/bwfor/work/ws/fr_os136-immune_escape")

SUFFIX = "_mc5"
for a in sys.argv[1:]:
    if a.startswith("--mincov="):
        SUFFIX = "_mc" + a.split("=", 1)[1]

RES = f"{WS}/simulation/step9_lohhla_panel/results"
PATH = f"{RES}/step9_all_loci{SUFFIX}.tsv"

THRESHOLDS = [0.001, 0.005, 0.01, 0.025, 0.05, 0.10, 0.15, 0.20, 0.50]
N_BOOT = 2000
RNG = np.random.default_rng(7)

if not os.path.exists(PATH):
    print(f"missing: {PATH}")
    print(f"run collect_lohhla_mc5.py first")
    sys.exit(1)

d = pd.read_csv(PATH, sep="\t")
ok = d[d.status == "ok"].copy()
ok["pval"] = pd.to_numeric(ok.pval, errors="coerce")

# a locus with no p-value cannot be called at any threshold; keeping it
# would make every rate depend on how many such loci happened to appear
scored = ok[ok.pval.notna()].copy()
targets = scored[scored.role == "TARGET"]
controls = scored[scored.role == "control"]

print(f"reading {os.path.basename(PATH)}\n")
print(f"  loci with a prediction   {len(ok)}")
print(f"  of those, with a p-value {len(scored)}")
print(f"    targets                {len(targets)}")
print(f"    controls               {len(controls)}")

# Controls almost never carry a p-value. LOHHLA computes one by testing
# whether the two allele copy numbers differ, and at an untouched locus
# they do not differ at all - the test has nothing to work with and
# returns NA. That is the correct answer rather than a failure, but it
# means specificity here is not a threshold question: a control either
# produced a significant p or produced none, and no cutoff moves it.
n_ctrl_na = int((ok.role == "control").sum()) - len(controls)
if len(controls) < 5:
    print(f"\n  {n_ctrl_na} controls produced a table but no p-value.")
    print(f"  At an untouched locus the two copy numbers are equal, so the")
    print(f"  test has no difference to evaluate. Specificity therefore does")
    print(f"  not vary with the threshold, and the curve below describes")
    print(f"  sensitivity alone.")

if len(targets) < 5:
    print("\ntoo few scored targets for a curve")
    sys.exit(1)


def boot_ci(hits, n, reps=N_BOOT):
    """Percentile bootstrap interval for a proportion."""
    if n == 0:
        return (np.nan, np.nan)
    draws = RNG.binomial(n, hits / n, reps) / n
    return (100 * np.percentile(draws, 2.5),
            100 * np.percentile(draws, 97.5))


# =====================================================================
print()
print("=" * 96)
print(" SENSITIVITY AND SPECIFICITY ACROSS THRESHOLDS")
print("=" * 96)
print(f"\n  {'p <':<9}{'sens':>8}{'95% CI':>16}{'spec':>9}{'95% CI':>16}"
      f"{'Youden':>9}")

rows = []
for th in THRESHOLDS:
    tp = int((targets.pval < th).sum())
    fn = len(targets) - tp
    fp = int((controls.pval < th).sum()) if len(controls) else 0
    tn = len(controls) - fp

    sens = tp / len(targets)
    spec = (tn / len(controls)) if len(controls) else float("nan")
    lo_s, hi_s = boot_ci(tp, len(targets))
    lo_p, hi_p = boot_ci(tn, len(controls)) if len(controls) else (float("nan"),) * 2

    rows.append({
        "threshold": th, "tp": tp, "fn": fn, "fp": fp, "tn": tn,
        "sensitivity": 100 * sens, "specificity": 100 * spec,
        "sens_lo": lo_s, "sens_hi": hi_s,
        "spec_lo": lo_p, "spec_hi": hi_p,
        "youden": sens + spec - 1,
        "ppv": 100 * tp / (tp + fp) if (tp + fp) else np.nan,
    })

    print(f"  {th:<9.3f}{100*sens:>7.1f}%{f'{lo_s:.0f}-{hi_s:.0f}':>16}"
          f"{100*spec:>8.1f}%{f'{lo_p:.0f}-{hi_p:.0f}':>16}"
          f"{sens+spec-1:>9.3f}")

R = pd.DataFrame(rows)
# Youden's J needs both halves. With no scored controls it is undefined,
# and the threshold cannot be chosen by balancing two rates when only one
# of them varies.
if R.youden.notna().any():
    best = R.loc[R.youden.idxmax()]
    print(f"\n  Youden's J peaks at p < {best.threshold:.3f}: "
          f"sensitivity {best.sensitivity:.1f}%, "
          f"specificity {best.specificity:.1f}%")
else:
    print(f"\n  Youden's J is undefined here: specificity does not vary")
    print(f"  with the threshold, so there is no balance point to find.")
    print(f"  The threshold has to be argued from the sensitivity curve")
    print(f"  and from what a false positive would cost.")

conv = R[R.threshold == 0.05]
if len(conv):
    c = conv.iloc[0]
    print(f"  the conventional 0.05 gives sensitivity {c.sensitivity:.1f}%")
    if R.youden.notna().any() and abs(best.threshold - 0.05) < 1e-9:
        print(f"\n  The convention and the optimum coincide, which is worth")
        print(f"  stating: the threshold was not chosen to flatter the result.")

# =====================================================================
print()
print("=" * 96)
print(" WHAT THE INTERVALS MEAN")
print("=" * 96)
c = R[R.threshold == 0.05].iloc[0]
print(f"\n  At p < 0.05 sensitivity is {c.sensitivity:.1f}% "
      f"on {len(targets)} scored targets.")
print(f"\n  The bootstrap intervals are {c.sens_lo:.0f}-{c.sens_hi:.0f}% and "
      f"{c.spec_lo:.0f}-{c.spec_hi:.0f}%.")
if len(controls) and c.specificity == 100:
    # rule of three: zero events in n trials bounds the true rate at 3/n
    print(f"  Zero false positives in {len(controls)} controls bounds the")
    print(f"  true rate below {300/len(controls):.1f}% with 95% confidence,")
    print(f"  by the rule of three. It does not establish that the rate is")
    print(f"  zero.")

# =====================================================================
print()
print("=" * 96)
print(" SENSITIVITY PER COLLECTION")
print("=" * 96)
print(f"\n  The loss collections share loci, so pooling them counts the same")
print(f"  sample three times. Each column below is independent within itself.\n")

cols = sorted(targets.collection.unique())
print(f"  {'p <':<9}", end="")
for c_ in cols:
    print(f"{c_:>16}", end="")
print()
for th in THRESHOLDS:
    line = f"  {th:<9.3f}"
    for c_ in cols:
        g = targets[targets.collection == c_]
        if len(g):
            line += f"{int((g.pval < th).sum()):>8} of {len(g):<7}"
        else:
            line += f"{'-':>16}"
    print(line)

# =====================================================================
print()
print("=" * 96)
print(" WHERE THE P-VALUES SIT")
print("=" * 96)
print(f"\n  If the test is calibrated, control p-values should be spread")
print(f"  across the unit interval and target p-values pushed toward zero.\n")
print(f"  {'range':<14}{'targets':>10}{'controls':>11}")
for lo, hi in [(0, 0.01), (0.01, 0.05), (0.05, 0.1),
               (0.1, 0.25), (0.25, 0.5), (0.5, 1.01)]:
    nt = int(((targets.pval >= lo) & (targets.pval < hi)).sum())
    nc = int(((controls.pval >= lo) & (controls.pval < hi)).sum())
    print(f"  {f'{lo}-{hi}':<14}{nt:>10}{nc:>11}")

print(f"\n  median p-value    targets {targets.pval.median():.4f}, "
      f"controls {controls.pval.median():.4f}")

from scipy import stats
if len(targets) > 3 and len(controls) > 3:  # needs both to be scored
    u, p = stats.mannwhitneyu(targets.pval, controls.pval,
                              alternative="less")
    print(f"  targets below controls: Mann-Whitney p = {p:.4g}")
    auc = 1 - u / (len(targets) * len(controls))
    print(f"  AUC = {auc:.3f}")
    print(f"\n  The AUC is threshold-free: it is the probability that a")
    print(f"  randomly chosen target scores lower than a randomly chosen")
    print(f"  control. 0.5 would mean the test carries no information.")

R.to_csv(f"{RES}/lohhla_roc{SUFFIX}.tsv", sep="\t", index=False)
print(f"\nwritten to {RES}/lohhla_roc{SUFFIX}.tsv")
