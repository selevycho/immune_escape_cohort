#!/usr/bin/env python3
"""
Step 10 analysis - what loosening the Mutect2 filters actually buys.

Step 3 found that Mutect2 detects 91.7% of the verified mutations but
reports only 73.0% as PASS: 286 variants were seen and then removed, at a
median VAF of 0.098 against 0.242 for those that survived. Precision was
perfect - no PASS call in forty samples fell outside the truth set - and a
caller with no false positives has room to be asked for more.

Five filter settings were applied to the same unfiltered VCFs. Since
Mutect2 itself was not re-run, every setting sees exactly the same
underlying calls and the differences are the filters alone.

What matters is the trade. Recall rises as the filters loosen; the
question is whether precision falls with it, and by how much. A setting
that adds twenty true variants and one false one is worth taking; one
that adds twenty of each is not.

Note on what counts as a false positive. Any PASS call at a position not
in the truth set is counted as one, which overstates the rate: the panel
also contains real germline variants the paired normal did not fully
suppress, and those are indistinguishable here from a caller error. The
number is therefore an upper bound.

Usage:
  python analyse_step10_filters.py [workspace]
"""
import sys, os, gzip
import numpy as np
import pandas as pd

args = [a for a in sys.argv[1:] if not a.startswith("--")]
WS = args[0] if args else os.environ.get(
    "WS", "/gpfs/bwfor/work/ws/fr_os136-immune_escape")

COHORT = f"{WS}/simulation/cohort"
SWEEP = f"{WS}/simulation/step10_filter_sweep"
MANIFEST = f"{COHORT}/manifest.tsv"
SETTINGS = f"{SWEEP}/settings.tsv"
OUT = f"{SWEEP}/results"

BINS = [0, 0.05, 0.10, 0.15, 0.20, 0.30, 1.01]
LABELS = ["<5%", "5-10%", "10-15%", "15-20%", "20-30%", ">30%"]

os.makedirs(OUT, exist_ok=True)


def read_calls(path):
    """Every record in a VCF, with its filter status."""
    out = {}
    try:
        with gzip.open(path, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 8:
                    continue
                out[(f[0], int(f[1]))] = {
                    "filter": f[6],
                    "pass": f[6] == "PASS",
                    "is_snv": len(f[3]) == 1 and len(f[4]) == 1,
                }
    except Exception:
        pass
    return out


settings = pd.read_csv(SETTINGS, sep="\t")
order = settings.setting.tolist()
desc = dict(zip(settings.setting, settings.description))

man = pd.read_csv(MANIFEST, sep="\t")
print(f"{len(order)} settings, {len(man)} samples\n", flush=True)

rows, per_mut = [], []

for setting in order:
    print(f"  {setting} ...", end="", flush=True)
    n_done = 0

    for _, m in man.iterrows():
        sid = m.sample_id
        vcf = f"{SWEEP}/{setting}/{sid}/{sid}.filtered.vcf.gz"
        truth_p = f"{COHORT}/{sid}/truth_set.tsv"

        if not (os.path.exists(vcf) and os.path.exists(truth_p)):
            continue

        t = pd.read_csv(truth_p, sep="\t")
        if "Variant_Type" in t.columns:
            # step 8 injects the indels into a separate BAM, so scoring
            # them here would report failures for mutations this VCF was
            # never asked about
            t = t[t.Variant_Type == "SNP"]
        if not len(t):
            continue

        calls = read_calls(vcf)
        truth_pos = set(zip(t.Chromosome_hg38, t.Start_Position_hg38))

        found = seen = 0
        for _, r in t.iterrows():
            key = (r.Chromosome_hg38, int(r.Start_Position_hg38))
            c = calls.get(key)
            if c:
                seen += 1
                if c["pass"]:
                    found += 1
            per_mut.append({
                "setting": setting, "sample": sid, "cohort": m.cohort,
                "chrom": key[0], "pos": key[1], "VAF": r.get("VAF"),
                "seen": c is not None,
                "passed": bool(c and c["pass"]),
                "filters": c["filter"] if c else "",
            })

        extra = [k for k, v in calls.items()
                 if v["pass"] and k not in truth_pos]

        rows.append({
            "setting": setting, "sample": sid, "cohort": m.cohort,
            "truth": len(t),
            "seen": seen,
            "passed": found,
            "recall": round(100 * found / len(t), 2),
            "seen_pct": round(100 * seen / len(t), 2),
            "filtered_away": seen - found,
            "vcf_pass": sum(1 for v in calls.values() if v["pass"]),
            "not_in_truth": len(extra),
            "median_vaf": round(t.VAF.median(), 3) if "VAF" in t else None,
        })
        n_done += 1

    print(f" {n_done} samples")

S = pd.DataFrame(rows)
M = pd.DataFrame(per_mut)

if S.empty:
    print("nothing to analyse")
    sys.exit(1)

pd.set_option("display.width", 230)

# =====================================================================
print()
print("=" * 100)
print(" RECALL AND PRECISION BY SETTING")
print("=" * 100)

summary = []
for setting in order:
    g = S[S.setting == setting]
    if not g.empty:
        n_truth = int(g.truth.sum())
        n_pass = int(g.passed.sum())
        n_seen = int(g.seen.sum())
        n_extra = int(g.not_in_truth.sum())
        n_calls = int(g.vcf_pass.sum())
        summary.append({
            "setting": setting,
            "samples": len(g),
            "truth": n_truth,
            "recovered": n_pass,
            "recall": 100 * n_pass / n_truth,
            "seen": 100 * n_seen / n_truth,
            "filtered": n_seen - n_pass,
            "pass_calls": n_calls,
            "not_in_truth": n_extra,
            "precision": 100 * n_pass / max(1, n_calls),
        })

D = pd.DataFrame(summary)
print(f"\n  {'setting':<16}{'recall':>9}{'seen':>8}{'filtered':>10}"
      f"{'recovered':>11}{'extra calls':>13}{'precision':>11}")
for _, r in D.iterrows():
    print(f"  {r.setting:<16}{r.recall:>8.1f}%{r.seen:>7.1f}%"
          f"{int(r.filtered):>10}{int(r.recovered):>11}"
          f"{int(r.not_in_truth):>13}{r.precision:>10.1f}%")

base = D[D.setting == "default"]
if len(base):
    b = base.iloc[0]
    print(f"\n  relative to default:")
    print(f"  {'setting':<16}{'extra true':>12}{'extra false':>13}"
          f"{'true per false':>16}")
    for _, r in D.iterrows():
        if r.setting == "default":
            continue
        dt = int(r.recovered - b.recovered)
        df_ = int(r.not_in_truth - b.not_in_truth)
        ratio = f"{dt/df_:.1f}" if df_ > 0 else ("all true" if dt > 0 else "-")
        print(f"  {r.setting:<16}{dt:>+12}{df_:>+13}{ratio:>16}")

    print(f"\n  A setting is worth taking if it adds true variants faster")
    print(f"  than false ones. The last column is that rate.")

# =====================================================================
print()
print("=" * 100)
print(" WHERE THE GAIN COMES FROM")
print("=" * 100)

M["bin"] = pd.cut(M.VAF, bins=BINS, labels=LABELS, right=False)
print(f"\n  recall by VAF bin\n")
header = f"  {'VAF':<10}{'n':>7}"
for s in order:
    header += f"{s[:9]:>11}"
print(header)

for lab in LABELS:
    g = M[M.bin == lab]
    if g.empty:
        continue
    n = len(g[g.setting == order[0]])
    if n == 0:
        continue
    line = f"  {lab:<10}{n:>7}"
    for s in order:
        gs = g[g.setting == s]
        line += f"{100*gs.passed.mean() if len(gs) else 0:>10.1f}%"
    print(line)

print(f"\n  Loosening the filters helps where the evidence was thin and")
print(f"  changes nothing where it was already strong.")

# =====================================================================
print()
print("=" * 100)
print(" WHICH FILTERS STOP FIRING")
print("=" * 100)

for setting in order:
    g = M[(M.setting == setting) & (M.seen) & (~M.passed)]
    if g.empty:
        print(f"\n  {setting}: nothing seen and then filtered")
        continue
    reasons = {}
    for f in g.filters.dropna():
        for tag in str(f).split(";"):
            if tag and tag != "PASS":
                reasons[tag] = reasons.get(tag, 0) + 1
    top = sorted(reasons.items(), key=lambda x: -x[1])[:5]
    print(f"\n  {setting}  ({len(g)} seen but filtered)")
    for tag, n in top:
        print(f"    {tag:<28}{n:>6}")

# =====================================================================
print()
print("=" * 100)
print(" PER SAMPLE, RECALL BY SETTING")
print("=" * 100)
piv = S.pivot_table(index="sample", columns="setting",
                    values="recall", aggfunc="first")
piv = piv[[c for c in order if c in piv.columns]]
print()
print(piv.round(1).to_string())

print(f"\n  spread between the strictest and loosest setting:")
if "strict" in piv.columns and order[-1] in piv.columns:
    loosest = "verypermissive" if "verypermissive" in piv.columns else order[-2]
    delta = piv[loosest] - piv["strict"]
    print(f"    median {delta.median():.1f} points, "
          f"range {delta.min():.1f} to {delta.max():.1f}")

# =====================================================================
print()
print("=" * 100)
print(" WHAT TO CONCLUDE")
print("=" * 100)
if len(D) > 1 and len(base):
    b = base.iloc[0]
    best = D.loc[D.recall.idxmax()]
    print(f"\n  default recovers {int(b.recovered)} of {int(b.truth)} "
          f"({b.recall:.1f}%) with {int(b.not_in_truth)} calls outside the truth set")
    print(f"  {best.setting} recovers {int(best.recovered)} "
          f"({best.recall:.1f}%) with {int(best.not_in_truth)}")
    gain = best.recall - b.recall
    print(f"\n  The most permissive setting tested moves recall by "
          f"{gain:+.1f} points.")
    print(f"  Mutect2 sees {b.seen:.1f}% of the truth set at every setting -")
    print(f"  the ceiling is detection, not filtering, and no filter")
    print(f"  configuration reaches past it.")

S.to_csv(f"{OUT}/filter_sweep_per_sample.tsv", sep="\t", index=False)
D.to_csv(f"{OUT}/filter_sweep_summary.tsv", sep="\t", index=False)
M.to_csv(f"{OUT}/filter_sweep_per_mutation.tsv", sep="\t", index=False)
print(f"\nwritten to {OUT}/")
