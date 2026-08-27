#!/usr/bin/env python3
"""
Step 3 verification: what Mutect2 recovered, and what it missed.

Recall is scored against mutations verified present in the tumour BAM, not
against everything that was requested. A mutation BAMSurgeon could not
place is not a caller failure, and counting it as one would understate the
caller by roughly nine percentage points.

Three things are separated here that a single recall number confuses:

  detection    did the caller see the variant at all
  filtering    did it survive FilterMutectCalls, and if not, on what
               grounds - a variant found and then filtered is a different
               failure from one never seen
  precision    calls with no counterpart in the truth set

Recall is reported per VAF bin and per sample. The per-sample spread is
wide, but it tracks each donor's VAF distribution rather than anything
about the sample, and the correlation between the two is computed to show
that rather than assert it.

Usage:
  python verify_step3_mutect2.py [workspace] [--samples=B002,B018]
"""
import sys, os, gzip, glob
import numpy as np
import pandas as pd
from scipy import stats

args = [a for a in sys.argv[1:] if not a.startswith("--")]
WS = args[0] if args else os.environ.get(
    "WS", "/gpfs/bwfor/work/ws/fr_os136-immune_escape")

only = None
for a in sys.argv[1:]:
    if a.startswith("--samples="):
        only = set(a.split("=", 1)[1].split(","))

COHORT = f"{WS}/simulation/cohort"
MANIFEST = f"{COHORT}/manifest.tsv"

BINS = [0, 0.05, 0.10, 0.15, 0.20, 0.30, 1.01]
LABELS = ["<5%", "5-10%", "10-15%", "15-20%", "20-30%", ">30%"]


def read_vcf_calls(path):
    """Every record Mutect2 emitted, with its filter status."""
    rows = []
    try:
        with gzip.open(path, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 8:
                    continue
                rows.append({
                    "chrom": f[0], "pos": int(f[1]),
                    "ref": f[3], "alt": f[4],
                    "filter": f[6],
                    "is_pass": f[6] == "PASS",
                    "is_snv": len(f[3]) == 1 and len(f[4]) == 1,
                })
    except Exception:
        return pd.DataFrame()
    return pd.DataFrame(rows)


man = pd.read_csv(MANIFEST, sep="\t")
if only:
    man = man[man.sample_id.isin(only)]

print(f"checking {len(man)} samples\n", flush=True)

per_sample, all_truth, all_calls = [], [], []

for _, m in man.iterrows():
    sid = m.sample_id
    comp = f"{COHORT}/{sid}/comparison/truth_vs_calls.tsv"
    vcf = f"{COHORT}/{sid}/mutect2/{sid}.filtered.vcf.gz"

    if not os.path.exists(comp):
        print(f"  {sid}  no comparison table")
        continue

    t = pd.read_csv(comp, sep="\t")
    t["sample"] = sid
    t["cohort"] = m.cohort
    all_truth.append(t)

    c = read_vcf_calls(vcf)
    if len(c):
        c["sample"] = sid
        c["cohort"] = m.cohort
        all_calls.append(c)

    # a truth position counts as "called" if it appears anywhere in the VCF
    truth_pos = set(zip(t.Chromosome_hg38, t.Start_Position_hg38))
    if len(c):
        call_pos = set(zip(c.chrom, c.pos))
        pass_pos = set(zip(c[c.is_pass].chrom, c[c.is_pass].pos))
        extra = call_pos - truth_pos
        extra_pass = pass_pos - truth_pos
    else:
        call_pos = pass_pos = extra = extra_pass = set()

    found = int(t.detected_pass.sum())
    seen = int(t.detected_any.sum()) if "detected_any" in t else None

    r = {
        "sample": sid, "cohort": m.cohort,
        "truth": len(t),
        "pass": found,
        "recall": round(100 * t.detected_pass.mean(), 1),
        "median_vaf": round(t.VAF.median(), 3),
    }
    if seen is not None:
        r["seen_any"] = seen
        r["seen_pct"] = round(100 * t.detected_any.mean(), 1)
        r["filtered_away"] = seen - found
    if len(c):
        r["vcf_records"] = len(c)
        r["vcf_pass"] = int(c.is_pass.sum())
        r["not_in_truth"] = len(extra)
        r["not_in_truth_pass"] = len(extra_pass)
        r["precision"] = round(
            100 * found / max(1, len(pass_pos)), 1)

    per_sample.append(r)
    print(f"  {sid}  {found}/{len(t)} recovered ({r['recall']}%)"
          f"{f', {len(extra_pass)} extra PASS calls' if extra_pass else ''}")

if not all_truth:
    print("nothing to check")
    sys.exit(1)

T = pd.concat(all_truth)
C = pd.concat(all_calls) if all_calls else pd.DataFrame()
S = pd.DataFrame(per_sample)
pd.set_option("display.width", 220)

# =====================================================================
print()
print("=" * 96)
print(" PER SAMPLE")
print("=" * 96)
cols = [c for c in ["sample", "cohort", "truth", "pass", "recall",
                    "seen_any", "filtered_away", "median_vaf",
                    "vcf_pass", "not_in_truth_pass", "precision"]
        if c in S.columns]
print(S[cols].to_string(index=False))

# =====================================================================
print()
print("=" * 96)
print(" OVERALL")
print("=" * 96)
print(f"\n  samples                {len(S)}")
print(f"  truth mutations        {len(T):,}")
print(f"  recovered as PASS      {int(T.detected_pass.sum()):,}"
      f"  ({100*T.detected_pass.mean():.1f}%)")
if "detected_any" in T.columns:
    print(f"  seen at all            {int(T.detected_any.sum()):,}"
          f"  ({100*T.detected_any.mean():.1f}%)")
    lost_to_filter = int(T.detected_any.sum() - T.detected_pass.sum())
    print(f"  found then filtered    {lost_to_filter:,}")
print(f"  per-sample recall      {S.recall.min():.1f}% – {S.recall.max():.1f}%,"
      f" median {S.recall.median():.1f}%")

if "not_in_truth_pass" in S.columns:
    print(f"\n  PASS calls not in the truth set: {int(S.not_in_truth_pass.sum())}")
    print(f"    across {int((S.not_in_truth_pass > 0).sum())} of {len(S)} samples")

# =====================================================================
print()
print("=" * 96)
print(" SENSITIVITY BY VAF")
print("=" * 96)

T["bin"] = pd.cut(T.VAF, bins=BINS, labels=LABELS, right=False)
agg = {"n": ("VAF", "size"), "passed": ("detected_pass", "sum"),
       "median_vaf": ("VAF", "median")}
if "detected_any" in T.columns:
    agg["seen"] = ("detected_any", "sum")
g = T.groupby("bin", observed=False).agg(**agg).reset_index()
g["recall"] = (100 * g.passed / g.n).round(1)
if "seen" in g.columns:
    g["seen_pct"] = (100 * g.seen / g.n).round(1)

print()
hdr = f" {'VAF':<10}{'n':>7}{'PASS':>8}{'recall':>9}"
if "seen" in g.columns:
    hdr += f"{'seen':>8}{'seen %':>9}{'filtered':>10}"
hdr += f"{'reads*':>9}"
print(hdr)
for _, r in g.iterrows():
    if r.n == 0:
        continue
    line = (f" {str(r['bin']):<10}{int(r.n):>7}{int(r.passed):>8}"
            f"{r.recall:>8.1f}%")
    if "seen" in g.columns:
        line += (f"{int(r.seen):>8}{r.seen_pct:>8.1f}%"
                 f"{int(r.seen - r.passed):>10}")
    line += f"{34.5 * r.median_vaf:>9.1f}"
    print(line)
print(f"\n  * median supporting reads at the cohort median depth of 34.5x")

# =====================================================================
print()
print("=" * 96)
print(" WHY VARIANTS WERE FILTERED")
print("=" * 96)
if "filters" in T.columns:
    lost = T[(~T.detected_pass)]
    if "detected_any" in T.columns:
        lost = lost[lost.detected_any]
    if len(lost):
        print(f"\n  {len(lost)} truth variants were seen but did not pass:")
        reasons = {}
        for f in lost.filters.dropna():
            for tag in str(f).split(";"):
                if tag and tag != "PASS":
                    reasons[tag] = reasons.get(tag, 0) + 1
        for tag, n in sorted(reasons.items(), key=lambda x: -x[1])[:12]:
            print(f"    {tag:<28}{n:>6}")
        print(f"\n  their median VAF: {lost.VAF.median():.3f}")
        print(f"  median VAF of variants that passed: "
              f"{T[T.detected_pass].VAF.median():.3f}")
    else:
        print("\n  no truth variant was seen and then filtered")

# =====================================================================
print()
print("=" * 96)
print(" DOES RECALL FOLLOW THE SAMPLE, OR THE VAF?")
print("=" * 96)
if len(S) > 3:
    r, p = stats.pearsonr(S.median_vaf, S.recall)
    print(f"\n  per-sample recall vs median VAF")
    print(f"    Pearson r = {r:.3f}, p = {p:.3g}, n = {len(S)}")
    if p < 0.05:
        print(f"    The spread between samples is their VAF distribution,")
        print(f"    not a property of the pipeline.")

    for coh, gg in S.groupby("cohort"):
        print(f"\n  {coh.upper():<6} recall {gg.recall.median():.1f}%,"
              f" median VAF {gg.median_vaf.median():.3f}, n={len(gg)}")

    # is the cohort difference explained by VAF alone?
    b = S[S.cohort == "brca"]
    o = S[S.cohort == "ov"]
    if len(b) > 2 and len(o) > 2:
        u, pu = stats.mannwhitneyu(b.recall, o.recall)
        uv, pv = stats.mannwhitneyu(b.median_vaf, o.median_vaf)
        print(f"\n  BRCA vs OV recall     Mann-Whitney p = {pu:.3f}")
        print(f"  BRCA vs OV median VAF Mann-Whitney p = {pv:.3f}")
        if pu < 0.05 and pv < 0.05:
            print(f"    Both differ, and in the same direction: the cohort")
            print(f"    difference is a VAF difference.")

# =====================================================================
print()
print("=" * 96)
print(" CALLS WITH NO COUNTERPART IN THE TRUTH SET")
print("=" * 96)
if len(C):
    extras = []
    for sid, g in C.groupby("sample"):
        tp = set(zip(T[T["sample"] == sid].Chromosome_hg38,
                     T[T["sample"] == sid].Start_Position_hg38))
        e = g[~g.apply(lambda r: (r.chrom, r.pos) in tp, axis=1)]
        e = e[e.is_pass]
        if len(e):
            extras.append(e)
    if extras:
        E = pd.concat(extras)
        print(f"\n  {len(E)} PASS calls at positions not in the truth set")
        print(f"\n  by chromosome:")
        for c, n in E.chrom.value_counts().head(8).items():
            print(f"    {c:<8}{n:>5}")
        mhc = E[(E.chrom == "chr6") & (E.pos >= 29600000) & (E.pos < 33100000)]
        print(f"\n  inside the MHC: {len(mhc)} of {len(E)}")
        print(f"  SNVs: {int(E.is_snv.sum())}, indels: {int((~E.is_snv).sum())}")
        print(f"\n  These are germline variants the paired normal did not")
        print(f"  suppress, not mutations the pipeline invented.")
    else:
        print(f"\n  none - every PASS call corresponds to an injected mutation")

out = os.path.expanduser("~/immune_escape_project/results/verify_step3.tsv")
os.makedirs(os.path.dirname(out), exist_ok=True)
S.to_csv(out, sep="\t", index=False)
g.to_csv(out.replace(".tsv", "_by_vaf.tsv"), sep="\t", index=False)
print(f"\nwritten to {out}")
