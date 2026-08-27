#!/usr/bin/env python3
"""
Step 5 verification: do the neoantigen predictions behave like predictions?

A binding predictor that has silently failed still returns numbers. What
distinguishes a real run from a broken one is the shape of the output:
percentile ranks should span their full range with most peptides scoring
poorly, binder counts should scale with the number of mutations and with
how many alleles were typed, and different HLA genotypes should produce
different peptides.

The checks:

  coverage      how many injected mutations reached prediction at all,
                and where the rest were lost
  distribution  are percentile ranks spread as a rank statistic must be,
                or piled at one value
  scaling       binders against mutation count, and against the number of
                alleles typed - the three samples with four alleles
                instead of six should predict proportionally fewer
  peptides      lengths within 8-11, and how much overlap there is
                between samples with different genotypes
  per allele    which alleles dominate the binder set, since one allele
                carrying everything would suggest the others were dropped

Nothing here validates the biology of the predictions. It establishes
that mhcflurry ran on the intended input and produced output with the
structure a working predictor produces.

Usage:
  python verify_step5_mhcflurry.py [workspace]
"""
import sys, os, glob
from collections import Counter
import numpy as np
import pandas as pd

args = [a for a in sys.argv[1:] if not a.startswith("--")]
WS = args[0] if args else os.environ.get(
    "WS", "/gpfs/bwfor/work/ws/fr_os136-immune_escape")

COHORT = f"{WS}/simulation/cohort"
MANIFEST = f"{COHORT}/manifest.tsv"

STRONG = 0.5
WEAK = 2.0

man = pd.read_csv(MANIFEST, sep="\t")
print(f"checking {len(man)} samples\n", flush=True)

rows, problems = [], []
skip_reasons = {}
all_peptides = {}
allele_hits = Counter()
rank_pool = []

for _, m in man.iterrows():
    sid = m.sample_id
    d = f"{COHORT}/{sid}/neoantigens"
    r = {"sample": sid, "cohort": m.cohort}

    per_mut = f"{d}/neoantigens_per_mutation.tsv"
    if not os.path.exists(per_mut):
        problems.append((sid, "missing", "no per-mutation table"))
        rows.append(r)
        continue

    t = pd.read_csv(per_mut, sep="\t")
    r["mutations_predicted"] = len(t)
    for col, key in [("peptides", "peptides"), ("strong", "strong"),
                     ("weak", "weak")]:
        if col in t.columns:
            r[key] = int(t[col].sum())

    if "strong" in t.columns:
        r["mut_with_binder"] = int((t.strong > 0).sum())
        r["pct_with_binder"] = round(100 * (t.strong > 0).mean(), 1)

    # how many injected mutations reached this step
    truth = f"{COHORT}/{sid}/truth_set.tsv"
    if os.path.exists(truth):
        tr = pd.read_csv(truth, sep="\t")
        n_snv = int((tr.Variant_Type == "SNP").sum()) if "Variant_Type" in tr else len(tr)
        n_mis = int((tr.Variant_Classification == "Missense_Mutation").sum()) \
            if "Variant_Classification" in tr else None
        r["truth_snv"] = n_snv
        r["truth_missense"] = n_mis
        if n_mis:
            r["reached_pct"] = round(100 * len(t) / n_mis, 1)
            if len(t) > n_mis:
                problems.append((sid, "count",
                                 f"{len(t)} predicted from {n_mis} missense"))

    # how many alleles were available
    hla = f"{COHORT}/{sid}/optitype/{sid}_result.tsv"
    if os.path.exists(hla):
        h = pd.read_csv(hla, sep="\t")
        if len(h):
            row = h.iloc[0]
            alleles = [str(row.get(c, "")) for c in
                       ["A1", "A2", "B1", "B2", "C1", "C2"]]
            alleles = [a for a in alleles if "*" in a]
            r["n_alleles"] = len(alleles)
            r["n_unique_alleles"] = len(set(alleles))

    # the full peptide table, if it survives
    # why mutations dropped out before prediction
    skip = f"{d}/skipped_mutations.tsv"
    if os.path.exists(skip):
        sk = pd.read_csv(skip, sep="\t")
        r["skipped"] = len(sk)
        if "reason" in sk.columns:
            for why in sk.reason:
                skip_reasons[str(why)] = skip_reasons.get(str(why), 0) + 1

    full = None
    for cand in [f"{d}/neoantigens_all.tsv", f"{d}/all_peptides.tsv"]:
        if os.path.exists(cand):
            full = pd.read_csv(cand, sep="\t", low_memory=False)
            break

    if full is not None and len(full):
        pep_col = next((c for c in full.columns
                        if c.lower() in ("peptide", "mt_epitope_seq", "sequence")), None)
        # presentation percentile is what defines a binder here, not the
        # affinity percentile: both are present and they disagree
        rank_col = ("mhcflurry_presentation_percentile"
                    if "mhcflurry_presentation_percentile" in full.columns
                    else next((c for c in full.columns
                               if "percentile" in c.lower()), None))
        allele_col = next((c for c in full.columns
                           if c.lower() in ("allele", "hla", "mhc")), None)

        if pep_col:
            peps = full[pep_col].dropna().astype(str)
            r["peptide_rows"] = len(peps)
            r["unique_peptides"] = peps.nunique()
            lens = peps.str.len()
            r["len_min"] = int(lens.min())
            r["len_max"] = int(lens.max())
            bad = int(((lens < 8) | (lens > 11)).sum())
            if bad:
                problems.append((sid, "length", f"{bad} peptides outside 8-11"))
            all_peptides[sid] = set(peps[
                pd.to_numeric(full[rank_col], errors="coerce") < STRONG
            ]) if rank_col else set()

        if rank_col:
            v = pd.to_numeric(full[rank_col], errors="coerce").dropna()
            if len(v):
                rank_pool.append(v.values)
                r["rank_median"] = round(float(v.median()), 2)
                r["rank_min"] = round(float(v.min()), 4)
                r["rank_distinct"] = int(v.nunique())
                if v.nunique() < 20:
                    problems.append((sid, "flat",
                                     f"only {v.nunique()} distinct rank values"))

        if allele_col and rank_col:
            v = pd.to_numeric(full[rank_col], errors="coerce")
            for a in full.loc[v < STRONG, allele_col].dropna():
                allele_hits[str(a)] += 1

    rows.append(r)
    print(f"  {sid}  {r.get('mutations_predicted', 0)} mutations, "
          f"{r.get('strong', 0)} strong binders")

t = pd.DataFrame(rows)
pd.set_option("display.width", 230)
pd.set_option("display.max_columns", 30)

# =====================================================================
print()
print("=" * 96)
print(" PER SAMPLE")
print("=" * 96)
cols = [c for c in ["sample", "cohort", "truth_missense", "mutations_predicted",
                    "reached_pct", "n_alleles", "peptides", "strong", "weak",
                    "mut_with_binder", "pct_with_binder"] if c in t.columns]
print(t[cols].to_string(index=False))

# =====================================================================
print()
print("=" * 96)
print(" OVERALL")
print("=" * 96)
for label, col in [("samples with predictions", "mutations_predicted"),
                   ("mutations predicted on", "mutations_predicted"),
                   ("peptides scored", "peptides"),
                   ("strong binders (<0.5%)", "strong"),
                   ("weak binders (0.5-2%)", "weak")]:
    if col not in t.columns:
        continue
    if label.startswith("samples"):
        print(f"\n  {label:<28}{int((t[col] > 0).sum())} of {len(t)}")
    else:
        print(f"  {label:<28}{int(t[col].sum()):>8,}")

if "strong" in t.columns and "mutations_predicted" in t.columns:
    ok = t[t.mutations_predicted > 0]
    per = ok.strong / ok.mutations_predicted
    print(f"\n  strong binders per mutation  median {per.median():.2f}, "
          f"range {per.min():.2f} – {per.max():.2f}")

if "pct_with_binder" in t.columns:
    v = t.pct_with_binder.dropna()
    print(f"  mutations yielding a binder  median {v.median():.0f}%, "
          f"range {v.min():.0f} – {v.max():.0f}%")

# =====================================================================
print()
print("=" * 96)
print(" DOES THE OUTPUT SCALE THE WAY IT SHOULD?")
print("=" * 96)
from scipy import stats

ok = t[(t.get("mutations_predicted", pd.Series(dtype=float)) > 0)]
if len(ok) > 3 and "strong" in ok.columns:
    r_, p_ = stats.pearsonr(ok.mutations_predicted, ok.strong)
    print(f"\n  strong binders vs mutation count")
    print(f"    r = {r_:.3f}, p = {p_:.3g}")
    print(f"    A predictor returning a constant would show no correlation.")

if "n_alleles" in t.columns and t.n_alleles.notna().any():
    print(f"\n  by number of alleles typed:")
    g = ok.groupby("n_alleles").agg(
        n=("sample", "size"),
        mutations=("mutations_predicted", "median"),
        strong=("strong", "median"))
    g["per_mutation"] = (g.strong / g.mutations).round(2)
    print(g.to_string())
    if len(g) > 1:
        print(f"\n    Samples typed at four alleles should yield fewer")
        print(f"    binders per mutation than those typed at six.")

# =====================================================================
print()
print("=" * 96)
print(" PERCENTILE RANK DISTRIBUTION")
print("=" * 96)
if rank_pool:
    v = np.concatenate(rank_pool)
    print(f"\n  {len(v):,} peptide-allele scores")
    print(f"\n  {'range':<14}{'n':>10}{'share':>9}")
    for lo, hi, lab in [(0, 0.5, "<0.5% strong"), (0.5, 2, "0.5-2% weak"),
                        (2, 10, "2-10%"), (10, 50, "10-50%"), (50, 101, ">50%")]:
        n = int(((v >= lo) & (v < hi)).sum())
        print(f"  {lab:<14}{n:>10,}{100*n/len(v):>8.1f}%")
    print(f"\n  median {np.median(v):.1f}, "
          f"distinct values {len(np.unique(v)):,}")
    print(f"\n  A rank statistic must span its range with most peptides")
    print(f"  scoring poorly. A pile-up at one value would mean the")
    print(f"  predictor returned a default.")
else:
    print(f"\n  no per-peptide table found — only the per-mutation summary")
    print(f"  is available, so the rank distribution cannot be checked")

# =====================================================================
if allele_hits:
    print()
    print("=" * 96)
    print(" WHICH ALLELES CARRY THE BINDERS")
    print("=" * 96)
    total = sum(allele_hits.values())
    print(f"\n  {total:,} strong binder calls across "
          f"{len(allele_hits)} alleles\n")
    for a, n in allele_hits.most_common(12):
        print(f"    {a:<14}{n:>6}  ({100*n/total:.1f}%)")
    top = allele_hits.most_common(1)[0]
    if top[1] / total > 0.4:
        print(f"\n  {top[0]} carries {100*top[1]/total:.0f}% of all binders,")
        print(f"  which would suggest the other alleles were not scored.")

# =====================================================================
if len(all_peptides) > 2:
    print()
    print("=" * 96)
    print(" PEPTIDE OVERLAP BETWEEN SAMPLES")
    print("=" * 96)
    ids = [s for s in all_peptides if all_peptides[s]]
    if len(ids) > 2:
        shares = []
        for i in range(len(ids)):
            for j in range(i + 1, len(ids)):
                a, b = all_peptides[ids[i]], all_peptides[ids[j]]
                if a and b:
                    shares.append(len(a & b) / len(a | b))
        print(f"\n  pairwise Jaccard, {len(shares)} pairs")
        print(f"    median {np.median(shares):.3f}, max {max(shares):.3f}")
        print(f"\n  Different mutations and different HLA genotypes should")
        print(f"  produce largely different peptide sets. High overlap")
        print(f"  would mean the genotype was not being used.")

# =====================================================================
print()
print("=" * 96)
print(f" PROBLEMS: {len(problems)}")
print("=" * 96)
if problems:
    by_kind = {}
    for sid, kind, msg in problems:
        by_kind.setdefault(kind, []).append((sid, msg))
    for kind in sorted(by_kind):
        print(f"\n  {kind} ({len(by_kind[kind])}):")
        for sid, msg in by_kind[kind][:12]:
            print(f"    {sid:<7} {msg}")
else:
    print("\n  none")

if skip_reasons:
    print()
    print("=" * 96)
    print(" WHY MUTATIONS WERE SKIPPED BEFORE PREDICTION")
    print("=" * 96)
    total = sum(skip_reasons.values())
    print(f"\n  {total} mutations across the cohort\n")
    for why, n in sorted(skip_reasons.items(), key=lambda x: -x[1])[:12]:
        print(f"    {n:>5}  {why[:70]}")

out = os.path.expanduser("~/immune_escape_project/results/verify_step5.tsv")
os.makedirs(os.path.dirname(out), exist_ok=True)
t.to_csv(out, sep="\t", index=False)
print(f"\nwritten to {out}")
