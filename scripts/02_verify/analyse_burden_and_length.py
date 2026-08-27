#!/usr/bin/env python3
"""
Two questions the existing data can answer without recomputing anything.

First: does mutational burden affect recall? Samples here carry between 15
and 342 panel mutations, and if a caller performed worse on crowded
samples that would be a property of the pipeline worth knowing. The
expectation is that it does not - Mutect2 assembles haplotypes locally and
one mutation should not know about another 10 kb away - but the samples
also differ in VAF, so burden and fraction have to be separated rather
than compared raw.

Second: does peptide length change what the two predictors do? Steps 5 and
7 pooled 8- to 11-mers. The lengths are not equivalent: 9-mers dominate
the class I repertoire biologically, and a predictor calibrated mostly on
them may behave differently at the edges. If mhcflurry and NetMHCpan
disagree, the disagreement may sit at particular lengths.

Neither needs new computation. Both come out of tables already written.

Usage:
  python analyse_burden_and_length.py [workspace]
"""
import sys, os, glob
from collections import Counter
import numpy as np
import pandas as pd
from scipy import stats

args = [a for a in sys.argv[1:] if not a.startswith("--")]
WS = args[0] if args else os.environ.get(
    "WS", "/gpfs/bwfor/work/ws/fr_os136-immune_escape")

COHORT = f"{WS}/simulation/cohort"
PVAC = f"{WS}/simulation/pvacseq/cohort"
MANIFEST = f"{COHORT}/manifest.tsv"
OUT = os.path.expanduser("~/immune_escape_project/results")

STRONG = 0.5
os.makedirs(OUT, exist_ok=True)
man = pd.read_csv(MANIFEST, sep="\t")
pd.set_option("display.width", 220)

# =====================================================================
#  PART 1 - mutational burden against recall
# =====================================================================
print("=" * 96)
print(" DOES MUTATIONAL BURDEN AFFECT RECALL?")
print("=" * 96)

rows, per_mut = [], []
for _, m in man.iterrows():
    sid = m.sample_id
    comp = f"{COHORT}/{sid}/comparison/truth_vs_calls.tsv"
    if not os.path.exists(comp):
        continue
    t = pd.read_csv(comp, sep="\t")
    if not len(t):
        continue
    t["sample"] = sid
    t["cohort"] = m.cohort
    t["burden"] = len(t)
    per_mut.append(t)
    rows.append({
        "sample": sid, "cohort": m.cohort,
        "burden": len(t),
        "recall": 100 * t.detected_pass.mean(),
        "seen": 100 * t.detected_any.mean() if "detected_any" in t else None,
        "median_vaf": t.VAF.median(),
    })

S = pd.DataFrame(rows)
M = pd.concat(per_mut) if per_mut else pd.DataFrame()

if len(S) > 5:
    r_b, p_b = stats.pearsonr(S.burden, S.recall)
    r_v, p_v = stats.pearsonr(S.median_vaf, S.recall)
    rs_b, ps_b = stats.spearmanr(S.burden, S.recall)

    print(f"\n  {len(S)} samples, burden {S.burden.min():.0f} to "
          f"{S.burden.max():.0f} mutations\n")
    print(f"  recall vs burden      Pearson r = {r_b:+.3f}  (p = {p_b:.3g})")
    print(f"                        Spearman r = {rs_b:+.3f}  (p = {ps_b:.3g})")
    print(f"  recall vs median VAF  Pearson r = {r_v:+.3f}  (p = {p_v:.3g})")

    # burden and VAF are not independent - hypermutated donors also tend
    # to have lower fractions - so the raw correlation with burden could
    # be VAF wearing a different hat
    r_bv, p_bv = stats.pearsonr(S.burden, S.median_vaf)
    print(f"  burden vs median VAF  Pearson r = {r_bv:+.3f}  (p = {p_bv:.3g})")

    try:
        import statsmodels.api as sm
        X = sm.add_constant(S[["burden", "median_vaf"]])
        fit = sm.OLS(S.recall, X).fit()
        print(f"\n  both together, ordinary least squares:")
        print(f"    burden      coef {fit.params['burden']:+.4f}  "
              f"p = {fit.pvalues['burden']:.3g}")
        print(f"    median VAF  coef {fit.params['median_vaf']:+.2f}  "
              f"p = {fit.pvalues['median_vaf']:.3g}")
        print(f"    R-squared   {fit.rsquared:.3f}")
    except ImportError:
        # partial correlation by hand: regress each on VAF, correlate the
        # residuals
        b_res = S.burden - np.poly1d(np.polyfit(S.median_vaf, S.burden, 1))(S.median_vaf)
        r_res = S.recall - np.poly1d(np.polyfit(S.median_vaf, S.recall, 1))(S.median_vaf)
        rp, pp = stats.pearsonr(b_res, r_res)
        print(f"\n  burden vs recall, holding VAF constant:")
        print(f"    partial r = {rp:+.3f}  (p = {pp:.3g})")

    print(f"\n  by burden group:")
    S["group"] = pd.cut(S.burden, bins=[0, 20, 30, 60, 10000],
                        labels=["<20", "20-29", "30-59", "60+"])
    g = S.groupby("group", observed=False).agg(
        n=("sample", "size"), recall=("recall", "median"),
        vaf=("median_vaf", "median")).round(3)
    print(g.to_string())

# the same question at the level of individual mutations, where VAF can
# be held constant exactly rather than approximately
if len(M) > 50:
    print(f"\n  per mutation, within VAF bands:")
    M["vaf_bin"] = pd.cut(M.VAF, bins=[0, 0.10, 0.20, 0.30, 1.01],
                          labels=["<10%", "10-20%", "20-30%", ">30%"])
    M["burden_group"] = pd.cut(M.burden, bins=[0, 30, 10000],
                               labels=["low burden", "high burden"])
    print(f"\n  {'VAF':<10}{'low burden':>16}{'high burden':>16}{'p':>10}")
    for vb in ["<10%", "10-20%", "20-30%", ">30%"]:
        g = M[M.vaf_bin == vb]
        lo = g[g.burden_group == "low burden"]
        hi = g[g.burden_group == "high burden"]
        if len(lo) > 5 and len(hi) > 5:
            tab = [[int(lo.detected_pass.sum()), len(lo) - int(lo.detected_pass.sum())],
                   [int(hi.detected_pass.sum()), len(hi) - int(hi.detected_pass.sum())]]
            _, p = stats.fisher_exact(tab)
            print(f"  {vb:<10}{100*lo.detected_pass.mean():>15.1f}%"
                  f"{100*hi.detected_pass.mean():>15.1f}%{p:>10.3f}")

    print(f"\n  Comparing within a VAF band removes the confound: if burden")
    print(f"  mattered on its own, crowded samples would do worse here too.")

# =====================================================================
#  PART 2 - peptide length
# =====================================================================
print()
print("=" * 96)
print(" DOES PEPTIDE LENGTH CHANGE WHAT THE PREDICTORS DO?")
print("=" * 96)

mf_len = {}      # length -> list of ranks
pv_len = {}
mf_binders = Counter()
pv_binders = Counter()
mf_total = Counter()
pv_total = Counter()

print("\n  reading tables ...", flush=True)
for i, m in enumerate(man.itertuples()):
    sid = m.sample_id

    f = f"{COHORT}/{sid}/neoantigens/neoantigens_all.tsv"
    if os.path.exists(f):
        d = pd.read_csv(f, sep="\t", low_memory=False,
                        usecols=lambda c: c in {
                            "peptide", "mhcflurry_presentation_percentile"})
        d = d.rename(columns={"mhcflurry_presentation_percentile": "rank"})
        d["rank"] = pd.to_numeric(d["rank"], errors="coerce")
        d = d.dropna(subset=["rank", "peptide"])
        d["len"] = d.peptide.astype(str).str.len()
        for L, g in d.groupby("len"):
            mf_total[L] += len(g)
            mf_binders[L] += int((g["rank"] < STRONG).sum())
            mf_len.setdefault(L, []).append(g["rank"].values)

    g2 = glob.glob(f"{PVAC}/{sid}/pvacseq_out/MHC_Class_I/*.all_epitopes.tsv")
    if g2:
        d = pd.read_csv(g2[0], sep="\t", low_memory=False,
                        usecols=lambda c: c in {"MT Epitope Seq",
                                                "Best MT Percentile"})
        d = d.rename(columns={"MT Epitope Seq": "peptide",
                              "Best MT Percentile": "rank"})
        d["rank"] = pd.to_numeric(d["rank"], errors="coerce")
        d = d.dropna(subset=["rank", "peptide"])
        d["len"] = d.peptide.astype(str).str.len()
        for L, g in d.groupby("len"):
            pv_total[L] += len(g)
            pv_binders[L] += int((g["rank"] < STRONG).sum())
            pv_len.setdefault(L, []).append(g["rank"].values)

lengths = sorted(set(mf_total) | set(pv_total))

print(f"\n  {'length':<9}{'mhcflurry':>32}{'NetMHCpan':>32}")
print(f"  {'':<9}{'peptides':>12}{'binders':>10}{'rate':>10}"
      f"{'peptides':>12}{'binders':>10}{'rate':>10}")
for L in lengths:
    mt, mb = mf_total.get(L, 0), mf_binders.get(L, 0)
    pt, pb = pv_total.get(L, 0), pv_binders.get(L, 0)
    print(f"  {L}mer{'':<4}{mt:>12,}{mb:>10}"
          f"{100*mb/mt if mt else 0:>9.2f}%"
          f"{pt:>12,}{pb:>10}{100*pb/pt if pt else 0:>9.2f}%")

tot_m = sum(mf_binders.values())
tot_p = sum(pv_binders.values())
if tot_m and tot_p:
    print(f"\n  share of all strong binders:")
    print(f"  {'length':<9}{'mhcflurry':>13}{'NetMHCpan':>13}{'difference':>13}")
    for L in lengths:
        sm = 100 * mf_binders.get(L, 0) / tot_m
        sp = 100 * pv_binders.get(L, 0) / tot_p
        print(f"  {L}mer{'':<4}{sm:>12.1f}%{sp:>12.1f}%{sp-sm:>+12.1f}")

    print(f"\n  9-mers dominate the class I repertoire biologically, so a")
    print(f"  predictor weighting them heavily is behaving as expected.")
    print(f"  A large difference between the two would localise their")
    print(f"  disagreement to particular lengths.")

print(f"\n  rank distribution by length:")
print(f"  {'length':<9}{'mhcflurry median':>20}{'NetMHCpan median':>20}")
for L in lengths:
    mv = np.concatenate(mf_len[L]) if L in mf_len else np.array([])
    pv = np.concatenate(pv_len[L]) if L in pv_len else np.array([])
    print(f"  {L}mer{'':<4}"
          f"{np.median(mv) if len(mv) else 0:>19.1f}"
          f"{np.median(pv) if len(pv) else 0:>19.1f}")

# does the binder rate differ by length more than chance
if len(lengths) > 2:
    tab = [[mf_binders.get(L, 0), mf_total.get(L, 0) - mf_binders.get(L, 0)]
           for L in lengths if mf_total.get(L, 0)]
    if len(tab) > 2:
        chi2, p, _, _ = stats.chi2_contingency(tab)
        print(f"\n  mhcflurry, binder rate across lengths: "
              f"chi2 p = {p:.3g}")
    tab = [[pv_binders.get(L, 0), pv_total.get(L, 0) - pv_binders.get(L, 0)]
           for L in lengths if pv_total.get(L, 0)]
    if len(tab) > 2:
        chi2, p, _, _ = stats.chi2_contingency(tab)
        print(f"  NetMHCpan, binder rate across lengths: "
              f"chi2 p = {p:.3g}")

S.to_csv(f"{OUT}/burden_vs_recall.tsv", sep="\t", index=False)
pd.DataFrame({
    "length": lengths,
    "mf_peptides": [mf_total.get(L, 0) for L in lengths],
    "mf_binders": [mf_binders.get(L, 0) for L in lengths],
    "pv_peptides": [pv_total.get(L, 0) for L in lengths],
    "pv_binders": [pv_binders.get(L, 0) for L in lengths],
}).to_csv(f"{OUT}/binders_by_length.tsv", sep="\t", index=False)
print(f"\nwritten to {OUT}/")
