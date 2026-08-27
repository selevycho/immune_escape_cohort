#!/usr/bin/env python3
"""
Every number that appears on a slide, in one place.

Numbers drift between drafts: a denominator changes, a figure is rebuilt,
a slide keeps the old value. This prints the current answer for each one
so the deck can be checked against a single page.

Usage:
  python all_numbers.py [workspace]
"""
import os, sys, glob, gzip
import pandas as pd, numpy as np

WS = sys.argv[1] if len(sys.argv) > 1 else os.environ["WS"]
RES = os.path.expanduser("~/immune_escape_project/results")
SIM = f"{WS}/simulation"
MHC = ("chr6", 29600000, 33100000)
num = lambda s: pd.to_numeric(s, errors="coerce")

def head(n, t):
    print(f"\n{'─'*66}\n SLIDE {n} — {t}\n{'─'*66}")

# ---------------------------------------------------------------- panel
head("3-4", "the panel")
p = pd.read_csv(f"{SIM}/panel/panel.bed", sep="\t", header=None,
                names=["chrom","start","end","gene"])
p["len"] = p.end - p.start
inm = (p.chrom=="chr6") & (p.end>MHC[1]) & (p.start<MHC[2])
print(f"  genes                    {p.gene.nunique()}")
print(f"  intervals                {len(p)}")
print(f"  total                    {p['len'].sum()/1e6:.2f} Mb")
print(f"  MHC                      {p[inm]['len'].sum()/1e6:.2f} Mb"
      f"   ({100*p[inm]['len'].sum()/p['len'].sum():.1f}%)")

man = pd.read_csv(f"{SIM}/cohort/manifest.tsv", sep="\t")
print(f"  samples                  {len(man)}")
tr = pd.concat([pd.read_csv(f, sep="\t")
                for f in glob.glob(f"{SIM}/cohort/*/truth_set.tsv")])
print(f"  mutations in truth sets  {len(tr)}")
print(f"    substitutions          {int((tr.Variant_Type=='SNP').sum())}")
print(f"    indels                 {int((tr.Variant_Type!='SNP').sum())}")
g = tr.groupby("Hugo_Symbol").size()
per = tr.groupby(tr.index // 1).size()
print(f"  genes hit, median/sample "
      f"{tr.groupby('source_patient').Hugo_Symbol.nunique().median():.0f}")

# --------------------------------------------------------------- depth
head(6, "coverage")
print(f"  see cohort_qc.tsv — median 34.5x, 31-46x, 99.87% mapped")

# ------------------------------------------------------------ injection
head("15-17", "what landed")
v = pd.read_csv(f"{RES}/verify_step2_per_mutation.tsv", sep="\t")
i = pd.read_csv(f"{RES}/verify_step8_per_indel.tsv", sep="\t")
n_sub, n_ok = len(v), int(v.landed.sum())
print(f"  substitutions attempted  {n_sub}")
print(f"  landed                   {n_ok}   ({100*n_ok/n_sub:.1f}%)")
print(f"  failed                   {n_sub-n_ok}")
if "in_mhc" in v.columns:
    f_mhc = int(v[~v.landed].in_mhc.sum())
    n_mhc = int(v.in_mhc.sum())
    print(f"    of those, in MHC       {f_mhc}"
          f"   ({100*f_mhc/(n_sub-n_ok):.0f}% of losses)")
    print(f"  MHC positions            {n_mhc}"
          f"   ({100*n_mhc/n_sub:.1f}% of all)")
    out = v[~v.in_mhc]
    print(f"  outside the MHC          {int(out.landed.sum())} of {len(out)}"
          f"   ({100*out.landed.mean():.1f}%)")
print(f"  indels attempted         {len(i)}")
print(f"  landed                   {int(i.landed.sum()) if 'landed' in i else len(i)}")
if "observed_vaf" in v.columns and "target_vaf" in v.columns:
    ok = v[v.landed & (v.target_vaf>0)]
    print(f"  fidelity, substitutions  "
          f"{(ok.observed_vaf/ok.target_vaf).median():.2f}")
    oi = i[i.observed_vaf>0] if "observed_vaf" in i else None
    if oi is not None and len(oi):
        print(f"  fidelity, indels         "
              f"{(oi.observed_vaf/oi.target_vaf).median():.2f}")

# ----------------------------------------------------------- optitype
head("19-20", "HLA typing")
a = pd.read_csv(f"{RES}/verify_step4_accuracy.tsv", sep="\t")
if "in_reference" in a.columns:
    a = a[a.in_reference]
tot = ok = 0
for loc in "ABC":
    if f"{loc}_correct" in a.columns:
        k, n = int(a[f"{loc}_correct"].sum()), int(a[f"{loc}_n_truth"].sum())
        ok += k; tot += n
        print(f"  HLA-{loc}                    {k} of {n}   ({100*k/n:.1f}%)")
print(f"  overall                  {ok} of {tot}   ({100*ok/tot:.1f}%)")

# ------------------------------------------------------------- mutect2
head("21-24", "variant calling")
d = pd.read_csv(f"{RES}/step3_rescored_full.tsv", sep="\t")
print(f"  verified in the reads    {len(d)}")
print(f"  emitted                  {int(d.emitted.sum())}"
      f"   ({100*d.emitted.mean():.1f}%)")
print(f"  PASS                     {int(d['pass'].sum())}"
      f"   ({100*d['pass'].mean():.1f}%)")
miss = d[~d['pass']]
print(f"  missed                   {len(miss)}")
print(f"    called then filtered   {int(miss.emitted.sum())}")
print(f"    never called           {len(miss)-int(miss.emitted.sum())}")
m = d[d.in_mhc]
print(f"  MHC                      {int(m['pass'].sum())} of {len(m)}")
print(f"  outside the MHC          {100*d[~d.in_mhc]['pass'].mean():.1f}%")

ind = pd.read_csv(f"{RES}/indel_recall.tsv", sep="\t")
print(f"  indels recovered         {int(ind.found.sum())} of {len(ind)}"
      f"   ({100*ind.found.mean():.1f}%)")

c = pd.read_csv(f"{RES}/step3_controls.tsv", sep="\t")
print(f"  PASS calls, all samples  {int(c.pass_total.sum())}")
print(f"  outside the truth set    {int(c.outside_truth.sum())}")

NC = f"{SIM}/step17b_split_control"
tot_f = mhc_f = n_f = 0
for f in glob.glob(f"{NC}/*/*.filtered.vcf.gz"):
    n_f += 1
    with gzip.open(f,'rt') as fh:
        for l in fh:
            if l.startswith('#'): continue
            x = l.split('\t')
            if len(x)<8 or x[6]!='PASS': continue
            tot_f += 1
            if x[0]=='chr6' and MHC[1]<=int(x[1])<MHC[2]: mhc_f += 1
if n_f:
    print(f"  negative control         {tot_f} false, {tot_f/n_f:.0f}/sample")
    print(f"    in the MHC             {mhc_f}   ({100*mhc_f/tot_f:.0f}%)")

# --------------------------------------------------------- predictors
head("25-28", "neoantigens")
mf = pd.read_csv(f"{RES}/mhcflurry_summary.tsv", sep="\t")
print(f"  missense in truth sets   "
      f"{int((tr.Variant_Classification=='Missense_Mutation').sum())}")
print(f"  used for prediction      {int(mf.mutations.sum())}")
print(f"  peptides scored          {int(mf.peptides.sum()):,}")
print(f"  strong binders           {int(mf.strong.sum())}")
print(f"  weak binders             {int(mf.weak.sum())}")
print(f"  mutations with a binder  {int(mf.with_binder.sum())}"
      f"   ({100*mf.with_binder.sum()/mf.mutations.sum():.1f}%)")
for n, gg in mf.groupby("n_alleles"):
    print(f"    {int(n)} alleles            "
          f"{(gg.strong/gg.mutations).median():.2f} per mutation")

pv = pd.read_csv(f"{RES}/pvacseq_summary.tsv", sep="\t")
print(f"  pVACseq epitopes         {int(pv.epitopes.sum()):,}")
print(f"  pVACseq strong           {int(pv.strong.sum())}")
j = pv.merge(mf[["sample","strong"]], on="sample", suffixes=("_pv","_mf"))
print(f"  correlation              "
      f"r = {np.corrcoef(j.strong_mf, j.strong_pv)[0,1]:.3f}")
print(f"  ratio                    "
      f"{j.strong_pv.sum()/j.strong_mf.sum():.2f}")

# -------------------------------------------------------------- lohhla
head("29-31", "HLA loss")
r = pd.read_csv(f"{RES}/lohhla_report.tsv", sep="\t")
r["pval"] = num(r.pval)
s9 = pd.read_csv(f"{SIM}/step9_lohhla_panel/results/step9_all_loci_mc5.tsv",
                 sep="\t")
s9["pval"] = num(s9.pval)
print(f"  locus attempts           {len(s9)}")
print(f"  reached the test         {int(s9.pval.notna().sum())}")
print(f"  significant p<0.05       {int((s9.pval<0.05).sum())}")
sg = s9[s9.pval<0.05]
if "role" in s9.columns:
    print(f"    on target              {int((sg.role=='TARGET').sum())}")
    print(f"    on controls            {int((sg.role!='TARGET').sum())}")
print(f"\n  ATTENTION — check what the collection names mean:")
for c, gg in s9.groupby("collection"):
    k = gg[gg.pval<0.05]
    print(f"    {c:<14}{len(gg):>5} loci, {int(gg.pval.notna().sum()):>3} tested,"
          f" {len(k):>2} significant")
if "keep" in s9.columns:
    print(f"\n  keep fraction per collection:")
    for c, gg in s9.groupby("collection"):
        print(f"    {c:<14}keep = {gg.keep.iloc[0]}")

roc = pd.read_csv(f"{SIM}/step9_lohhla_panel/results/lohhla_roc_mc5.tsv",
                  sep="\t")
print(f"\n  sensitivity at p<0.05    {roc[roc.threshold==0.05].sensitivity.iloc[0]:.1f}%")
print(f"  false positives, any p   {int(roc.fp.sum())}")

# -------------------------------------------------------------- sweeps
head("32-34", "sweeps")
f = pd.read_csv(f"{SIM}/step10_filter_sweep/results/filter_sweep_summary.tsv",
                sep="\t")
for _, x in f.iterrows():
    print(f"  {x.setting:<16}{int(x.recovered):>6} recovered,"
          f" {int(x.not_in_truth)} false")
print(f"  NOTE: those recalls use 1 528 as denominator, the slides use 1 497")

cv = pd.read_csv(f"{RES}/coverage_sweep_summary.tsv", sep="\t")
print()
for _, x in cv.sort_values("depth").iterrows():
    print(f"  {int(x.depth)}x{'':<14}{int(x.found):>6} of {int(x.verified)}"
          f"   {x.recall:.1f}%")

gg = d.groupby("sample").agg(n=("pass","size"), rec=("pass","mean"),
                             vaf=("observed_vaf","median"))
print(f"\n  recall vs count          r = {np.corrcoef(gg.n, gg.rec)[0,1]:.2f}")
print(f"  recall vs VAF            r = {np.corrcoef(gg.vaf, gg.rec)[0,1]:.2f}")
print()
