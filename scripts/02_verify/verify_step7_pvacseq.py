#!/usr/bin/env python3
"""
Step 7 verification: the second route to the same answer.

Step 5 took mutations from the lifted MAF with the protein consequence
already annotated by TCGA, cut peptides out of GENCODE sequences, and
scored them with mhcflurry. This step starts from the Mutect2 VCF, lets
VEP recompute the consequence from the genome, and lets pVACseq build the
peptides before handing them to NetMHCpan. The two share no code and
begin from different files.

That makes the comparison between them the point of the exercise. Where
they agree, neither is likely to carry a systematic error; where they
diverge, the divergence localises to a step only one of them performs.

Note on what counts as a binder. pVACseq writes a filtered.tsv that is
routinely empty - its filter chain shortlists vaccine candidates and cut
49 690 epitopes to three on the first sample tested. Binders are counted
here from all_epitopes by percentile rank, the same threshold applied to
mhcflurry, so the two routes are measured the same way.

Usage:
  python verify_step7_pvacseq.py [workspace]
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

STRONG = 0.5
WEAK = 2.0
RANK_COL = "Best MT Percentile"

man = pd.read_csv(MANIFEST, sep="\t")
print(f"checking {len(man)} samples\n", flush=True)

rows, problems = [], []
rank_pool, allele_hits = [], Counter()
gene_sets, peptide_sets = {}, {}
length_counts = Counter()

for _, m in man.iterrows():
    sid = m.sample_id
    d = f"{PVAC}/{sid}"
    r = {"sample": sid, "cohort": m.cohort}

    # ---------------- VEP ----------------
    vep = f"{d}/{sid}.vep.vcf"
    if os.path.exists(vep):
        n_rec = n_csq = 0
        with open(vep, errors="ignore") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                n_rec += 1
                if "CSQ=" in line:
                    n_csq += 1
        r["vep_records"] = n_rec
        r["vep_annotated"] = n_csq
        if n_rec and n_csq < n_rec:
            problems.append((sid, "vep",
                             f"{n_rec - n_csq} records without CSQ"))
    else:
        problems.append((sid, "vep", "no annotated VCF"))

    # ---------------- pVACseq ----------------
    allep = glob.glob(f"{d}/pvacseq_out/MHC_Class_I/*.all_epitopes.tsv")
    aggr = glob.glob(f"{d}/pvacseq_out/MHC_Class_I/*.aggregated.tsv")
    filt = glob.glob(f"{d}/pvacseq_out/MHC_Class_I/*.filtered.tsv")

    if not allep:
        problems.append((sid, "missing", "no all_epitopes table"))
        rows.append(r)
        print(f"  {sid}  no epitopes")
        continue

    t = pd.read_csv(allep[0], sep="\t", low_memory=False)
    r["epitopes"] = len(t)

    if RANK_COL in t.columns:
        v = pd.to_numeric(t[RANK_COL], errors="coerce")
        r["strong"] = int((v < STRONG).sum())
        r["weak"] = int(((v >= STRONG) & (v < WEAK)).sum())
        r["scored"] = int(v.notna().sum())
        rank_pool.append(v.dropna().values)
        if v.notna().sum() and v.nunique() < 50:
            problems.append((sid, "flat",
                             f"only {v.nunique()} distinct rank values"))
    else:
        problems.append((sid, "column", f"{RANK_COL} absent"))
        v = pd.Series(dtype=float)

    # mutations reached, from the aggregated report
    if aggr:
        a = pd.read_csv(aggr[0], sep="\t", low_memory=False)
        r["mutations"] = len(a)
        gcol = next((c for c in a.columns if c.lower() == "gene"), None)
        if gcol:
            r["genes"] = a[gcol].nunique()

    if filt:
        r["pvac_filtered"] = max(0, sum(1 for _ in open(filt[0])) - 1)

    # peptide properties
    pcol = next((c for c in ["MT Epitope Seq", "Epitope Seq", "MT_Epitope_Seq"]
                 if c in t.columns), None)
    if pcol and len(v):
        peps = t.loc[v < STRONG, pcol].dropna().astype(str)
        peptide_sets[sid] = set(peps)
        r["unique_binder_peptides"] = peps.nunique()
        lens = t[pcol].dropna().astype(str).str.len()
        length_counts.update(lens.tolist())
        out_of_range = int(((lens < 8) | (lens > 11)).sum())
        if out_of_range:
            problems.append((sid, "length",
                             f"{out_of_range} peptides outside 8-11"))

    acol = next((c for c in ["HLA Allele", "Allele", "MHC"] if c in t.columns), None)
    if acol and len(v):
        for al in t.loc[v < STRONG, acol].dropna():
            allele_hits[str(al)] += 1
        r["alleles_used"] = t[acol].nunique()

    gcol2 = next((c for c in ["Gene Name", "Gene"] if c in t.columns), None)
    if gcol2 and len(v):
        gene_sets[sid] = set(t.loc[v < STRONG, gcol2].dropna())

    # ---------------- mhcflurry, for comparison ----------------
    mf = f"{COHORT}/{sid}/neoantigens/neoantigens_per_mutation.tsv"
    if os.path.exists(mf):
        n = pd.read_csv(mf, sep="\t")
        r["mf_mutations"] = len(n)
        if "strong" in n.columns:
            r["mf_strong"] = int(n.strong.sum())
        mfb = f"{COHORT}/{sid}/neoantigens/neoantigens_binders.tsv"
        if os.path.exists(mfb) and sid in gene_sets:
            b = pd.read_csv(mfb, sep="\t", low_memory=False)
            if "gene" in b.columns:
                mfg = set(b.gene.dropna())
                pvg = gene_sets[sid]
                r["genes_shared"] = len(mfg & pvg)
                r["genes_mf_only"] = len(mfg - pvg)
                r["genes_pv_only"] = len(pvg - mfg)
                union = mfg | pvg
                r["jaccard"] = round(len(mfg & pvg) / len(union), 3) if union else None

    rows.append(r)
    print(f"  {sid}  {r.get('epitopes', 0):,} epitopes, "
          f"{r.get('strong', 0)} strong")

t = pd.DataFrame(rows)
pd.set_option("display.width", 240)
pd.set_option("display.max_columns", 40)

# =====================================================================
print()
print("=" * 100)
print(" PER SAMPLE")
print("=" * 100)
cols = [c for c in ["sample", "cohort", "vep_records", "mutations", "genes",
                    "epitopes", "strong", "weak", "pvac_filtered",
                    "mf_strong", "jaccard"] if c in t.columns]
print(t[cols].to_string(index=False))

# =====================================================================
print()
print("=" * 100)
print(" OVERALL")
print("=" * 100)
for lab, col in [("samples with output", "epitopes"),
                 ("VEP records annotated", "vep_annotated"),
                 ("epitopes scored", "epitopes"),
                 ("strong binders (<0.5%)", "strong"),
                 ("weak binders (0.5-2%)", "weak")]:
    if col not in t.columns:
        continue
    if lab.startswith("samples"):
        print(f"\n  {lab:<26}{int((t[col] > 0).sum())} of {len(t)}")
    else:
        print(f"  {lab:<26}{int(t[col].sum()):>10,}")

if "pvac_filtered" in t.columns:
    print(f"\n  pVACseq's own filtered.tsv retained "
          f"{int(t.pvac_filtered.sum())} epitopes in total")
    print(f"    empty in {int((t.pvac_filtered == 0).sum())} of {len(t)} samples")
    print(f"    Its filter chain shortlists vaccine candidates, which is a")
    print(f"    different question from what this pipeline is measuring.")

# =====================================================================
print()
print("=" * 100)
print(" PERCENTILE RANK DISTRIBUTION")
print("=" * 100)
if rank_pool:
    v = np.concatenate(rank_pool)
    print(f"\n  {len(v):,} scored epitopes\n")
    print(f"  {'range':<16}{'n':>12}{'share':>9}")
    for lo, hi, lab in [(0, 0.5, "<0.5% strong"), (0.5, 2, "0.5-2% weak"),
                        (2, 10, "2-10%"), (10, 50, "10-50%"), (50, 101, ">50%")]:
        n = int(((v >= lo) & (v < hi)).sum())
        print(f"  {lab:<16}{n:>12,}{100*n/len(v):>8.1f}%")
    print(f"\n  median {np.median(v):.1f}, distinct values {len(np.unique(v)):,}")

if length_counts:
    print(f"\n  peptide lengths:")
    for L in sorted(length_counts):
        print(f"    {L}mer  {length_counts[L]:>9,}")

if allele_hits:
    total = sum(allele_hits.values())
    print(f"\n  {total:,} strong binder calls across {len(allele_hits)} alleles")
    for a, n in allele_hits.most_common(8):
        print(f"    {a:<16}{n:>7}  ({100*n/total:.1f}%)")

# =====================================================================
print()
print("=" * 100)
print(" AGREEMENT WITH mhcflurry")
print("=" * 100)
x = t.dropna(subset=["strong", "mf_strong"]) if "mf_strong" in t.columns else pd.DataFrame()
if len(x) > 3:
    rp, pp = stats.pearsonr(x.strong, x.mf_strong)
    rs, ps = stats.spearmanr(x.strong, x.mf_strong)
    print(f"\n  samples with both routes  {len(x)}")
    print(f"\n  strong binders per sample")
    print(f"    Pearson  r = {rp:.3f}  (p = {pp:.3g})")
    print(f"    Spearman r = {rs:.3f}  (p = {ps:.3g})")
    print(f"\n    NetMHCpan  total {int(x.strong.sum()):>6}, "
          f"median {x.strong.median():.0f}")
    print(f"    mhcflurry  total {int(x.mf_strong.sum()):>6}, "
          f"median {x.mf_strong.median():.0f}")
    ratio = x.strong.sum() / max(1, x.mf_strong.sum())
    print(f"\n    NetMHCpan calls {ratio:.2f}x as many strong binders")
    print(f"    The two calibrate the 0.5 percentile differently; the")
    print(f"    ranking of samples is what the correlation tests.")

    for coh, g in x.groupby("cohort"):
        if len(g) > 2:
            rr, pv2 = stats.pearsonr(g.strong, g.mf_strong)
            print(f"\n    {coh.upper():<6} r = {rr:.3f} (p = {pv2:.3g}), n = {len(g)}")

if "jaccard" in t.columns and t.jaccard.notna().any():
    print(f"\n  gene-level overlap")
    print(f"    shared           {int(t.genes_shared.sum()):>6}")
    print(f"    mhcflurry only   {int(t.genes_mf_only.sum()):>6}")
    print(f"    NetMHCpan only   {int(t.genes_pv_only.sum()):>6}")
    print(f"    median Jaccard   {t.jaccard.median():.3f}")

# =====================================================================
if len(peptide_sets) > 2:
    print()
    print("=" * 100)
    print(" PEPTIDE OVERLAP BETWEEN SAMPLES")
    print("=" * 100)
    ids = [s for s in peptide_sets if peptide_sets[s]]
    shares = []
    for i in range(len(ids)):
        for j in range(i + 1, len(ids)):
            a, b = peptide_sets[ids[i]], peptide_sets[ids[j]]
            u = a | b
            if u:
                shares.append(len(a & b) / len(u))
    if shares:
        print(f"\n  pairwise Jaccard over {len(shares)} pairs")
        print(f"    median {np.median(shares):.4f}, max {max(shares):.4f}")
        print(f"\n  Different mutations and different genotypes should give")
        print(f"  different peptides; overlap near zero is what confirms")
        print(f"  each sample was scored against its own HLA.")

# =====================================================================
print()
print("=" * 100)
print(f" PROBLEMS: {len(problems)}")
print("=" * 100)
if problems:
    by_kind = {}
    for sid, kind, msg in problems:
        by_kind.setdefault(kind, []).append((sid, msg))
    for kind in sorted(by_kind):
        print(f"\n  {kind} ({len(by_kind[kind])}):")
        for sid, msg in by_kind[kind][:12]:
            print(f"    {sid:<7} {msg}")
        if len(by_kind[kind]) > 12:
            print(f"    ... and {len(by_kind[kind])-12} more")
else:
    print("\n  none")

out = os.path.expanduser("~/immune_escape_project/results/verify_step7.tsv")
os.makedirs(os.path.dirname(out), exist_ok=True)
t.to_csv(out, sep="\t", index=False)
print(f"\nwritten to {out}")
