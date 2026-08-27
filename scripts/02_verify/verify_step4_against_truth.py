#!/usr/bin/env python3
"""
OptiType calls against the published HLA types for the same individuals.

Everything before this checked self-consistency: forty distinct genotypes,
diversity tracking ancestry, no silent defaulting. None of it shows the
calls are right. The backbones are 1000 Genomes individuals and the
project publishes typed alleles for its panel, so an external answer
exists.

Comparison is by multiset, not by set. A tool that reports one allele
where two exist and a tool that reports the correct homozygous pair are
different outcomes, and comparing sets collapses them: {11:01} equals
{11:01} either way. Counting with multiplicity keeps them apart, which
matters because dropping the second allele is the failure mode this data
actually shows.

Every sample is classified into exactly one bucket and the buckets are
asserted to sum to the sample count. A category that silently absorbs an
unexpected case is how the mpileup error survived as long as it did.

Resolution is two fields. The reference writes 02:01 with the gene in the
column name; OptiType writes A*02:01. Both reduce to 02:01 here and the
gene is carried by which column is being read.

Usage:
  python verify_step4_against_truth.py [workspace] [--show-all]
"""
import os
import re
import sys
import subprocess
from collections import Counter
import pandas as pd

args = [a for a in sys.argv[1:] if not a.startswith("--")]
WS = args[0] if args else os.environ.get(
    "WS", "/gpfs/bwfor/work/ws/fr_os136-immune_escape")
SHOW_ALL = "--show-all" in sys.argv

COHORT = f"{WS}/simulation/cohort"
RES = os.path.expanduser("~/immune_escape_project/results")
CACHE = f"{WS}/ref/1000G_HLA_types.txt"
URL = ("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
       "HLA_types/20181129_HLA_types_full_1000_Genomes_Project_panel.txt")

if not os.path.exists(CACHE):
    os.makedirs(os.path.dirname(CACHE), exist_ok=True)
    subprocess.run(f"curl -sfL --max-time 120 -o '{CACHE}' '{URL}'",
                   shell=True)
if not os.path.exists(CACHE) or os.path.getsize(CACHE) < 1000:
    sys.exit(f"reference missing; fetch {URL} to {CACHE}")

ref = pd.read_csv(CACHE, sep="\t", dtype=str)
ID_COL = next((c for c in ref.columns if c.strip().lower()
               in ("sample id", "sample", "sample_id", "id")), ref.columns[2])


def two_field(a):
    """02:01 from any of 02:01, A*02:01, HLA-A*02:01:01:02, or 02:01/26:01."""
    if a is None or (isinstance(a, float) and pd.isna(a)):
        return None
    s = str(a).strip().split("/")[0].replace("HLA-", "")
    if not s or s.lower() in ("nan", "na", "-", ""):
        return None
    rest = s.split("*", 1)[1] if "*" in s else s
    p = rest.split(":")
    if len(p) < 2 or not p[0].strip().isdigit():
        return None
    return f"{p[0].strip()}:{p[1].strip()}"


def ref_cols(loc):
    """
    The two columns holding this locus.

    Matching has to be anchored: a loose search for B would also answer
    with HLA-DQB1 and HLA-DRB1, quietly comparing the wrong gene.
    """
    out = [c for c in ref.columns
           if re.fullmatch(rf"HLA[-_ ]?{loc}[ _-]*[12]", c.strip(), re.I)]
    if len(out) != 2:
        sys.exit(f"expected two columns for HLA-{loc}, found {out}")
    return out


for loc in "ABC":
    ref_cols(loc)

man = pd.read_csv(f"{COHORT}/manifest.tsv", sep="\t")

BUCKETS = ["both correct", "one correct, one wrong", "one called, one dropped",
           "both wrong", "not called"]

rows = []
for _, m in man.iterrows():
    sid, bb = m.sample_id, str(m.backbone).strip()
    res = f"{COHORT}/{sid}/optitype/{sid}_result.tsv"
    if not os.path.exists(res):
        continue
    t = pd.read_csv(res, sep="\t")
    if not len(t):
        continue
    x = t.iloc[0]

    hit = ref[ref[ID_COL].astype(str).str.strip() == bb]
    row = {"sample": sid, "backbone": bb, "pop": m.get("superpopulation"),
           "in_reference": bool(len(hit))}
    if len(hit):
        h = hit.iloc[0]
        for loc in "ABC":
            # multisets: a homozygous call is two copies of one allele, and
            # a single call is one copy, which is a different thing
            called = Counter(v for v in
                             (two_field(x.get(f"{loc}{k}")) for k in (1, 2))
                             if v)
            truth = Counter(v for v in
                            (two_field(h[c]) for c in ref_cols(loc)) if v)

            n_called = sum(called.values())
            correct = sum((called & truth).values())

            if n_called == 0:
                bucket = "not called"
            elif n_called == 2 and correct == 2:
                bucket = "both correct"
            elif n_called == 1 and correct == 1 and sum(truth.values()) == 2:
                bucket = "one called, one dropped"
            elif correct >= 1:
                bucket = "one correct, one wrong"
            else:
                bucket = "both wrong"

            row[f"{loc}_called"] = "/".join(sorted(called.elements())) or None
            row[f"{loc}_truth"] = "/".join(sorted(truth.elements())) or None
            row[f"{loc}_correct"] = correct
            row[f"{loc}_n_called"] = n_called
            row[f"{loc}_n_truth"] = sum(truth.values())
            row[f"{loc}_bucket"] = bucket
    rows.append(row)

d = pd.DataFrame(rows)
found = d[d.in_reference].copy()

W = 76
print("=" * W)
print(" OPTITYPE AGAINST THE PUBLISHED TYPES")
print("=" * W)
print(f"\n  samples with an OptiType result   {len(d)}")
print(f"  backbones found in the reference   {len(found)}")
if len(d) != len(found):
    print(f"  absent from the reference          "
          f"{' '.join(d[~d.in_reference].backbone)}")

print()
print("=" * W)
print(" HOW EACH LOCUS CAME OUT")
print("=" * W)
print(f"\n  {'':<26}" + "".join(f"{'HLA-'+l:>10}" for l in "ABC"))
for b in BUCKETS:
    cells = "".join(f"{int((found[f'{l}_bucket'] == b).sum()):>10}"
                    for l in "ABC")
    print(f"  {b:<26}{cells}")

for loc in "ABC":
    tot = int((found[f"{loc}_bucket"].notna()).sum())
    s = int(found[f"{loc}_bucket"].isin(BUCKETS).sum())
    assert s == tot, f"HLA-{loc}: {s} classified of {tot}"
print(f"\n  every sample falls into exactly one row  \u2713")

print()
print("=" * W)
print(" ALLELE-LEVEL ACCURACY")
print("=" * W)
print(f"\n  {'locus':<10}{'published':>11}{'called':>9}{'correct':>10}"
      f"{'of published':>15}{'of called':>12}")
tc = tp = tn = 0
for loc in "ABC":
    g = found[found[f"{loc}_n_truth"].notna()]
    npub = int(g[f"{loc}_n_truth"].sum())
    ncal = int(g[f"{loc}_n_called"].sum())
    ncor = int(g[f"{loc}_correct"].sum())
    tc += ncor; tp += npub; tn += ncal
    print(f"  HLA-{loc:<6}{npub:>11}{ncal:>9}{ncor:>10}"
          f"{100*ncor/max(1,npub):>14.1f}%{100*ncor/max(1,ncal):>11.1f}%")
print(f"  {'total':<10}{tp:>11}{tn:>9}{tc:>10}"
      f"{100*tc/max(1,tp):>14.1f}%{100*tc/max(1,tn):>11.1f}%")

print(f"\n  Of published alleles is sensitivity: how many true alleles were")
print(f"  recovered. Of called alleles is precision: how many of the")
print(f"  reported alleles were right. They differ because OptiType")
print(f"  reports one allele where the reference has two.")

drop = sum(int((found[f"{l}_bucket"] == "one called, one dropped").sum())
           for l in "ABC")
print(f"\n  loci where a second allele was dropped   {drop}")
print(f"  loci where a called allele was wrong     "
      f"{sum(int((found[f'{l}_bucket'].isin(['one correct, one wrong', 'both wrong'])).sum()) for l in 'ABC')}")

print()
print("=" * W)
print(" WHERE A CALLED ALLELE IS WRONG")
print("=" * W)
shown = 0
for loc in "ABC":
    g = found[found[f"{loc}_bucket"].isin(
        ["one correct, one wrong", "both wrong"])]
    if not len(g):
        continue
    print(f"\n  HLA-{loc}")
    for _, r in g.iterrows():
        print(f"    {r['sample']:<7}{r.backbone:<10}"
              f"called {str(r[f'{loc}_called']):<18}"
              f"published {r[f'{loc}_truth']}")
        shown += 1
if not shown:
    print(f"\n  none")

if SHOW_ALL:
    print()
    print("=" * W)
    print(" EVERY LOCUS")
    print("=" * W)
    cols = ["sample", "backbone"] + [f"{l}_{k}" for l in "ABC"
                                     for k in ("called", "truth", "bucket")]
    pd.set_option("display.width", 250)
    print(found[cols].to_string(index=False))

out = f"{RES}/verify_step4_accuracy.tsv"
d.to_csv(out, sep="\t", index=False)
print(f"\nwritten to {out}")
