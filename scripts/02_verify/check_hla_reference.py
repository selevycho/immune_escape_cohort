#!/usr/bin/env python3
"""
Is the reference table the right table, and are we reading the right rows?

An accuracy of 78% is low enough to be a real limitation of typing at this
depth and also low enough to be a lookup error, and those look identical
from the accuracy figure alone. The mpileup episode came from trusting a
measurement without asking what it was measuring, so the same question is
asked here before the number reaches a slide.

Four things are checked.

Does the population recorded for each backbone in our manifest agree with
the population recorded in the reference? If we matched the wrong row, the
alleles would be a stranger's and the population would usually disagree.
This is the strongest single test available without new data.

Does any sample ID appear twice in the reference? A duplicate would mean
the first match wins arbitrarily.

Are the reference values well formed for every one of our forty, or are
some blank, ambiguous, or truncated in a way that makes a mismatch
inevitable?

And are the alleles we call plausible at all? An allele absent from the
whole reference panel of 2500 people is more likely a typing artefact than
a rare discovery.

Usage:
  python check_hla_reference.py [workspace]
"""
import os
import re
import sys
from collections import Counter
import pandas as pd

WS = sys.argv[1] if len(sys.argv) > 1 else os.environ.get(
    "WS", "/gpfs/bwfor/work/ws/fr_os136-immune_escape")
COHORT = f"{WS}/simulation/cohort"
CACHE = f"{WS}/ref/1000G_HLA_types.txt"

ref = pd.read_csv(CACHE, sep="\t", dtype=str)
man = pd.read_csv(f"{COHORT}/manifest.tsv", sep="\t")

W = 74
print("=" * W)
print(" 1. WHAT THIS FILE IS")
print("=" * W)
print(f"\n  path        {CACHE}")
print(f"  rows        {len(ref)}")
print(f"  columns     {len(ref.columns)}")
print(f"  populations {ref.Population.nunique()}")
print(f"  regions     {sorted(ref.Region.dropna().unique())}")
print(f"\n  first two rows:")
print(ref.head(2).to_string(index=False))

dups = ref[ID := ref.columns[2]].duplicated().sum() if False else \
    ref[ref.columns[2]].duplicated().sum()
ID_COL = ref.columns[2]
print(f"\n  duplicate sample IDs: {int(dups)}")
if dups:
    print(f"    {ref[ref[ID_COL].duplicated(keep=False)][ID_COL].unique()[:10]}")

print()
print("=" * W)
print(" 2. DO WE MATCH THE RIGHT PERSON")
print("=" * W)
print(f"\n  If a backbone were matched to the wrong row, the population")
print(f"  recorded there would usually differ from ours.\n")

rows = []
for _, m in man.iterrows():
    bb = str(m.backbone).strip()
    hit = ref[ref[ID_COL].astype(str).str.strip() == bb]
    rows.append({
        "sample": m.sample_id, "backbone": bb,
        "our_pop": str(m.get("population", "")).strip(),
        "our_super": str(m.get("superpopulation", "")).strip(),
        "ref_pop": str(hit.iloc[0].Population).strip() if len(hit) else None,
        "ref_region": str(hit.iloc[0].Region).strip() if len(hit) else None,
        "found": bool(len(hit)),
    })
c = pd.DataFrame(rows)
c["pop_agrees"] = c.our_pop.str.upper() == c.ref_pop.str.upper()
c["region_agrees"] = c.our_super.str.upper() == c.ref_region.str.upper()

print(f"  found in the reference     {int(c.found.sum())} of {len(c)}")
print(f"  population agrees          {int(c.pop_agrees.sum())} of {len(c)}")
print(f"  superpopulation agrees     {int(c.region_agrees.sum())} of {len(c)}")

bad = c[c.found & ~c.pop_agrees]
if len(bad):
    print(f"\n  disagreements:")
    print(bad[["sample", "backbone", "our_pop", "ref_pop",
               "our_super", "ref_region"]].to_string(index=False))
else:
    print(f"\n  Every backbone matches a row of the same population. The")
    print(f"  lookup is right; the alleles belong to the person we sliced.")

print()
print("=" * W)
print(" 3. ARE THE REFERENCE VALUES USABLE")
print("=" * W)


def cols_for(loc):
    return [x for x in ref.columns
            if re.fullmatch(rf"HLA[-_ ]?{loc}[ _-]*[12]", x.strip(), re.I)]


ours = ref[ref[ID_COL].astype(str).str.strip().isin(
    man.backbone.astype(str).str.strip())]
print(f"\n  {'locus':<8}{'blank':>8}{'ambiguous':>12}{'one field':>12}"
      f"{'usable':>9}")
for loc in "ABC":
    blank = amb = short = ok = 0
    for cc in cols_for(loc):
        for v in ours[cc]:
            s = str(v).strip()
            if not s or s.lower() in ("nan", "na", "-"):
                blank += 1
            elif "/" in s or ";" in s:
                amb += 1
            elif s.count(":") == 0:
                short += 1
            else:
                ok += 1
    print(f"  HLA-{loc:<4}{blank:>8}{amb:>12}{short:>12}{ok:>9}")

print()
print("=" * W)
print(" 4. ARE OUR CALLS PLAUSIBLE ALLELES")
print("=" * W)
print(f"\n  An allele we report that appears nowhere in {len(ref)} typed")
print(f"  individuals is more likely an artefact than a discovery.\n")


def tf(a):
    if a is None or (isinstance(a, float) and pd.isna(a)):
        return None
    s = str(a).strip().split("/")[0].replace("HLA-", "")
    if not s or s.lower() in ("nan", "na", "-"):
        return None
    rest = s.split("*", 1)[1] if "*" in s else s
    p = rest.split(":")
    return f"{p[0].strip()}:{p[1].strip()}" if len(p) >= 2 \
        and p[0].strip().isdigit() else None


for loc in "ABC":
    pool = Counter()
    for cc in cols_for(loc):
        for v in ref[cc]:
            a = tf(v)
            if a:
                pool[a] += 1
    total = sum(pool.values())

    called, unseen, rare = Counter(), [], []
    for _, m in man.iterrows():
        p = f"{COHORT}/{m.sample_id}/optitype/{m.sample_id}_result.tsv"
        if not os.path.exists(p):
            continue
        t = pd.read_csv(p, sep="\t")
        if not len(t):
            continue
        for k in (1, 2):
            a = tf(t.iloc[0].get(f"{loc}{k}"))
            if not a:
                continue
            called[a] += 1
            if pool[a] == 0:
                unseen.append((m.sample_id, a))
            elif pool[a] / total < 0.002:
                rare.append((m.sample_id, a, pool[a]))

    print(f"  HLA-{loc}")
    print(f"    distinct alleles we call        {len(called)}")
    print(f"    never seen in the reference     {len(unseen)}"
          f"{'   ' + ', '.join(f'{s}:{a}' for s, a in unseen[:6]) if unseen else ''}")
    print(f"    seen but below 0.2% frequency   {len(rare)}"
          f"{'   ' + ', '.join(f'{s}:{a}({n})' for s, a, n in rare[:5]) if rare else ''}")

out = os.path.expanduser(
    "~/immune_escape_project/results/hla_reference_check.tsv")
c.to_csv(out, sep="\t", index=False)
print(f"\nwritten to {out}")
