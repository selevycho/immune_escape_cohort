#!/usr/bin/env python3
"""
What the panel contains, and what each patient's mutations do inside it.

The panel is one fixed BED applied to all forty samples, as a clinical
panel would be. What differs between patients is only which of their
mutations fall inside it. This script demonstrates both halves of that
claim: the panel is identical everywhere, and the per-patient mutation
counts follow from their own mutation load rather than from any per-patient
selection.

For each patient it also reports how many mutations were lost outside the
panel, since that fraction is the price of using a panel at all and belongs
in the methods rather than in a footnote.

Usage:
  python check_panel.py <workspace>
"""
import sys, os
import numpy as np
import pandas as pd

WS = sys.argv[1] if len(sys.argv) > 1 else os.environ.get(
    "WS", "/gpfs/bwfor/work/ws/fr_os136-immune_escape")

PANEL = f"{WS}/simulation/panel/panel.bed"
MANIFEST = f"{WS}/simulation/cohort/manifest.tsv"
COHORT = f"{WS}/simulation/cohort"

NONSYN = {"Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
          "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
          "Splice_Site", "Nonstop_Mutation", "Translation_Start_Site"}
MIN_VAF = 0.05

# =====================================================================
print("=" * 78)
print(" THE PANEL")
print("=" * 78)

panel = pd.read_csv(PANEL, sep="\t", header=None,
                    names=["chrom", "start", "end", "name"])
panel["length"] = panel.end - panel.start

genes = set()
for n in panel.name:
    genes.update(str(n).split("+"))
genes.discard("")

print(f"\n  file      {PANEL}")
print(f"  intervals {len(panel):,}")
print(f"  genes     {len(genes):,}")
print(f"  total     {panel.length.sum()/1e6:.2f} Mb")
print(f"  interval  median {panel.length.median():.0f} bp, "
      f"range {panel.length.min()}-{panel.length.max()}")

print(f"\n  by chromosome:")
bych = panel.groupby("chrom").agg(
    n=("length", "size"), mb=("length", lambda x: x.sum()/1e6))
bych = bych.reindex(sorted(bych.index,
                    key=lambda c: (len(c), c)))
for c, r in bych.iterrows():
    print(f"    {c:<7} {r.n:>5} intervals   {r.mb:>6.2f} Mb")

# antigen presentation genes, checked by name
MACHINERY = ["B2M", "HLA-A", "HLA-B", "HLA-C", "TAP1", "TAP2",
             "TAPBP", "NLRC5", "PSMB8", "PSMB9", "CIITA"]
# The MHC is one contiguous block, not a set of per-gene intervals:
# OptiType and LOHHLA both realign whole regions against allele sequences
# and need flanking sequence, which exon-level intervals would not provide.
# Reporting its 3.5 Mb against each gene it contains would count the same
# bases a dozen times, so the block is reported once.
print(f"\n  antigen presentation genes in the panel:")
for g in MACHINERY:
    rows = panel[panel.name.astype(str).apply(lambda n: g in n.split("+"))]
    if not len(rows):
        print(f"    {g:<8} NO")
        continue
    in_block = rows.name.astype(str).str.contains("HLA_REGION_BLOCK").any()
    if in_block:
        print(f"    {g:<8} yes   inside the 3.5 Mb MHC block")
    else:
        print(f"    {g:<8} yes   {rows.length.sum():,} bp "
              f"in {len(rows)} interval(s)")

blk = panel[panel.name.astype(str).str.contains("HLA_REGION_BLOCK")]
if len(blk):
    r = blk.iloc[0]
    print(f"\n  MHC block: {r.chrom}:{r.start:,}-{r.end:,}  "
          f"({r.length/1e6:.2f} Mb)")
    print(f"    covers {len(str(r['name']).split('+')) - 1} named genes")
    print(f"    the remaining {len(panel)-1:,} intervals have a median "
          f"length of {panel[panel.length < 100000].length.median():.0f} bp")

with open(f"{WS}/simulation/panel/panel_genes.txt", "w") as fh:
    fh.write("\n".join(sorted(genes)) + "\n")
print(f"\n  full gene list written to "
      f"{WS}/simulation/panel/panel_genes.txt")

# =====================================================================
print()
print("=" * 78)
print(" IS THE PANEL THE SAME FOR EVERY PATIENT?")
print("=" * 78)

man = pd.read_csv(MANIFEST, sep="\t")
per_sample_panels = []
for sid in man.sample_id:
    p = f"{COHORT}/{sid}/panel.bed"
    if os.path.exists(p):
        per_sample_panels.append(sid)

if per_sample_panels:
    print(f"\n  WARNING: per-sample panel files found for "
          f"{len(per_sample_panels)} samples")
    print(f"    {' '.join(per_sample_panels[:10])}")
else:
    print(f"\n  No per-sample panel files exist. All forty samples were")
    print(f"  sliced against the single BED above.")

# confirm from the slicing logs, if they survive
import glob
logs = glob.glob(f"{WS}/logs/step1_*.log")
if logs:
    used = set()
    for lg in logs:
        for line in open(lg, errors="ignore"):
            if "panel.bed" in line:
                used.add(line.strip())
    print(f"  slicing logs reference {len(used) or 1} distinct panel path(s)")

# =====================================================================
print()
print("=" * 78)
print(" WHAT EACH PATIENT CONTRIBUTES")
print("=" * 78)

by_chrom = {c: g[["start", "end"]].to_numpy()
            for c, g in panel.groupby("chrom")}


def in_panel(chroms, positions):
    out = np.zeros(len(chroms), dtype=bool)
    for i, (c, p) in enumerate(zip(chroms, positions)):
        iv = by_chrom.get(c)
        if iv is None:
            continue
        pp = int(p) - 1
        out[i] = bool(((iv[:, 0] <= pp) & (pp < iv[:, 1])).any())
    return out


rows = []
mafs = {}
for cohort in ["brca", "ov"]:
    mafs[cohort] = pd.read_csv(f"{WS}/liftover/out/{cohort}.hg38.maf.tsv",
                               sep="\t", low_memory=False)

for _, m in man.iterrows():
    sid, cohort, barcode = m.sample_id, m.cohort, m.tcga_barcode
    d = mafs[cohort]
    x = d[d.Tumor_Sample_Barcode == barcode]
    n_all = len(x)

    x = x[x.Variant_Classification.isin(NONSYN)]
    n_nonsyn = len(x)

    den = x.t_ref_count.fillna(0) + x.t_alt_count.fillna(0)
    x = x[den > 0].copy()
    x["VAF"] = x.t_alt_count / den
    x = x[x.VAF >= MIN_VAF]
    n_vaf = len(x)

    n_panel = 0
    genes_hit = set()
    if n_vaf:
        mask = in_panel(x.Chromosome_hg38.values, x.Start_Position_hg38.values)
        inside = x[mask]
        n_panel = len(inside)
        genes_hit = set(inside.Hugo_Symbol.dropna())

    # what was actually injected
    truth = f"{COHORT}/{sid}/truth_set.tsv"
    n_injected = 0
    if os.path.exists(truth):
        n_injected = len(pd.read_csv(truth, sep="\t"))

    rows.append({
        "sample": sid, "cohort": cohort, "donor": barcode,
        "donor_mutations": n_all,
        "nonsyn": n_nonsyn,
        "vaf_filtered": n_vaf,
        "in_panel": n_panel,
        "pct_captured": round(100 * n_panel / n_vaf, 1) if n_vaf else 0,
        "genes_hit": len(genes_hit),
        "injected": n_injected,
    })

t = pd.DataFrame(rows)
pd.set_option("display.width", 200)
print()
print(t.drop(columns=["donor"]).to_string(index=False))

print()
print(f"  panel captures a median of {t.pct_captured.median():.1f}% "
      f"of each donor's usable mutations")
print(f"  range {t.pct_captured.min():.1f}% to {t.pct_captured.max():.1f}%")
print(f"  mutations in panel: {t.in_panel.sum():,} of "
      f"{t.vaf_filtered.sum():,} donor mutations")
print(f"  genes hit per patient: median {t.genes_hit.median():.0f}, "
      f"range {t.genes_hit.min()}-{t.genes_hit.max()}")

# =====================================================================
print()
print("=" * 78)
print(" WHICH PANEL GENES ARE ACTUALLY MUTATED, AND IN HOW MANY PATIENTS")
print("=" * 78)

hit_count = {}
for _, m in man.iterrows():
    d = mafs[m.cohort]
    x = d[d.Tumor_Sample_Barcode == m.tcga_barcode]
    x = x[x.Variant_Classification.isin(NONSYN)]
    den = x.t_ref_count.fillna(0) + x.t_alt_count.fillna(0)
    x = x[den > 0].copy()
    x["VAF"] = x.t_alt_count / den
    x = x[x.VAF >= MIN_VAF]
    if not len(x):
        continue
    mask = in_panel(x.Chromosome_hg38.values, x.Start_Position_hg38.values)
    for g in set(x[mask].Hugo_Symbol.dropna()):
        hit_count[g] = hit_count.get(g, 0) + 1

hc = pd.Series(hit_count).sort_values(ascending=False)
print(f"\n  panel genes carrying a mutation in at least one patient: "
      f"{len(hc)} of {len(genes)}")
print(f"  panel genes never mutated in this cohort: {len(genes) - len(hc)}")
print(f"\n  most frequently hit:")
print(f"    {'gene':<12}{'patients':>10}")
for g, n in hc.head(20).items():
    print(f"    {g:<12}{n:>10}")

out = os.path.expanduser("~/immune_escape_project/results/panel_audit.tsv")
os.makedirs(os.path.dirname(out), exist_ok=True)
t.to_csv(out, sep="\t", index=False)
hc.to_csv(out.replace(".tsv", "_gene_hits.tsv"), sep="\t", header=["patients"])
print(f"\nwritten to {out}")
