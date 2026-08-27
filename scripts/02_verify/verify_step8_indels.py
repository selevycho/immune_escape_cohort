#!/usr/bin/env python3
"""
Step 8 verification: did the indels land in the reads?

Same correction as step 2: bases come from the alignment records rather
than from mpileup, which in samtools 1.23 silently requires the
proper-pair flag that realignment strips from edited reads.

Indels need more than a base lookup. An insertion adds query bases that
consume no reference, so it is found by looking for an I operation in the
CIGAR at the right offset. A deletion removes reference bases, so it is
found as a D operation spanning the position. Neither appears as a letter
at the position the way a substitution does, which is why counting
characters finds nothing.

Insertions and deletions are reported separately: their failure modes
differ, and pooling 16 insertions with 59 deletions hides whichever is
worse.

Usage:
  python verify_step8_indels.py [workspace] [--samples=B002,B018]
"""
import os
import re
import sys
import subprocess
import numpy as np
import pandas as pd

args = [a for a in sys.argv[1:] if not a.startswith("--")]
WS = args[0] if args else os.environ.get(
    "WS", "/gpfs/bwfor/work/ws/fr_os136-immune_escape")

only = None
for a in sys.argv[1:]:
    if a.startswith("--samples="):
        only = set(a.split("=", 1)[1].split(","))

INDELS = f"{WS}/simulation/indels"
TRUTH = f"{INDELS}/indel_truth_all.tsv"
MHC = ("chr6", 29600000, 33100000)

SAMTOOLS = None
for cand in ["samtools",
             os.path.expanduser("~/miniconda3/envs/bio_work/bin/samtools"),
             "/home/fr/fr_fr/fr_os136/miniconda3/envs/bio_work/bin/samtools"]:
    if subprocess.run(f"{cand} --version", shell=True,
                      capture_output=True).returncode == 0:
        SAMTOOLS = cand
        break
if SAMTOOLS is None:
    sys.exit("samtools not found; put its directory on PATH")

CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def indel_support(bam, chrom, pos, kind, window=6):
    """
    Reads carrying an indel at or near a position.

    A window is necessary rather than an exact match: bwa places an indel
    at the leftmost equivalent position, and in a repeat several positions
    are equivalent. Requiring the exact coordinate would score a correctly
    placed indel as absent whenever the repeat allowed a shift.
    """
    out = subprocess.run(
        f"{SAMTOOLS} view -F 0x404 {bam} {chrom}:{max(1,pos-window)}-{pos+window}",
        shell=True, capture_output=True, text=True).stdout
    depth = support = 0
    for line in out.splitlines():
        f = line.split("\t")
        if len(f) < 11:
            continue
        ref = int(f[3])
        covers = False
        hit = False
        for n, op in CIGAR.findall(f[5]):
            n = int(n)
            if op in "M=X":
                if ref <= pos < ref + n:
                    covers = True
                ref += n
            elif op == "I":
                if kind == "INS" and abs(ref - pos) <= window:
                    hit = True
            elif op in "DN":
                if kind == "DEL" and ref - window <= pos < ref + n + window:
                    hit = True
                if ref <= pos < ref + n:
                    covers = True
                ref += n
            elif op == "S":
                pass
        if covers or hit:
            depth += 1
        if hit:
            support += 1
    return depth, support


if not os.path.exists(TRUTH):
    sys.exit(f"missing {TRUTH}")

t = pd.read_csv(TRUTH, sep="\t")
if only:
    t = t[t["sample"].isin(only)]

print(f"counting indels from CIGAR operations, not mpileup")
print(f"{len(t)} indels across {t['sample'].nunique()} samples\n", flush=True)

rows = []
for sid, g in t.groupby("sample"):
    bam = None
    # step 8 writes into a copy that already carries the substitutions,
    # so the name records both rounds of editing
    for cand in [f"{INDELS}/{sid}/{sid}_tumor_snv_indel.bam",
                 f"{INDELS}/{sid}/{sid}_tumor_indel.bam",
                 f"{INDELS}/{sid}/{sid}_indel.bam"]:
        if os.path.exists(cand):
            bam = cand
            break
    if bam is None:
        print(f"  {sid}  no indel BAM found")
        continue

    print(f"  {sid}  {len(g)} indels ...", end="", flush=True)
    n_ok = 0
    for _, r in g.iterrows():
        chrom = r.get("chrom", r.get("Chromosome_hg38"))
        pos = int(r.get("pos", r.get("Start_Position_hg38")))
        kind = str(r.get("type", r.get("Variant_Type", ""))).upper()
        kind = "INS" if "INS" in kind else "DEL"
        want = float(r.get("vaf", r.get("VAF", np.nan)))

        d, s = indel_support(bam, chrom, pos, kind)
        landed = s > 0
        n_ok += landed
        rows.append({
            "sample": sid, "chrom": chrom, "pos": pos, "type": kind,
            "gene": r.get("gene", r.get("Hugo_Symbol")),
            "consequence": r.get("consequence", r.get("Variant_Classification")),
            "target_vaf": want, "depth": d, "support": s,
            "observed_vaf": s / d if d else 0.0,
            "in_mhc": (chrom == MHC[0] and MHC[1] <= pos < MHC[2]),
            "landed": landed,
        })
    print(f" {n_ok}/{len(g)}")

d = pd.DataFrame(rows)
if d.empty:
    sys.exit("nothing verified")

pd.set_option("display.width", 220)


def wilson(k, m):
    if m == 0:
        return (np.nan, np.nan)
    z, p = 1.96, k / m
    den = 1 + z * z / m
    c = (p + z * z / (2 * m)) / den
    h = z * np.sqrt(p * (1 - p) / m + z * z / (4 * m * m)) / den
    return (100 * max(0, c - h), 100 * min(1, c + h))


print()
print("=" * 84)
print(" OVERALL")
print("=" * 84)
n, landed = len(d), int(d.landed.sum())
lo, hi = wilson(landed, n)
print(f"\n  indels injected      {n}")
print(f"  present in the reads {landed}   ({100*landed/n:.1f}%, "
      f"95% CI {lo:.0f}-{hi:.0f}%)")
print(f"  samples with any     {d['sample'].nunique()}")

print(f"\n  {'type':<8}{'n':>6}{'landed':>9}{'rate':>9}{'95% CI':>16}")
for ty, g in d.groupby("type"):
    k = int(g.landed.sum())
    l2, h2 = wilson(k, len(g))
    print(f"  {ty:<8}{len(g):>6}{k:>9}{100*k/len(g):>8.1f}%"
          f"{f'{l2:.0f}-{h2:.0f}%':>16}")

if "consequence" in d.columns and d.consequence.notna().any():
    print(f"\n  by consequence:")
    for c, g in d.groupby("consequence"):
        print(f"    {str(c):<26}{int(g.landed.sum()):>4} of {len(g)}")

ok = d[d.landed & (d.target_vaf > 0) & (d.observed_vaf > 0)]
if len(ok) > 3:
    r = float(np.corrcoef(ok.target_vaf, ok.observed_vaf)[0, 1])
    ratio = float((ok.observed_vaf / ok.target_vaf).median())
    print(f"\n  allele fraction")
    print(f"    correlation            r = {r:.3f}   (n = {len(ok)})")
    print(f"    median observed/target {ratio:.3f}")
    print(f"    mean deviation         "
          f"{(ok.observed_vaf - ok.target_vaf).mean():+.4f}")

print()
print("=" * 84)
print(" WHERE THEY FAILED")
print("=" * 84)
lost = d[~d.landed]
if not len(lost):
    print("\n  none failed")
else:
    print(f"\n  {len(lost)} did not land\n")
    print(f"  {'sample':<8}{'position':<24}{'type':<6}{'gene':<12}"
          f"{'depth':>7}{'target VAF':>12}")
    for _, r in lost.iterrows():
        print(f"  {r['sample']:<8}{r.chrom + ':' + str(int(r.pos)):<24}"
              f"{r.type:<6}{str(r.gene):<12}{int(r.depth):>7}"
              f"{r.target_vaf:>12.3f}")
    genes = sorted(set(lost.gene.dropna()))
    if genes:
        print(f"\n  genes involved: {', '.join(genes)}")
    print(f"\n  depth where lost   median {lost.depth.median():.0f}x")
    print(f"  depth where landed median {d[d.landed].depth.median():.0f}x")

outp = os.path.expanduser("~/immune_escape_project/results/verify_step8.tsv")
os.makedirs(os.path.dirname(outp), exist_ok=True)
d.groupby("sample").agg(
    injected=("landed", "size"), landed=("landed", "sum"),
    ins=("type", lambda x: int((x == "INS").sum())),
    dels=("type", lambda x: int((x == "DEL").sum())),
    median_depth=("depth", "median")).reset_index().to_csv(
        outp, sep="\t", index=False)
d.to_csv(outp.replace(".tsv", "_per_indel.tsv"), sep="\t", index=False)
print(f"\nwritten to {outp}")
