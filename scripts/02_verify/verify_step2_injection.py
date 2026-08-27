#!/usr/bin/env python3
"""
Step 2 verification: did the substitutions land in the reads?

Bases are counted from the alignment records directly, not through
mpileup. samtools 1.23 mpileup silently requires the proper-pair flag and
gives no option to disable it, which is fatal here: BAMSurgeon edits reads
and realigns them through bwa, and realignment strips proper-pair status
from precisely the reads it rewrote. mpileup therefore discards
mutation-carrying reads preferentially. At one position checked by hand,
four of forty reads carried the alternate base, none of the four was
flagged proper-pair, and mpileup reported a depth of one.

Reading the alignments instead means walking each read's CIGAR to find
which query base sits at the reference position, since soft clips shift
the query without shifting the reference.

Duplicates are excluded, unmapped and secondary records are excluded, and
nothing else is filtered — a read that lost its proper-pair flag during
realignment still carries the base it was given.

The MHC is reported separately throughout. Realignment against a reference
carrying alternate haplotypes there behaves differently from the rest of
the panel, and pooling the two hides it.

Usage:
  python verify_step2_injection.py [workspace] [--samples=B002,B018]
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

COHORT = f"{WS}/simulation/cohort"
MANIFEST = f"{COHORT}/manifest.tsv"
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


def base_at(start, cigar, seq, target):
    """
    The query base aligned to a reference position, or None if the read
    does not reach it. A deletion covering the position returns '*'.
    """
    ref, qry = start, 0
    for n, op in CIGAR.findall(cigar):
        n = int(n)
        if op in "M=X":
            if ref <= target < ref + n:
                return seq[qry + (target - ref)]
            ref += n
            qry += n
        elif op in "IS":
            qry += n
        elif op in "DN":
            if ref <= target < ref + n:
                return "*"
            ref += n
    return None


def pileup_direct(bam, chrom, pos, alt):
    """Depth and alternate-supporting reads, read from the alignments."""
    out = subprocess.run(
        f"{SAMTOOLS} view -F 0x404 {bam} {chrom}:{pos}-{pos}",
        shell=True, capture_output=True, text=True).stdout
    depth = n_alt = 0
    for line in out.splitlines():
        f = line.split("\t")
        if len(f) < 11:
            continue
        b = base_at(int(f[3]), f[5], f[9], pos)
        if b is None:
            continue
        depth += 1
        if str(alt) in ("-", ""):
            if b == "*":
                n_alt += 1
        elif b.upper() == str(alt).upper()[0]:
            n_alt += 1
    return depth, n_alt


man = pd.read_csv(MANIFEST, sep="\t")
if only:
    man = man[man.sample_id.isin(only)]

print(f"counting bases from alignments, not mpileup")
print(f"checking {len(man)} samples\n", flush=True)

per_mut, per_sample = [], []

for i, m in man.iterrows():
    sid = m.sample_id
    truth = f"{COHORT}/{sid}/truth_set.tsv"
    tbam = f"{COHORT}/{sid}/{sid}_tumor.bam"
    nbam = f"{COHORT}/{sid}/{sid}_normal.bam"
    if not (os.path.exists(truth) and os.path.exists(tbam)):
        print(f"  {sid}  truth set or tumour BAM missing")
        continue

    t = pd.read_csv(truth, sep="\t")
    # step 2 places substitutions only; indels go through addindel in step 8
    n_all = len(t)
    if "Variant_Type" in t.columns:
        t = t[t.Variant_Type == "SNP"]
    n_ind = n_all - len(t)

    print(f"  [{i+1:>2}/{len(man)}] {sid}  {len(t)} SNVs"
          f"{f' ({n_ind} indels deferred)' if n_ind else ''} ...",
          end="", flush=True)

    rows = []
    for _, r in t.iterrows():
        chrom = r.get("Chromosome_hg38", r.get("chrom"))
        pos = int(r.get("Start_Position_hg38", r.get("pos")))
        alt = str(r.get("Tumor_Seq_Allele2", r.get("alt")))
        want = float(r.get("VAF", np.nan))

        d_t, a_t = pileup_direct(tbam, chrom, pos, alt)
        d_n, a_n = pileup_direct(nbam, chrom, pos, alt)

        rows.append({
            "sample": sid, "cohort": m.cohort,
            "chrom": chrom, "pos": pos, "alt": alt,
            "gene": r.get("Hugo_Symbol"),
            "target_vaf": want,
            "tumour_depth": d_t, "alt_reads": a_t,
            "observed_vaf": a_t / d_t if d_t else 0.0,
            "normal_depth": d_n, "normal_alt": a_n,
            "in_mhc": (chrom == MHC[0] and MHC[1] <= pos < MHC[2]),
            "landed": a_t > 0,
        })

    d = pd.DataFrame(rows)
    per_mut.append(d)

    ok = d[d.landed & (d.target_vaf > 0)]
    s = {"sample": sid, "cohort": m.cohort,
         "injected": len(d), "landed": int(d.landed.sum()),
         "landed_pct": round(100 * d.landed.mean(), 1),
         "in_mhc": int(d.in_mhc.sum()),
         "mhc_landed": int(d[d.in_mhc].landed.sum()) if d.in_mhc.any() else 0,
         "median_target_vaf": round(d.target_vaf.median(), 3),
         "median_observed_vaf": round(ok.observed_vaf.median(), 3) if len(ok) else None,
         "median_depth": round(d.tumour_depth.median(), 1),
         "normal_contam": int((d.normal_alt > 0).sum())}
    if len(ok) > 2:
        s["vaf_r"] = round(float(np.corrcoef(ok.target_vaf, ok.observed_vaf)[0, 1]), 3)
        s["vaf_ratio"] = round(float((ok.observed_vaf / ok.target_vaf).median()), 3)
    per_sample.append(s)
    print(f" {s['landed']}/{s['injected']} ({s['landed_pct']}%)")

mut = pd.concat(per_mut)
smp = pd.DataFrame(per_sample)
pd.set_option("display.width", 220)

print()
print("=" * 88)
print(" PER SAMPLE")
print("=" * 88)
cols = [c for c in ["sample", "cohort", "injected", "landed", "landed_pct",
                    "in_mhc", "mhc_landed", "median_target_vaf",
                    "median_observed_vaf", "vaf_r", "vaf_ratio",
                    "normal_contam"] if c in smp.columns]
print(smp[cols].to_string(index=False))

print()
print("=" * 88)
print(" OVERALL")
print("=" * 88)
n, landed = len(mut), int(mut.landed.sum())
print(f"\n  substitutions injected   {n}")
print(f"  present in the reads     {landed}   ({100*landed/n:.1f}%)")

out_m = mut[~mut.in_mhc]
in_m = mut[mut.in_mhc]
print(f"\n  outside the MHC          {int(out_m.landed.sum())} of {len(out_m)}"
      f"   ({100*out_m.landed.mean():.1f}%)")
if len(in_m):
    print(f"  inside the MHC           {int(in_m.landed.sum())} of {len(in_m)}"
          f"   ({100*in_m.landed.mean():.1f}%)")

ok = mut[mut.landed & (mut.target_vaf > 0)]
if len(ok) > 5:
    r = float(np.corrcoef(ok.target_vaf, ok.observed_vaf)[0, 1])
    ratio = float((ok.observed_vaf / ok.target_vaf).median())
    print(f"\n  allele fraction, landed mutations only")
    print(f"    correlation            r = {r:.3f}   (n = {len(ok)})")
    print(f"    median observed/target {ratio:.3f}")
    print(f"    median target          {ok.target_vaf.median():.3f}")
    print(f"    median observed        {ok.observed_vaf.median():.3f}")
    print(f"    mean deviation         {(ok.observed_vaf - ok.target_vaf).mean():+.4f}")

print(f"\n  positions where the normal also carries the alt base: "
      f"{int((mut.normal_alt > 0).sum())}")

print()
print("=" * 88)
print(" LANDING RATE BY REQUESTED FRACTION")
print("=" * 88)
b = [0, .05, .10, .15, .20, .30, 1.01]
l = ["<5%", "5-10%", "10-15%", "15-20%", "20-30%", ">30%"]
mut["vaf_bin"] = pd.cut(mut.target_vaf, bins=b, labels=l, right=False)
g = out_m.assign(vaf_bin=pd.cut(out_m.target_vaf, bins=b, labels=l,
                                right=False)) \
         .groupby("vaf_bin", observed=False).agg(
             n=("landed", "size"), landed=("landed", "sum"),
             alt=("alt_reads", "median"))
print(f"\n  outside the MHC only\n")
print(f"  {'requested':<11}{'n':>7}{'landed':>9}{'rate':>9}{'alt reads':>12}")
for idx, r in g.iterrows():
    if r.n == 0:
        continue
    print(f"  {str(idx):<11}{int(r.n):>7}{int(r.landed):>9}"
          f"{100*r.landed/r.n:>8.1f}%{r.alt:>12.0f}")

print()
print("=" * 88)
print(" WHERE MUTATIONS WERE LOST")
print("=" * 88)
lost = mut[~mut.landed]
print(f"\n  {len(lost)} of {n} did not land\n")
cats = {
    "inside the MHC": lost.in_mhc,
    "outside, depth below 5": (~lost.in_mhc) & (lost.tumour_depth < 5),
    "outside, covered": (~lost.in_mhc) & (lost.tumour_depth >= 5),
}
for name, mask in cats.items():
    k = int(mask.sum())
    print(f"  {name:<26}{k:>5}   {100*k/max(1,len(lost)):>5.0f}% of losses")

rest = lost[(~lost.in_mhc) & (lost.tumour_depth >= 5)]
if len(rest):
    print(f"\n  the covered group:")
    print(f"    median depth           {rest.tumour_depth.median():.0f}x")
    print(f"    median requested VAF   {rest.target_vaf.median():.3f}")
    print(f"    by chromosome:")
    for c, k in rest.chrom.value_counts().head(6).items():
        tot = int((mut.chrom == c).sum())
        print(f"      {c:<8}{k:>4} of {tot:>5}")

if len(lost):
    print(f"\n  depth where lost   median {lost.tumour_depth.median():.0f}x")
    print(f"  depth where landed median {mut[mut.landed].tumour_depth.median():.0f}x")

outp = os.path.expanduser("~/immune_escape_project/results/verify_step2.tsv")
os.makedirs(os.path.dirname(outp), exist_ok=True)
smp.to_csv(outp, sep="\t", index=False)
mut.to_csv(outp.replace(".tsv", "_per_mutation.tsv"), sep="\t", index=False)
print(f"\nwritten to {outp}")
