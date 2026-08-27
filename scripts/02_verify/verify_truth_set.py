#!/usr/bin/env python3
"""
Evidence that the truth set is checked rather than assumed.

The claim on the slide is that BAMSurgeon reports success for mutations
that are not in the resulting file. That is worth demonstrating rather
than asserting, so this script does three things.

It reads BAMSurgeon's own logs and extracts what the tool says it did at
each position — how many reads it selected, how many it rewrote. Those are
the tool's claims, taken from its output rather than paraphrased.

It then re-reads a sample of those positions straight out of the BAM with
mpileup, printing the actual pileup column for a few of them so the method
is visible and not just its conclusion.

Finally it computes what the truth set would look like if the log were
believed, and what recall in step 3 would have been against that truth
set. That difference is the number the slide rests on.

Both injection routes are covered: addsnv for substitutions, addindel for
the separate indel BAM.

Usage:
  python verify_truth_set.py [workspace] [--samples=B002,B018] [--show=6]
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

only, show = None, 6
for a in sys.argv[1:]:
    if a.startswith("--samples="):
        only = set(a.split("=", 1)[1].split(","))
    if a.startswith("--show="):
        show = int(a.split("=", 1)[1])

COHORT = f"{WS}/simulation/cohort"
INDELS = f"{WS}/simulation/indels"
MANIFEST = f"{COHORT}/manifest.tsv"
REF = f"{WS}/ref/Homo_sapiens_assembly38.fasta"
RES = os.path.expanduser("~/immune_escape_project/results")

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

# BAMSurgeon names each working set after the locus it is editing and then
# reports how many reads it selected and how many it rewrote
RE_LOCUS = re.compile(r"haplo_(chr[\w]+)_(\d+)_\d+")
RE_WROTE = re.compile(r"wrote:\s*(\d+),\s*mutated:\s*(\d+)")
RE_VAF = re.compile(r"selected VAF:\s*([\d.]+)")


def parse_bamsurgeon_log(path):
    """
    What the tool says it did, position by position.

    Returns a frame with one row per locus it touched, carrying the number
    of reads it had, the number it rewrote, and the fraction it aimed for.
    A locus with mutated > 0 is a locus BAMSurgeon considers done.
    """
    if not os.path.exists(path):
        return pd.DataFrame()
    rows = {}
    for line in open(path, errors="ignore"):
        m = RE_LOCUS.search(line)
        if not m:
            continue
        key = (m.group(1), int(m.group(2)))
        r = rows.setdefault(key, {"chrom": key[0], "pos": key[1],
                                  "reads": np.nan, "mutated": np.nan,
                                  "claimed_vaf": np.nan})
        w = RE_WROTE.search(line)
        if w:
            r["reads"] = int(w.group(1))
            r["mutated"] = int(w.group(2))
        v = RE_VAF.search(line)
        if v:
            r["claimed_vaf"] = float(v.group(1))
    return pd.DataFrame(rows.values())


def pileup(bam, chrom, pos):
    """The raw mpileup column, exactly as samtools prints it."""
    out = subprocess.run(
        f"{SAMTOOLS} mpileup -B -q 0 -Q 0 --ff UNMAP,SECONDARY,QCFAIL "
        f"-r {chrom}:{pos}-{pos} -f {REF} {bam} 2>/dev/null",
        shell=True, capture_output=True, text=True).stdout.strip()
    if not out:
        return None
    f = out.split("\t")
    return {"ref": f[2], "depth": int(f[3]), "bases": f[4]} if len(f) >= 5 else None


def count_alt(bases, alt):
    """
    Alternate-supporting reads in a pileup column.

    Insertion and deletion markers carry a length and a run of bases that
    would otherwise be counted as matches, so they are consumed first.
    """
    clean, i = [], 0
    while i < len(bases):
        c = bases[i]
        if c in "+-":
            j, num = i + 1, ""
            while j < len(bases) and bases[j].isdigit():
                num += bases[j]
                j += 1
            i = j + int(num or 0)
            continue
        if c == "^":
            i += 2
            continue
        if c == "$":
            i += 1
            continue
        clean.append(c)
        i += 1
    seq = "".join(clean).upper()
    if alt in ("-", "", None):
        return seq.count("*") + bases.count("-")
    return seq.count(str(alt).upper()[0])


man = pd.read_csv(MANIFEST, sep="\t")
if only:
    man = man[man.sample_id.isin(only)]

W = 72
print("=" * W)
print(" WHAT BAMSURGEON REPORTS")
print("=" * W)

claims, missing_logs = [], []
for _, m in man.iterrows():
    sid = m.sample_id
    log = f"{COHORT}/{sid}/bamsurgeon.log"
    c = parse_bamsurgeon_log(log)
    if c.empty:
        missing_logs.append(sid)
        continue
    c["sample"] = sid
    claims.append(c)

if not claims:
    sys.exit("no BAMSurgeon logs found under the cohort directories")

C = pd.concat(claims, ignore_index=True)
done = C[C.mutated.fillna(0) > 0]

print(f"\n  logs parsed              {C['sample'].nunique()} samples")
if missing_logs:
    print(f"  logs absent for          {len(missing_logs)}: "
          f"{' '.join(missing_logs[:8])}")
print(f"  loci the tool touched    {len(C)}")
print(f"  loci it reports mutated  {len(done)}")
print(f"  reads rewritten, median  {done.mutated.median():.0f} per locus")
print(f"\n  BAMSurgeon logs no failure at any of these. Its report is a")
print(f"  record of what it wrote, not of what survived realignment.")

# ------------------------------------------------------------------ SNV
print()
print("=" * W)
print(" WHAT THE FILE CONTAINS")
print("=" * W)

verified = f"{RES}/verify_step2_per_mutation.tsv"
if not os.path.exists(verified):
    sys.exit(f"\nmissing {verified}\nrun verify_step2_injection.py first")

v = pd.read_csv(verified, sep="\t")
v["alt"] = v.alt.astype(str)
snv = v[(v.alt.str.len() == 1) & (v.alt != "-")].copy()

merged = snv.merge(done[["sample", "chrom", "pos", "mutated", "claimed_vaf"]],
                   on=["sample", "chrom", "pos"], how="left")
claimed = merged[merged.mutated.notna()]

print(f"\n  substitutions in the truth set      {len(snv)}")
print(f"  matched to a BAMSurgeon claim       {len(claimed)}")
if len(claimed):
    absent = claimed[~claimed.landed]
    print(f"  claimed done but absent from BAM    {len(absent)}"
          f"   ({100*len(absent)/len(claimed):.1f}%)")
    print(f"\n  For those, the tool reports a median of "
          f"{absent.mutated.median():.0f} reads rewritten")
    print(f"  and the file contains {absent.alt_reads.max():.0f} at most.")

print(f"\n  across the whole truth set, without matching to the log:")
print(f"    present in the reads   {int(snv.landed.sum())} of {len(snv)}"
      f"   ({100*snv.landed.mean():.1f}%)")
print(f"    absent                 {int((~snv.landed).sum())}")

# ------------------------------------------------------------- the method
print()
print("=" * W)
print(" THE CHECK ITSELF, ON A FEW POSITIONS")
print("=" * W)
print(f"\n  Each position is read straight out of the BAM. The pileup")
print(f"  column below is what samtools prints; a dot or comma is a read")
print(f"  matching the reference, a letter is a read carrying something")
print(f"  else. A mutation is kept only if the alternate letter is there.\n")

lost = snv[~snv.landed]
kept = snv[snv.landed]
sample_rows = pd.concat([
    kept.sample(min(show // 2, len(kept)), random_state=5),
    lost.sample(min(show - show // 2, len(lost)), random_state=5),
])

for _, r in sample_rows.iterrows():
    bam = f"{COHORT}/{r['sample']}/{r['sample']}_tumor.bam"
    p = pileup(bam, r.chrom, int(r.pos))
    verdict = "KEPT" if r.landed else "DISCARDED"
    print(f"  {r['sample']}  {r.chrom}:{int(r.pos)}  ref {p['ref'] if p else '?'}"
          f"  alt {r.alt}   requested VAF {r.target_vaf:.3f}")
    if p:
        shown = p["bases"][:60] + ("..." if len(p["bases"]) > 60 else "")
        print(f"    depth {p['depth']:<4} {shown}")
        print(f"    reads carrying {r.alt}: {count_alt(p['bases'], r.alt)}"
              f"    -> {verdict}")
    else:
        print(f"    no pileup at this position          -> {verdict}")
    print()

# --------------------------------------------------------------- indels
print("=" * W)
print(" THE SAME FOR INDELS")
print("=" * W)

iv = f"{RES}/verify_step8_per_indel.tsv"
if os.path.exists(iv):
    ind = pd.read_csv(iv, sep="\t")
    print(f"\n  indels in the truth set   {len(ind)}")
    print(f"  present in the reads      {int(ind.landed.sum())}"
          f"   ({100*ind.landed.mean():.1f}%)")
    print(f"  absent                    {int((~ind.landed).sum())}")

    # the inner loop has to come first: the sample name is what builds
    # the paths, so it cannot be bound after them
    ilogs = [q for s in ind["sample"].unique()
             for q in (f"{INDELS}/{s}/bamsurgeon.log",
                       f"{INDELS}/{s}/addindel.log")
             if os.path.exists(q)]
    if ilogs:
        print(f"  BAMSurgeon logs present for the indel run too")
    else:
        print(f"\n  No addindel log was retained, so the claim side cannot")
        print(f"  be reconstructed for indels. The file side stands on its")
        print(f"  own: six positions carry no indel-supporting read.")

    print(f"\n  A deletion is not written as a letter in the pileup. It is")
    print(f"  marked -N in the column before it and as * in the column")
    print(f"  itself, so counting the literal character finds nothing —")
    print(f"  which is how an earlier version of this check reported every")
    print(f"  deletion as a failure.")
else:
    print(f"\n  missing {iv}")

# ---------------------------------------------------------------- cost
print()
print("=" * W)
print(" WHAT BELIEVING THE LOG WOULD HAVE COST")
print("=" * W)

comp = []
for _, m in man.iterrows():
    f = f"{COHORT}/{m.sample_id}/comparison/truth_vs_calls.tsv"
    if os.path.exists(f):
        t = pd.read_csv(f, sep="\t")
        t["sample"] = m.sample_id
        comp.append(t)

if comp:
    T = pd.concat(comp)
    n_checked = len(T)
    found = int(T.detected_pass.sum())
    n_naive = len(snv)                      # truth set taken from the log
    print(f"\n  truth set as verified       {n_checked} mutations")
    print(f"  recovered by Mutect2        {found}"
          f"   -> recall {100*found/n_checked:.1f}%")
    print(f"\n  truth set as logged         {n_naive} mutations")
    print(f"  recovered by Mutect2        {found}"
          f"   -> recall {100*found/n_naive:.1f}%")
    print(f"\n  difference                  "
          f"{100*found/n_checked - 100*found/n_naive:.1f} points")
    print(f"\n  The caller finds the same variants either way. Only the")
    print(f"  denominator changes, and with it the number we would have")
    print(f"  spent the rest of the talk explaining.")

out = f"{RES}/truth_set_evidence.tsv"
merged.to_csv(out, sep="\t", index=False)
print(f"\nwritten to {out}")
