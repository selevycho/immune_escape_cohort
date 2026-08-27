#!/usr/bin/env python3
"""
Why one indel was not called when the reads plainly carry it.

DOCK10 in B002 sits at 29% allele fraction with ten supporting reads and
Mutect2 did not report it. Eleven of the other twelve misses are below 16%
and need no explanation; this one does.

Three possibilities, and they are distinguishable.

The variant was never emitted. Mutect2 builds haplotypes over an active
region and decides what to consider; if the region was not made active, or
the assembly did not produce a haplotype carrying the indel, nothing
appears in the VCF at all - not even a filtered record.

It was emitted and filtered. Then a record exists with a filter tag, and
the tag says which model rejected it.

It was emitted somewhere else. A long indel in a repeat can be described
at a position further away than the matching window allows, or as part of
a larger event, and the comparison would score it missed.

The raw VCF is checked as well as the filtered one, because
FilterMutectCalls only annotates - a variant absent from the raw file was
never considered.

Usage:
  python inspect_missed_indel.py [workspace] [--all]
  python inspect_missed_indel.py [workspace] --pos=B002:chr2:224874130
"""
import os
import re
import sys
import gzip
import subprocess
import pandas as pd

args = [a for a in sys.argv[1:] if not a.startswith("--")]
WS = args[0] if args else os.environ.get(
    "WS", "/gpfs/bwfor/work/ws/fr_os136-immune_escape")
ALL = "--all" in sys.argv
ONE = None
for a in sys.argv[1:]:
    if a.startswith("--pos="):
        ONE = a.split("=", 1)[1]

INDELS = f"{WS}/simulation/indels"
RES = os.path.expanduser("~/immune_escape_project/results")
REC = f"{RES}/indel_recall.tsv"
REF = f"{WS}/ref/Homo_sapiens_assembly38.fasta"

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

if not os.path.exists(REC):
    sys.exit(f"missing {REC}; run compare_indel_calls.py first")

d = pd.read_csv(REC, sep="\t")
miss = d[~d.found].copy()

if ONE:
    sid, chrom, pos = ONE.split(":")
    targets = miss[(miss["sample"] == sid) & (miss.chrom == chrom) &
                   (miss.pos == int(pos))]
    if not len(targets):
        targets = pd.DataFrame([{"sample": sid, "chrom": chrom,
                                 "pos": int(pos), "type": "?", "gene": "?",
                                 "target_vaf": float("nan"),
                                 "observed_vaf": float("nan"),
                                 "depth": 0, "support": 0}])
elif ALL:
    targets = miss
else:
    # the unexplained one: enough support that low coverage is not the answer
    targets = miss[miss.support >= 8]
    if not len(targets):
        targets = miss.nlargest(1, "support")

CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def vcf_near(path, chrom, pos, span=200):
    out = []
    if not os.path.exists(path):
        return out
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[0] != chrom or abs(int(f[1]) - pos) > span:
                continue
            out.append({"pos": int(f[1]), "ref": f[3], "alt": f[4],
                        "qual": f[5], "filter": f[6],
                        "info": f[7][:120]})
    return out


def reads_at(bam, chrom, pos, span=3):
    """Reads covering the position, with what each carries there."""
    out = subprocess.run(
        f"{SAMTOOLS} view -F 0x404 {bam} {chrom}:{pos-span}-{pos+span}",
        shell=True, capture_output=True, text=True).stdout
    carriers, plain, clipped = 0, 0, 0
    mapqs = []
    for line in out.splitlines():
        f = line.split("\t")
        if len(f) < 11:
            continue
        ref = int(f[3])
        has_indel = False
        soft = False
        for n, op in CIGAR.findall(f[5]):
            n = int(n)
            if op in "M=X":
                ref += n
            elif op == "I":
                if abs(ref - pos) <= span:
                    has_indel = True
            elif op in "DN":
                if ref - span <= pos < ref + n + span:
                    has_indel = True
                ref += n
            elif op == "S":
                soft = True
        mapqs.append(int(f[4]))
        if has_indel:
            carriers += 1
        else:
            plain += 1
        if soft:
            clipped += 1
    return carriers, plain, clipped, mapqs


W = 74
for _, r in targets.iterrows():
    sid, chrom, pos = r["sample"], r.chrom, int(r.pos)
    print("=" * W)
    print(f" {sid}  {chrom}:{pos}  {r.type}  {r.gene}")
    print("=" * W)
    print(f"\n  requested VAF {r.target_vaf:.3f}, observed {r.observed_vaf:.3f}"
          f", depth {int(r.depth)}, supporting reads {int(r.support)}")

    tbam = f"{INDELS}/{sid}/{sid}_tumor_snv_indel.bam"
    nbam = f"{WS}/simulation/cohort/{sid}/{sid}_normal.bam"

    print(f"\n  --- what the tumour reads show ---")
    c, p, cl, mq = reads_at(tbam, chrom, pos)
    print(f"    reads covering            {c + p}")
    print(f"    carrying an indel here    {c}")
    print(f"    plain                     {p}")
    print(f"    soft-clipped              {cl}")
    if mq:
        print(f"    median MAPQ               {sorted(mq)[len(mq)//2]}")

    print(f"\n  --- what the normal shows ---")
    if os.path.exists(nbam):
        cn, pn, cln, mqn = reads_at(nbam, chrom, pos)
        print(f"    reads covering            {cn + pn}")
        print(f"    carrying an indel here    {cn}"
              f"{'   <- germline, correctly rejected' if cn > 1 else ''}")
    else:
        print(f"    normal BAM not found")

    for label, path in [("raw", f"{INDELS}/{sid}/calls/{sid}.raw.vcf.gz"),
                        ("filtered",
                         f"{INDELS}/{sid}/calls/{sid}.filtered.vcf.gz")]:
        near = vcf_near(path, chrom, pos)
        print(f"\n  --- {label} VCF, within 200 bp ---")
        if not os.path.exists(path):
            print(f"    file absent")
            continue
        if not near:
            print(f"    nothing emitted anywhere near this position")
            continue
        for v in near:
            d_ = v["pos"] - pos
            kind = "indel" if len(v["ref"]) != len(v["alt"]) else "SNV"
            print(f"    {v['pos']} ({d_:+d})  {v['ref']}>{v['alt']}"
                  f"  [{kind}]  {v['filter']}")

    print(f"\n  --- the reference around it ---")
    seq = subprocess.run(
        f"{SAMTOOLS} faidx {REF} {chrom}:{pos-25}-{pos+25}",
        shell=True, capture_output=True, text=True).stdout
    seq = "".join(l.strip() for l in seq.splitlines()[1:]).upper()
    if seq:
        print(f"    {seq[:25]} [{seq[25]}] {seq[26:]}")
        tail = seq[26:40]
        run = 1
        while run < len(tail) and tail[run] == tail[0]:
            run += 1
        if run >= 4:
            print(f"    a run of {run} {tail[0]}s starts right after the site")
            print(f"    — a homopolymer, where indel calling is least certain")
    print()

print("=" * W)
print(" READING THIS")
print("=" * W)
print("""
  Nothing in the raw VCF means Mutect2 never considered the variant: the
  active region was not triggered, or assembly produced no haplotype
  carrying it. That is a detection limit, not a filtering decision.

  A record in the raw file with a filter tag in the filtered one means it
  was considered and rejected, and the tag names the reason.

  Indel carriers in the normal mean the variant is germline at that
  position — the injection landed on a site the backbone already varies
  at, and rejecting it is correct behaviour.
""")
