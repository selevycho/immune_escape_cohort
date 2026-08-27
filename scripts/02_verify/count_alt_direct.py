#!/usr/bin/env python3
"""
Count alternate-supporting reads without mpileup.

samtools 1.23 mpileup silently requires the proper-pair flag and offers no
option to disable it. That is fatal here: BAMSurgeon edits reads and
realigns them through bwa, and realignment strips proper-pair status from
exactly the reads it rewrote. mpileup therefore discards the
mutation-carrying reads preferentially, and every position looks like a
failed injection.

At one checked position, four of forty reads carry the alternate base and
none of those four is flagged proper-pair. mpileup reported depth 1.

So bases are read from the alignment directly. The read's own CIGAR gives
the offset from its start to the position of interest; soft clips consume
query bases without consuming reference, insertions consume query only,
deletions and skips consume reference only.

Usage:
  python count_alt_direct.py <bam> <chrom> <pos> <alt>
  python count_alt_direct.py --check <workspace>
"""
import re
import subprocess
import sys

CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def base_at(pos_read, cigar, seq, target):
    """
    The query base aligned to a reference position, or None.

    Walking the CIGAR is necessary rather than subtracting coordinates:
    a soft clip shifts the query without shifting the reference, so a
    naive offset lands on the wrong base whenever a read is clipped.
    """
    ref, qry = pos_read, 0
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
                return "*"          # the position is deleted in this read
            ref += n
        elif op == "H":
            pass
    return None


def count(bam, chrom, pos, alt, samtools="samtools", exclude_dup=True):
    flags = "-F 0x400" if exclude_dup else ""
    out = subprocess.run(
        f"{samtools} view {flags} {bam} {chrom}:{pos}-{pos}",
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
    return depth, n_alt, (n_alt / depth if depth else 0.0)


if __name__ == "__main__":
    if sys.argv[1] == "--check":
        import os
        import pandas as pd
        WS = sys.argv[2]
        RES = os.path.expanduser("~/immune_escape_project/results")
        v = pd.read_csv(f"{RES}/verify_step2_per_mutation.tsv", sep="\t")
        v["alt"] = v.alt.astype(str)
        v = v[(v.alt.str.len() == 1) & (v.alt != "-")]
        lost = v[~v.landed].sample(min(25, int((~v.landed).sum())),
                                   random_state=3)
        print(" positions mpileup called absent, re-counted directly\n")
        print(f"  {'sample':<8}{'position':<24}{'alt':<5}"
              f"{'depth':>7}{'alt reads':>11}{'VAF':>8}")
        rescued = 0
        for _, r in lost.iterrows():
            bam = f"{WS}/simulation/cohort/{r['sample']}/{r['sample']}_tumor.bam"
            d, a, f_ = count(bam, r.chrom, int(r.pos), r.alt)
            if a > 0:
                rescued += 1
            print(f"  {r['sample']:<8}{r.chrom + ':' + str(int(r.pos)):<24}"
                  f"{r.alt:<5}{d:>7}{a:>11}{f_:>8.3f}")
        print(f"\n  {rescued} of {len(lost)} carry the alternate base after all")
    else:
        bam, chrom, pos, alt = sys.argv[1:5]
        d, a, f_ = count(bam, chrom, int(pos), alt)
        print(f"depth {d}, alt {a}, VAF {f_:.3f}")
