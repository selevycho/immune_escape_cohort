#!/usr/bin/env python3
"""
Can reads be assigned to one HLA allele or the other?

Step 9 removed reads from a locus without regard to which allele they came
from, so both copies lost the same share and the ratio between them stayed
near 0.6 whether 20% or 50% of the coverage was left. That is a loss of
coverage, not a loss of heterozygosity, and LOHHLA tests for the latter.

Simulating the real thing means removing reads from one allele only, which
requires knowing which allele each read came from. An earlier attempt gave
up: on panel data only 131 of 5082 reads could be assigned. But that
attempt used one IMGT sequence per allele, and the reference has since
been rebuilt with every subtype - the same change that took a sandbox
sample from 3 usable mismatch positions to 95.

This measures the assignable fraction as things now stand, before any
effort is spent on a simulation that may not be possible. Reads are
aligned to both alleles and assigned where the edit distance to one is
clearly lower; where the two are equally good the read is uninformative
and stays unassigned.

A locus needs enough assignable reads that removing a share of one
allele's reads produces a coverage difference larger than the sampling
noise. At 30x with 20% assignable, one allele has about 3 informative
reads - not enough. At 60% it has 9, which is workable.

Usage:
  python measure_allele_assignment.py [workspace] [--samples=O006,B018]
"""
import sys, os, subprocess, tempfile, shutil
from collections import Counter
import pandas as pd

args = [a for a in sys.argv[1:] if not a.startswith("--")]
WS = args[0] if args else os.environ.get(
    "WS", "/gpfs/bwfor/work/ws/fr_os136-immune_escape")

only = None
for a in sys.argv[1:]:
    if a.startswith("--samples="):
        only = set(a.split("=", 1)[1].split(","))

PANEL = f"{WS}/simulation/step9_lohhla_panel"
COHORT = f"{WS}/simulation/cohort"
DESIGN = f"{PANEL}/collections.tsv"
MANIFEST = f"{COHORT}/manifest.tsv"
LOHHLA = os.path.expanduser("~/immune_escape_project/soft/lohhla")
FASTA_ALL = f"{LOHHLA}/data/hla_all_lohhla.fasta"
OUT = os.path.expanduser("~/immune_escape_project/results")

WINDOWS = {
    "A": ("chr6", 29932000, 29956000),
    "B": ("chr6", 31343000, 31367000),
    "C": ("chr6", 31258000, 31282000),
}
MIN_MARGIN = 2      # edit distance difference needed to call a read


def run(cmd):
    return subprocess.run(cmd, shell=True, capture_output=True, text=True)


def build_ref(allele, out_fa):
    """One FASTA holding every IMGT subtype of a single allele."""
    n, write = 0, False
    with open(out_fa, "w") as fh:
        for line in open(FASTA_ALL):
            if line.startswith(">"):
                write = line[1:].strip().startswith(allele)
                if write:
                    n += 1
            if write:
                fh.write(line)
    return n


def best_nm(bam, ref_fa, tmp):
    """
    Smallest edit distance of each read against one allele's subtypes.

    A read may align to several subtypes of the same allele; the best of
    those is what matters, since they all represent the same allele.
    """
    run(f"bwa index {ref_fa} 2>/dev/null")
    sam = f"{tmp}/aln.sam"
    r = run(f"bwa mem -a -t 2 {ref_fa} {tmp}/reads.fq > {sam} 2>/dev/null")
    best = {}
    if not os.path.exists(sam):
        return best
    for line in open(sam):
        if line.startswith("@"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 12 or int(f[1]) & 4:
            continue
        nm = None
        for tag in f[11:]:
            if tag.startswith("NM:i:"):
                nm = int(tag[5:])
                break
        if nm is None:
            continue
        if f[0] not in best or nm < best[f[0]]:
            best[f[0]] = nm
    return best


man = pd.read_csv(MANIFEST, sep="\t")
if only:
    man = man[man.sample_id.isin(only)]

design = pd.read_csv(DESIGN, sep="\t") if os.path.exists(DESIGN) else pd.DataFrame()
target_of = {r["sample"]: r.target_locus
             for _, r in design[design.collection == "C_loss35"].iterrows()} \
            if len(design) else {}

print(f"checking {len(man)} samples\n", flush=True)
rows = []

for _, m in man.iterrows():
    sid = m.sample_id
    bam = f"{COHORT}/{sid}/{sid}_tumor.bam"
    hla = f"{COHORT}/{sid}/optitype/{sid}_result.tsv"
    if not (os.path.exists(bam) and os.path.exists(hla)):
        continue

    h = pd.read_csv(hla, sep="\t")
    if not len(h):
        continue
    x = h.iloc[0]

    tgt = str(target_of.get(sid, "")).replace("HLA-", "")
    loci = [tgt] if tgt in WINDOWS else list(WINDOWS)

    for loc in loci:
        a1, a2 = str(x.get(f"{loc}1", "")), str(x.get(f"{loc}2", ""))
        if "*" not in a1 or "*" not in a2 or a1 == a2:
            continue

        def imgt(a):
            g, rest = a.split("*")
            p = rest.split(":")
            return f"hla_{g.lower()}_{p[0]}_{p[1]}"

        k1, k2 = imgt(a1), imgt(a2)
        c, s, e = WINDOWS[loc]
        tmp = tempfile.mkdtemp(prefix=f"aa_{sid}_{loc}_")

        try:
            run(f"samtools view -b {bam} {c}:{s}-{e} 2>/dev/null | "
                f"samtools fastq - > {tmp}/reads.fq 2>/dev/null")
            n_reads = sum(1 for _ in open(f"{tmp}/reads.fq")) // 4
            if n_reads < 50:
                continue

            n1 = build_ref(k1, f"{tmp}/a1.fa")
            n2 = build_ref(k2, f"{tmp}/a2.fa")
            if n1 == 0 or n2 == 0:
                print(f"  {sid} HLA-{loc}: no IMGT sequences for "
                      f"{k1 if n1 == 0 else k2}")
                continue

            b1 = best_nm(bam, f"{tmp}/a1.fa", tmp)
            b2 = best_nm(bam, f"{tmp}/a2.fa", tmp)

            names = set(b1) | set(b2)
            to_1 = to_2 = tied = 0
            for nm in names:
                d1 = b1.get(nm, 999)
                d2 = b2.get(nm, 999)
                if d1 + MIN_MARGIN <= d2:
                    to_1 += 1
                elif d2 + MIN_MARGIN <= d1:
                    to_2 += 1
                else:
                    tied += 1

            assigned = to_1 + to_2
            rows.append({
                "sample": sid, "cohort": m.cohort, "locus": f"HLA-{loc}",
                "allele1": a1, "allele2": a2,
                "subtypes1": n1, "subtypes2": n2,
                "reads": n_reads, "aligned": len(names),
                "to_allele1": to_1, "to_allele2": to_2, "tied": tied,
                "assignable_pct": round(100 * assigned / max(1, len(names)), 1),
                "balance": round(min(to_1, to_2) / max(1, max(to_1, to_2)), 3),
            })
            print(f"  {sid} HLA-{loc}  {a1}/{a2}  "
                  f"{assigned}/{len(names)} assignable "
                  f"({100*assigned/max(1,len(names)):.0f}%), "
                  f"{to_1} vs {to_2}", flush=True)
        finally:
            shutil.rmtree(tmp, ignore_errors=True)

if not rows:
    print("nothing measured")
    sys.exit(1)

t = pd.DataFrame(rows)
pd.set_option("display.width", 220)

print()
print("=" * 92)
print(" PER LOCUS")
print("=" * 92)
print(t[["sample", "cohort", "locus", "allele1", "allele2", "reads",
         "aligned", "to_allele1", "to_allele2", "tied",
         "assignable_pct", "balance"]].to_string(index=False))

print()
print("=" * 92)
print(" SUMMARY")
print("=" * 92)
print(f"\n  loci measured        {len(t)}")
print(f"  assignable fraction  median {t.assignable_pct.median():.1f}%, "
      f"range {t.assignable_pct.min():.1f}-{t.assignable_pct.max():.1f}%")
print(f"  balance between the two alleles: median {t.balance.median():.2f}")
print(f"    (1.0 would mean equal numbers assigned each way, which is what")
print(f"     a heterozygous locus should give)")

print(f"\n  by locus:")
for loc, g in t.groupby("locus"):
    print(f"    {loc:<8}n={len(g):<3} assignable "
          f"{g.assignable_pct.median():.0f}%, "
          f"{g.to_allele1.median():.0f} vs {g.to_allele2.median():.0f} reads")

print()
print("=" * 92)
print(" IS AN ALLELE-SPECIFIC SIMULATION FEASIBLE?")
print("=" * 92)
t["informative_per_allele"] = t[["to_allele1", "to_allele2"]].min(axis=1)
workable = t[t.informative_per_allele >= 8]
print(f"\n  A locus needs enough reads assigned to the allele being removed")
print(f"  that taking a share of them changes coverage by more than the")
print(f"  sampling noise. Eight is a rough floor.\n")
print(f"  loci with 8 or more assignable reads on the weaker allele: "
      f"{len(workable)} of {len(t)}")
print(f"  median assignable per allele: "
      f"{t.informative_per_allele.median():.0f}")

if len(workable) >= 10:
    print(f"\n  Enough loci qualify to build the simulation.")
elif len(workable) > 0:
    print(f"\n  Only {len(workable)} loci qualify. A simulation on that many")
    print(f"  would not measure sensitivity to any useful precision.")
else:
    print(f"\n  No locus has enough assignable reads. The panel does not")
    print(f"  supply the information an allele-specific simulation needs,")
    print(f"  which is the same limit that stops LOHHLA itself.")

t.to_csv(f"{OUT}/allele_assignment.tsv", sep="\t", index=False)
print(f"\nwritten to {OUT}/allele_assignment.tsv")
