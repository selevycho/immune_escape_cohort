#!/usr/bin/env python3
"""
Step 6a verification: was the loss simulated as designed?

The simulation removes reads at one class I locus until a known fraction
of the coverage remains, and leaves the other two loci alone. Both halves
matter equally. The thinned locus is what LOHHLA is asked to find; the
untouched ones are what tells us a detection means anything, since a tool
that reports loss where none was made is worse than one that misses it.

This is arithmetic rather than inference, and it is worth checking
directly because the thinning is done by a random subsample: samtools
view -s draws reads independently, so the realised fraction will not
match the requested one exactly, and a systematic gap between them would
mean the design is not what the results assume.

Checked here:

  target        observed depth ratio against the fraction requested
  controls      did the untouched loci move at all
  choice        was the target genuinely heterozygous, and does the
                distribution of chosen loci follow from the genotypes
  files         are the BAMs readable, indexed, and paired with a normal

Usage:
  python verify_step6a_loh.py [workspace]
"""
import sys, os, subprocess
import numpy as np
import pandas as pd

args = [a for a in sys.argv[1:] if not a.startswith("--")]
WS = args[0] if args else os.environ.get(
    "WS", "/gpfs/bwfor/work/ws/fr_os136-immune_escape")

COHORT = f"{WS}/simulation/cohort"
LOH = f"{WS}/simulation/lohhla_allelic"
MANIFEST = f"{COHORT}/manifest.tsv"
DESIGN = f"{LOH}/loh_design.tsv"

HLA = {
    "A": ("chr6", 29932000, 29956000),
    "B": ("chr6", 31343000, 31367000),
    "C": ("chr6", 31258000, 31282000),
}
TOLERANCE = 0.05      # how far the realised fraction may drift
CONTROL_DRIFT = 0.5   # depth change that counts as a control having moved


def depth(bam, chrom, start, end):
    cmd = (f"samtools depth -a -r {chrom}:{start}-{end} {bam} 2>/dev/null | "
           f"awk '{{s+=$3;n++}} END {{if(n) printf \"%.2f\", s/n}}'")
    out = subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout
    try:
        return float(out.strip())
    except ValueError:
        return None


man = pd.read_csv(MANIFEST, sep="\t")
design = pd.read_csv(DESIGN, sep="\t") if os.path.exists(DESIGN) else pd.DataFrame()

if design.empty:
    print("no loh_design.tsv - step 6a has not run")
    sys.exit(1)

target_of = dict(zip(design["sample"], design.target_locus))
keep_of = dict(zip(design["sample"], design.keep_fraction)) \
    if "keep_fraction" in design.columns else {}

print(f"checking {len(man)} samples against the design\n", flush=True)

rows, problems = [], []

for i, m in man.iterrows():
    sid = m.sample_id
    loh_bam = f"{LOH}/{sid}/{sid}_tumor_LOH.bam"
    src_bam = f"{LOH}/{sid}/{sid}_tumor_original.bam"
    if not os.path.exists(src_bam):
        src_bam = f"{COHORT}/{sid}/{sid}_tumor.bam"
    norm_bam = f"{LOH}/{sid}/{sid}_normal.bam"

    r = {"sample": sid, "cohort": m.cohort,
         "target": target_of.get(sid),
         "requested": keep_of.get(sid)}

    print(f"  [{i+1:>2}/{len(man)}] {sid} ...", end="", flush=True)

    if not os.path.exists(loh_bam):
        problems.append((sid, "missing", "no LOH BAM"))
        print(" MISSING")
        rows.append(r)
        continue

    r["indexed"] = os.path.exists(loh_bam + ".bai")
    r["has_normal"] = os.path.exists(norm_bam)
    r["size_MB"] = round(os.path.getsize(loh_bam) / 1e6, 1)
    if not r["indexed"]:
        problems.append((sid, "index", "LOH BAM not indexed"))
    if not r["has_normal"]:
        problems.append((sid, "normal", "no paired normal in the LOH directory"))

    tgt = str(r["target"]).replace("HLA-", "") if r["target"] else None

    for loc, (c, s, e) in HLA.items():
        before = depth(src_bam, c, s, e)
        after = depth(loh_bam, c, s, e)
        r[f"{loc}_before"] = before
        r[f"{loc}_after"] = after
        if before and after is not None:
            ratio = after / before if before else None
            r[f"{loc}_ratio"] = round(ratio, 3) if ratio is not None else None

            if loc == tgt:
                want = r["requested"] or 0.35
                if ratio is not None and abs(ratio - want) > TOLERANCE:
                    problems.append((sid, "ratio",
                                     f"HLA-{loc} kept {ratio:.2f}, asked {want:.2f}"))
            else:
                if abs(after - before) > CONTROL_DRIFT:
                    problems.append((sid, "control",
                                     f"HLA-{loc} moved {before:.1f} -> {after:.1f}"))

    # was the target actually heterozygous
    hla = f"{COHORT}/{sid}/optitype/{sid}_result.tsv"
    if os.path.exists(hla) and tgt:
        h = pd.read_csv(hla, sep="\t")
        if len(h):
            x = h.iloc[0]
            a, b = str(x.get(f"{tgt}1", "")), str(x.get(f"{tgt}2", ""))
            r["target_alleles"] = f"{a}/{b}"
            r["target_het"] = ("*" in a and "*" in b and a != b)
            if not r["target_het"]:
                problems.append((sid, "homozygous",
                                 f"target HLA-{tgt} is {a}/{b}"))

            het = ""
            for loc in "ABC":
                p, q = str(x.get(f"{loc}1", "")), str(x.get(f"{loc}2", ""))
                if "*" in p and "*" in q and p != q:
                    het += loc
            r["het_loci"] = het or "-"
            if het and tgt != het[0]:
                problems.append((sid, "choice",
                                 f"target {tgt}, first heterozygous is {het[0]}"))

    rows.append(r)
    rr = r.get(f"{tgt}_ratio") if tgt else None
    print(f" HLA-{tgt} {rr if rr is not None else '?'}")

t = pd.DataFrame(rows)
pd.set_option("display.width", 240)
pd.set_option("display.max_columns", 40)

# =====================================================================
print()
print("=" * 100)
print(" PER SAMPLE")
print("=" * 100)
print(f"\n  {'sample':<7}{'target':<8}{'alleles':<20}"
      f"{'before':>9}{'after':>8}{'ratio':>8}   controls")
for _, r in t.iterrows():
    tgt = str(r.target).replace("HLA-", "") if pd.notna(r.target) else None
    if not tgt:
        continue
    b = r.get(f"{tgt}_before")
    a = r.get(f"{tgt}_after")
    ratio = r.get(f"{tgt}_ratio")
    ctrl = []
    for loc in "ABC":
        if loc == tgt:
            continue
        cb, ca = r.get(f"{loc}_before"), r.get(f"{loc}_after")
        if cb is None or ca is None:
            continue
        moved = abs(ca - cb) > CONTROL_DRIFT
        ctrl.append(f"{loc} {ca:.1f}{' MOVED' if moved else ''}")
    print(f"  {r['sample']:<7}{str(r.target):<8}{str(r.get('target_alleles','')):<20}"
          f"{b if b else 0:>9.1f}{a if a else 0:>8.1f}"
          f"{ratio if ratio else 0:>8.2f}   {', '.join(ctrl)}")

# =====================================================================
print()
print("=" * 100)
print(" TARGET THINNING")
print("=" * 100)
ratios = []
for _, r in t.iterrows():
    tgt = str(r.target).replace("HLA-", "") if pd.notna(r.target) else None
    v = r.get(f"{tgt}_ratio") if tgt else None
    if v is not None and not pd.isna(v):
        ratios.append(v)

if ratios:
    ratios = np.array(ratios)
    want = np.median([v for v in t.requested.dropna()]) if t.requested.notna().any() else 0.35
    print(f"\n  requested fraction    {want:.2f}")
    print(f"  realised              median {np.median(ratios):.3f}, "
          f"range {ratios.min():.3f} – {ratios.max():.3f}")
    print(f"  deviation             mean {np.mean(ratios) - want:+.4f}, "
          f"sd {ratios.std():.4f}")
    print(f"  within {TOLERANCE:.2f}         "
          f"{int((abs(ratios - want) <= TOLERANCE).sum())} of {len(ratios)}")
    print(f"\n  The spread comes from samtools drawing reads independently;")
    print(f"  a systematic offset would mean the thinning is not doing")
    print(f"  what the design records.")

# =====================================================================
print()
print("=" * 100)
print(" CONTROL LOCI")
print("=" * 100)
moved = 0
checked = 0
for _, r in t.iterrows():
    tgt = str(r.target).replace("HLA-", "") if pd.notna(r.target) else None
    for loc in "ABC":
        if loc == tgt:
            continue
        b, a = r.get(f"{loc}_before"), r.get(f"{loc}_after")
        if b is None or a is None:
            continue
        checked += 1
        if abs(a - b) > CONTROL_DRIFT:
            moved += 1
print(f"\n  control loci checked  {checked}")
print(f"  depth unchanged       {checked - moved}")
print(f"  moved by >{CONTROL_DRIFT}x         {moved}")
if moved == 0:
    print(f"\n  No control locus moved. Any loss LOHHLA reports at one of")
    print(f"  them is a false positive with no ambiguity about the input.")

# =====================================================================
print()
print("=" * 100)
print(" WHICH LOCUS WAS CHOSEN")
print("=" * 100)
print(f"\n  {t.target.value_counts().to_string()}")
if "het_loci" in t.columns:
    print(f"\n  heterozygous loci per sample:")
    print(f"  {t.het_loci.value_counts().to_string()}")
    print(f"\n  The target is the first heterozygous locus in the order")
    print(f"  A, B, C, so HLA-A dominates because it is heterozygous most")
    print(f"  often, not because it was preferred.")

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
        for sid, msg in by_kind[kind][:15]:
            print(f"    {sid:<7} {msg}")
else:
    print("\n  none")

out = os.path.expanduser("~/immune_escape_project/results/verify_step6a.tsv")
os.makedirs(os.path.dirname(out), exist_ok=True)
t.to_csv(out, sep="\t", index=False)
print(f"\nwritten to {out}")
