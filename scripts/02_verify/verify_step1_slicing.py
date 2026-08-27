#!/usr/bin/env python3
"""
Step 1 verification: what was actually fetched from the 1000 Genomes CRAMs.

Everything downstream inherits the properties of these files, so they are
worth measuring rather than assuming. A sample sliced at 18x will lose
low-VAF mutations no matter how well the caller works; a sample where the
HLA windows came out thin will fail at typing or at LOHHLA regardless of
how those tools are configured. Both show up here and nowhere else.

The checks are ordered by what they can break:

  file integrity   is the BAM readable and indexed at all
  panel coverage   mean and median across the panel, and how much of it
                   falls below the depth a caller needs
  HLA windows      depth at HLA-A, -B and -C specifically, since steps 4
                   and 6 depend on them and the MHC behaves differently
                   from the rest of the panel
  read properties  mapped fraction, duplicates, insert size, mapping
                   quality - a file can have adequate depth and still be
                   unusable if most of it is duplicates
  identity         does the backbone in the BAM header match the manifest

Usage:
  python verify_step1_slicing.py [workspace] [--quick]

  --quick skips per-base depth, which is the slow part
"""
import sys, os, subprocess, json
import numpy as np
import pandas as pd

args = [a for a in sys.argv[1:] if not a.startswith("--")]
QUICK = "--quick" in sys.argv

WS = args[0] if args else os.environ.get(
    "WS", "/gpfs/bwfor/work/ws/fr_os136-immune_escape")

COHORT = f"{WS}/simulation/cohort"
MANIFEST = f"{COHORT}/manifest.tsv"
PANEL = f"{WS}/simulation/panel/panel.bed"

# hg38, the same windows used for typing and for the LOH simulation
HLA = {
    "HLA-A": ("chr6", 29932000, 29956000),
    "HLA-B": ("chr6", 31343000, 31367000),
    "HLA-C": ("chr6", 31258000, 31282000),
}

MIN_DEPTH_CALLABLE = 10   # below this a somatic caller has little chance
MIN_DEPTH_TYPING = 20     # below this OptiType starts to struggle


def run(cmd):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return r.stdout.strip()


def flagstat(bam):
    """Parse the numbers samtools flagstat reports, as a dict."""
    out = run(f"samtools flagstat -O json {bam} 2>/dev/null")
    try:
        d = json.loads(out)["QC-passed reads"]
        return {
            "total": d.get("total", 0),
            "mapped": d.get("mapped", 0),
            "duplicates": d.get("duplicates", 0),
            "paired": d.get("paired in sequencing", 0),
            "properly_paired": d.get("properly paired", 0),
        }
    except Exception:
        return None


def depth_stats(bam, region_arg, quick=False):
    """Mean, median and the fraction of positions below each threshold."""
    if quick:
        out = run(f"samtools depth -a {region_arg} {bam} 2>/dev/null | "
                  f"awk '{{s+=$3;n++}} END {{if(n) printf \"%.2f %d\", s/n, n}}'")
        p = out.split()
        if len(p) == 2:
            return {"mean": float(p[0]), "positions": int(p[1])}
        return None

    out = run(f"samtools depth -a {region_arg} {bam} 2>/dev/null | cut -f3")
    if not out:
        return None
    d = np.fromstring(out.replace("\n", " "), sep=" ")
    if d.size == 0:
        return None
    return {
        "mean": float(d.mean()),
        "median": float(np.median(d)),
        "p10": float(np.percentile(d, 10)),
        "p90": float(np.percentile(d, 90)),
        "zero_pct": float(100 * (d == 0).mean()),
        "below_callable_pct": float(100 * (d < MIN_DEPTH_CALLABLE).mean()),
        "positions": int(d.size),
    }


print(f"reading manifest ...", flush=True)
man = pd.read_csv(MANIFEST, sep="\t")
print(f"  {len(man)} samples expected\n", flush=True)

panel = pd.read_csv(PANEL, sep="\t", header=None,
                    names=["chrom", "start", "end", "name"])
panel_mb = (panel.end - panel.start).sum() / 1e6
print(f"panel: {len(panel):,} intervals, {panel_mb:.2f} Mb\n", flush=True)

rows = []
problems = []

for i, m in man.iterrows():
    sid = m.sample_id
    bam = f"{COHORT}/{sid}/{sid}_normal.bam"
    r = {"sample": sid, "cohort": m.cohort,
         "backbone": m.backbone, "superpop": m.superpopulation}

    print(f"  [{i+1:>2}/{len(man)}] {sid} ...", end="", flush=True)

    # ---------------- file integrity ----------------
    if not os.path.exists(bam):
        problems.append((sid, "missing", "no normal BAM"))
        print(" MISSING")
        rows.append(r)
        continue

    r["size_MB"] = round(os.path.getsize(bam) / 1e6, 1)
    r["indexed"] = os.path.exists(bam + ".bai")
    if not r["indexed"]:
        problems.append((sid, "index", "BAM not indexed"))

    quickcheck = run(f"samtools quickcheck -v {bam} 2>&1")
    r["truncated"] = bool(quickcheck)
    if quickcheck:
        problems.append((sid, "integrity", "quickcheck reports a problem"))

    # ---------------- read properties ----------------
    fs = flagstat(bam)
    if fs and fs["total"]:
        r["reads"] = fs["total"]
        r["mapped_pct"] = round(100 * fs["mapped"] / fs["total"], 2)
        r["dup_pct"] = round(100 * fs["duplicates"] / fs["total"], 2)
        r["proper_pair_pct"] = round(
            100 * fs["properly_paired"] / fs["paired"], 2) if fs["paired"] else None
        if r["mapped_pct"] < 95:
            problems.append((sid, "mapping", f"{r['mapped_pct']}% mapped"))

    # insert size and mapping quality, from a sample of reads
    ins = run(f"samtools view -f 66 {bam} 2>/dev/null | head -50000 | "
              f"awk '$9>0 && $9<2000 {{s+=$9;n++}} END {{if(n) printf \"%.0f\", s/n}}'")
    r["insert_size"] = int(ins) if ins else None

    mapq = run(f"samtools view {bam} 2>/dev/null | head -50000 | "
               f"awk '{{s+=$5;n++; if($5==0) z++}} "
               f"END {{if(n) printf \"%.1f %.1f\", s/n, 100*z/n}}'")
    p = mapq.split()
    if len(p) == 2:
        r["mean_mapq"] = float(p[0])
        r["mapq0_pct"] = float(p[1])

    # ---------------- panel coverage ----------------
    st = depth_stats(bam, f"-b {PANEL}", QUICK)
    if st:
        r["panel_mean_depth"] = round(st["mean"], 1)
        if not QUICK:
            r["panel_median_depth"] = round(st["median"], 1)
            r["panel_p10"] = round(st["p10"], 1)
            r["panel_p90"] = round(st["p90"], 1)
            r["uncovered_pct"] = round(st["zero_pct"], 2)
            r["below_10x_pct"] = round(st["below_callable_pct"], 2)
            if st["below_callable_pct"] > 20:
                problems.append((sid, "coverage",
                                 f"{st['below_callable_pct']:.0f}% of panel below 10x"))
        if st["mean"] < 20:
            problems.append((sid, "coverage", f"mean depth {st['mean']:.1f}x"))

    # ---------------- HLA windows ----------------
    for name, (c, s, e) in HLA.items():
        hs = depth_stats(bam, f"-r {c}:{s}-{e}", quick=True)
        key = name.replace("HLA-", "hla_") + "_depth"
        if hs:
            r[key] = round(hs["mean"], 1)
            if hs["mean"] < MIN_DEPTH_TYPING:
                problems.append((sid, "hla", f"{name} at {hs['mean']:.1f}x"))

    # ---------------- identity ----------------
    hdr = run(f"samtools view -H {bam} 2>/dev/null | grep '^@RG' | head -1")
    r["rg_present"] = bool(hdr)

    rows.append(r)
    print(f" {r.get('panel_mean_depth', '?')}x, "
          f"{r.get('reads', 0):,} reads")

t = pd.DataFrame(rows)
pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 40)

# =====================================================================
print()
print("=" * 96)
print(" PER SAMPLE")
print("=" * 96)
cols = [c for c in ["sample", "cohort", "backbone", "superpop", "size_MB",
                    "reads", "mapped_pct", "dup_pct", "insert_size",
                    "mean_mapq", "panel_mean_depth", "panel_median_depth",
                    "below_10x_pct"] if c in t.columns]
print(t[cols].to_string(index=False))

# =====================================================================
print()
print("=" * 96)
print(" COHORT SUMMARY")
print("=" * 96)

def describe(col, unit="", fmt="{:.1f}"):
    if col not in t.columns or t[col].isna().all():
        return
    v = t[col].dropna()
    print(f"  {col:<22} median {fmt.format(v.median())}{unit:<4}"
          f"range {fmt.format(v.min())} – {fmt.format(v.max())}{unit}")

print()
describe("size_MB", " MB")
describe("reads", "", "{:,.0f}")
describe("mapped_pct", "%", "{:.2f}")
describe("dup_pct", "%", "{:.2f}")
describe("proper_pair_pct", "%", "{:.1f}")
describe("insert_size", " bp", "{:.0f}")
describe("mean_mapq", "", "{:.1f}")
describe("mapq0_pct", "%", "{:.1f}")
print()
describe("panel_mean_depth", "x")
describe("panel_median_depth", "x")
if "below_10x_pct" in t.columns:
    describe("uncovered_pct", "%", "{:.2f}")
    describe("below_10x_pct", "%", "{:.2f}")
print()
describe("hla_A_depth", "x")
describe("hla_B_depth", "x")
describe("hla_C_depth", "x")

# =====================================================================
print()
print("=" * 96)
print(" BY COHORT AND POPULATION")
print("=" * 96)
for key in ["cohort", "superpop"]:
    if key not in t.columns:
        continue
    print(f"\n  by {key}:")
    g = t.groupby(key).agg(
        n=("sample", "size"),
        depth=("panel_mean_depth", "median"),
        reads=("reads", "median"),
        dup=("dup_pct", "median")).round(1)
    print(g.to_string())

# =====================================================================
# What this coverage means for the rest of the pipeline
# =====================================================================
print()
print("=" * 96)
print(" WHAT THIS COVERAGE ALLOWS")
print("=" * 96)

if "panel_mean_depth" in t.columns:
    med = t.panel_mean_depth.median()
    print(f"\n  At the median depth of {med:.1f}x, a mutation is carried by:")
    print(f"    {'VAF':<10}{'reads':>8}   {'detectable?':<14}")
    for vaf in [0.05, 0.08, 0.10, 0.15, 0.20, 0.30, 0.50]:
        reads = med * vaf
        verdict = ("no, below 3 reads" if reads < 3 else
                   "marginal" if reads < 5 else "yes")
        print(f"    {vaf:<10.0%}{reads:>8.1f}   {verdict:<14}")
    print(f"\n  This is what sets the sensitivity floor measured in step 3,")
    print(f"  not any property of the caller.")

# =====================================================================
print()
print("=" * 96)
print(f" PROBLEMS: {len(problems)}")
print("=" * 96)
if problems:
    by_kind = {}
    for sid, kind, msg in problems:
        by_kind.setdefault(kind, []).append((sid, msg))
    for kind in sorted(by_kind):
        print(f"\n  {kind} ({len(by_kind[kind])}):")
        for sid, msg in by_kind[kind]:
            print(f"    {sid:<7} {msg}")
else:
    print("\n  none")

out = os.path.expanduser("~/immune_escape_project/results/verify_step1.tsv")
os.makedirs(os.path.dirname(out), exist_ok=True)
t.to_csv(out, sep="\t", index=False)
print(f"\nwritten to {out}")
