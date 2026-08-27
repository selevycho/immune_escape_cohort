#!/usr/bin/env python3
"""
Line 3 refined: separate FOCAL from ARM-LEVEL loss at antigen-presentation genes.

Logic (after McGranahan et al. 2017):
  - arm-level loss = the lost region containing the gene spans most of the arm
                     -> a passenger consequence of chromosomal chaos
  - focal loss     = a narrow lost region centred near the gene
                     -> cannot happen by accident; evidence of selection

We compare the observed FOCAL loss rate at each target gene against each
patient's own genome-wide focal loss background.
"""
import sys, os
import pandas as pd
import numpy as np
from scipy import stats

SEG_FILE = sys.argv[1]           # .../data_cna_hg19.seg
OUT_DIR  = sys.argv[2]
COHORT   = sys.argv[3]

LOSS_THR   = -0.20     # log2 ratio below this = loss
ARM_FRAC   = 0.50      # lost block >= 50% of arm  -> arm-level, else focal

# hg19 gene coordinates (TCGA .seg files are hg19)
GENES = {
    "HLA-A": ("6", 29910247, 29913661),
    "HLA-B": ("6", 31321649, 31324989),
    "HLA-C": ("6", 31236526, 31239913),
    "TAP1" : ("6", 32812986, 32821755),
    "TAP2" : ("6", 32789610, 32806547),
    "B2M"  : ("15", 45003675, 45011075),
    "NLRC5": ("16", 57023472, 57113451),
    # negative controls: neighbours, same arms, no immune role
    "TRIM69": ("15", 45028545, 45050010),
    "SPG11" : ("15", 44854781, 44955876),
    "CYLD"  : ("16", 50775961, 50835846),
    "NKD1"  : ("16", 50583491, 50666687),
}

# hg19 centromere positions and chromosome lengths
CENTRO = {"6": 61000000, "15": 20000000, "16": 38300000}
CHRLEN = {"6": 171115067, "15": 102531392, "16": 90354753}

def arm_of(chrom, start):
    c = CENTRO[chrom]
    return (chrom + "p", 0, c) if start < c else (chrom + "q", c, CHRLEN[chrom])

def norm_id(s):
    p = str(s).split("-")
    return "-".join(p[:4])[:15] if len(p) >= 4 else str(s)

os.makedirs(OUT_DIR, exist_ok=True)

print(f"[1/5] reading {SEG_FILE} ...", flush=True)
seg = pd.read_csv(SEG_FILE, sep="\t", low_memory=False)
seg.columns = [c.strip() for c in seg.columns]
idcol   = seg.columns[0]
chrcol  = [c for c in seg.columns if "chrom" in c.lower()][0]
startc  = [c for c in seg.columns if "start" in c.lower()][0]
endc    = [c for c in seg.columns if "end"   in c.lower()][0]
meanc   = [c for c in seg.columns if "mean"  in c.lower()][0]

seg = seg[[idcol, chrcol, startc, endc, meanc]]
seg.columns = ["sample", "chrom", "start", "end", "logr"]
seg["sample"] = seg["sample"].map(norm_id)
seg["chrom"]  = seg["chrom"].astype(str).str.replace("chr", "", regex=False)
seg = seg.dropna(subset=["logr"])
seg["len"] = seg["end"] - seg["start"]
samples = sorted(seg["sample"].unique())
print(f"      samples={len(samples)}  segments={len(seg)}", flush=True)

print("[2/5] computing per-sample focal background ...", flush=True)
# a lost segment is "focal" if shorter than ARM_FRAC of the arm it sits on
def is_focal_seg(row):
    ch = row["chrom"]
    if ch not in CENTRO:
        # approximate arms for other chromosomes using segment context
        return row["len"] < 50e6
    _, a0, a1 = arm_of(ch, row["start"])
    return row["len"] < ARM_FRAC * (a1 - a0)

lost = seg[seg["logr"] < LOSS_THR].copy()
lost["focal"] = lost.apply(is_focal_seg, axis=1)

total_bp = seg.groupby("sample")["len"].sum()
focal_bp = lost[lost["focal"]].groupby("sample")["len"].sum().reindex(total_bp.index, fill_value=0)
arm_bp   = lost[~lost["focal"]].groupby("sample")["len"].sum().reindex(total_bp.index, fill_value=0)

bg_focal = (focal_bp / total_bp).fillna(0)
bg_arm   = (arm_bp   / total_bp).fillna(0)
print(f"      genome-wide FOCAL loss fraction : median={bg_focal.median():.3f}", flush=True)
print(f"      genome-wide ARM   loss fraction : median={bg_arm.median():.3f}", flush=True)

print("[3/5] merging contiguous lost blocks per sample ...", flush=True)
# build, per sample+chrom, merged blocks of consecutive lost segments
blocks = {}
for (smp, ch), g in lost.groupby(["sample", "chrom"], sort=False):
    g = g.sort_values("start")
    cur_s, cur_e = None, None
    out = []
    for s, e in zip(g["start"].values, g["end"].values):
        if cur_s is None:
            cur_s, cur_e = s, e
        elif s - cur_e <= 1e6:          # allow 1 Mb gap
            cur_e = max(cur_e, e)
        else:
            out.append((cur_s, cur_e)); cur_s, cur_e = s, e
    if cur_s is not None:
        out.append((cur_s, cur_e))
    blocks[(smp, ch)] = out

print("[4/5] classifying loss at each target gene ...", flush=True)
rows = []
per_gene_calls = {}
for gene, (ch, gs, ge) in GENES.items():
    _, a0, a1 = arm_of(ch, gs)
    arm_len = a1 - a0
    calls = {}
    for smp in samples:
        blist = blocks.get((smp, ch), [])
        hit = None
        for (bs, be) in blist:
            if bs <= ge and be >= gs:
                hit = (bs, be); break
        if hit is None:
            calls[smp] = "none"
        else:
            blen = hit[1] - hit[0]
            calls[smp] = "arm" if blen >= ARM_FRAC * arm_len else "focal"
    per_gene_calls[gene] = calls
    v = pd.Series(calls)
    n = len(v)
    n_focal = int((v == "focal").sum())
    n_arm   = int((v == "arm").sum())
    exp_focal = bg_focal.reindex(v.index).mean()
    exp_arm   = bg_arm.reindex(v.index).mean()
    obs_focal = n_focal / n
    obs_arm   = n_arm / n
    # binomial test of focal loss against that gene's expected focal rate
    p_focal = stats.binomtest(n_focal, n, min(max(exp_focal, 1e-6), 0.999),
                              alternative="greater").pvalue
    rows.append({
        "gene": gene, "arm": ch + ("p" if gs < CENTRO[ch] else "q"),
        "n": n,
        "n_focal": n_focal, "obs_focal": obs_focal, "exp_focal": exp_focal,
        "enrich_focal": obs_focal / exp_focal if exp_focal > 0 else np.nan,
        "p_focal": p_focal,
        "n_arm": n_arm, "obs_arm": obs_arm, "exp_arm": exp_arm,
        "enrich_arm": obs_arm / exp_arm if exp_arm > 0 else np.nan,
    })

res = pd.DataFrame(rows).sort_values("enrich_focal", ascending=False)
res["p_focal_bonf"] = (res.p_focal * len(res)).clip(upper=1.0)
res["cohort"] = COHORT

print("[5/5] writing ...", flush=True)
res.to_csv(f"{OUT_DIR}/{COHORT}_line3_focal_vs_arm.tsv", sep="\t", index=False)
calls_df = pd.DataFrame(per_gene_calls)
calls_df.index.name = "sample"
calls_df["bg_focal"] = bg_focal.reindex(calls_df.index)
calls_df["cohort"] = COHORT
calls_df.to_csv(f"{OUT_DIR}/{COHORT}_line3_focal_calls.tsv", sep="\t")

pd.set_option("display.width", 200)
print()
print(f"=== {COHORT.upper()} — Line 3: focal vs arm-level ===")
print(res[["gene","arm","n","n_focal","obs_focal","exp_focal","enrich_focal",
           "p_focal_bonf","n_arm","obs_arm","enrich_arm"]]
      .to_string(index=False, float_format=lambda x: f"{x:.3f}"))
print()
print("How to read this:")
print("  enrich_arm   ~ 1  -> arm-level loss is just background chromosomal chaos")
print("  enrich_focal > 1 with p_focal_bonf < 0.05  -> narrow, targeted loss = SELECTION")
print("  the neighbour genes (TRIM69, SPG11, CYLD, NKD1) are negative controls:")
print("  they should NOT show focal enrichment. If they do, the signal is regional.")
