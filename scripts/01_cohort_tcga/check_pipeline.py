#!/usr/bin/env python3
"""
Audit every stage of the cohort pipeline against the manifest.

Presence of a file is not evidence that a stage worked. OptiType writes a
default-looking genotype when razers3 dies, LOHHLA writes a "no analysis
possible" file that looks like a result until opened, and pVACseq happily
produces a header-only TSV. This script therefore opens each output and
checks its contents, not just its existence.

Per sample it reports:
  BAM        normal and tumour present, indexed, with plausible read counts
  truth      injected mutations and how many survived verification
  VCF        Mutect2 calls, PASS count, recall against the truth set
  HLA        genotype, whether it looks like a real call
  mhcflurry  peptides and strong binders
  LOH        which locus was thinned and by how much
  LOHHLA     whether a real prediction came back, and what it says
  pVACseq    epitopes before and after filtering

Anything that fails a check is listed at the end with the reason.

Usage:
  python check_pipeline.py [workspace]
"""
import sys, os, glob, gzip
import pandas as pd

WS = sys.argv[1] if len(sys.argv) > 1 else os.environ.get(
    "WS", "/gpfs/bwfor/work/ws/fr_os136-immune_escape")

COHORT = f"{WS}/simulation/cohort"
LOH = f"{WS}/simulation/lohhla_allelic"
PVAC = f"{WS}/simulation/pvacseq/cohort"
MANIFEST = f"{COHORT}/manifest.tsv"

# the commonest European alleles - a genotype made only of these across
# many samples is the signature of OptiType falling back to defaults
SUSPICIOUS = {"A*02:01", "B*07:02", "C*07:02"}

problems = []


def note(sid, stage, msg):
    problems.append((sid, stage, msg))


def count_lines(path, gz=False):
    try:
        op = gzip.open if gz else open
        with op(path, "rt") as fh:
            return sum(1 for _ in fh)
    except Exception:
        return 0


def vcf_counts(path):
    """Return (records, PASS records) without loading the whole file."""
    n = p = 0
    try:
        with gzip.open(path, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                n += 1
                if line.split("\t")[6] == "PASS":
                    p += 1
    except Exception:
        pass
    return n, p


print("reading manifest ...", flush=True)
man = pd.read_csv(MANIFEST, sep="\t")
print(f"  {len(man)} samples\n", flush=True)

rows = []
genotypes = []

for _, m in man.iterrows():
    sid = m.sample_id
    d = f"{COHORT}/{sid}"
    r = {"sample": sid, "cohort": m.cohort, "expected_mut": m.n_panel_mut}

    # ---------- BAMs ----------
    for kind in ["normal", "tumor"]:
        bam = f"{d}/{sid}_{kind}.bam"
        ok = os.path.exists(bam) and os.path.exists(bam + ".bai")
        r[f"{kind}_bam"] = "ok" if ok else "MISSING"
        if ok:
            size_mb = os.path.getsize(bam) / 1e6
            r[f"{kind}_MB"] = round(size_mb, 0)
            if size_mb < 50:
                note(sid, "bam", f"{kind} only {size_mb:.0f} MB")
        else:
            note(sid, "bam", f"{kind} missing or unindexed")

    # ---------- truth set ----------
    truth = f"{d}/truth_set.tsv"
    if os.path.exists(truth):
        t = pd.read_csv(truth, sep="\t")
        r["injected"] = len(t)
        r["median_VAF"] = round(t.VAF.median(), 3) if "VAF" in t else None
        if len(t) < 5:
            note(sid, "truth", f"only {len(t)} mutations")
    else:
        r["injected"] = 0
        note(sid, "truth", "truth_set.tsv missing")

    # ---------- verification ----------
    ver = f"{d}/verification/injection_verification.tsv"
    if os.path.exists(ver):
        v = pd.read_csv(ver, sep="\t")
        usable = int((v.iloc[:, -1] == "OK").sum()) if len(v) else 0
        r["usable"] = usable
        if r["injected"]:
            pct = 100 * usable / r["injected"]
            r["usable_%"] = round(pct, 1)
            if pct < 70:
                note(sid, "verify", f"only {pct:.0f}% usable")

    # ---------- Mutect2 ----------
    vcf = f"{d}/mutect2/{sid}.filtered.vcf.gz"
    if os.path.exists(vcf):
        n, p = vcf_counts(vcf)
        r["vcf_records"] = n
        r["vcf_PASS"] = p
        if n == 0:
            note(sid, "vcf", "no records")
        elif p == 0:
            note(sid, "vcf", "no PASS calls")
    else:
        r["vcf_records"] = 0
        note(sid, "vcf", "missing")

    # recall against the truth
    comp = f"{d}/comparison/truth_vs_calls.tsv"
    if os.path.exists(comp):
        c = pd.read_csv(comp, sep="\t")
        if "detected_pass" in c and len(c):
            rec = 100 * c.detected_pass.mean()
            r["recall_%"] = round(rec, 1)
            if rec < 40:
                note(sid, "recall", f"{rec:.0f}%")

    # ---------- OptiType ----------
    hla = f"{d}/optitype/{sid}_result.tsv"
    if os.path.exists(hla):
        h = pd.read_csv(hla, sep="\t")
        if len(h):
            row = h.iloc[0]
            alleles = [str(row.get(k, "")) for k in ["A1", "A2", "B1", "B2", "C1", "C2"]]
            typed = [a for a in alleles if "*" in a]
            r["hla"] = " ".join(typed)
            r["n_alleles"] = len(typed)
            r["hla_reads"] = row.get("Reads")
            genotypes.append((sid, tuple(alleles)))
            if len(typed) < 4:
                note(sid, "hla", f"only {len(typed)} alleles typed")
            if row.get("Reads", 0) and float(row["Reads"]) < 100:
                note(sid, "hla", f"only {row['Reads']} reads used")
            het = "".join(
                x for x, (p_, q_) in zip("ABC", [(alleles[0], alleles[1]),
                                                 (alleles[2], alleles[3]),
                                                 (alleles[4], alleles[5])])
                if p_ != q_ and "*" in p_ and "*" in q_)
            r["het_loci"] = het or "-"
        else:
            note(sid, "hla", "result file empty")
    else:
        note(sid, "hla", "missing")

    # ---------- mhcflurry ----------
    neo = f"{d}/neoantigens/neoantigens_per_mutation.tsv"
    if os.path.exists(neo):
        n = pd.read_csv(neo, sep="\t")
        r["mhcflurry_mut"] = len(n)
        r["mhcflurry_strong"] = int(n.strong.sum()) if "strong" in n else 0
        if len(n) == 0:
            note(sid, "mhcflurry", "no mutations processed")
    else:
        note(sid, "mhcflurry", "missing")

    # ---------- LOH simulation ----------
    lohbam = f"{LOH}/{sid}/{sid}_tumor_LOH.bam"
    r["loh_bam"] = "ok" if os.path.exists(lohbam + ".bai") else "MISSING"
    if not os.path.exists(lohbam + ".bai"):
        note(sid, "loh_bam", "missing")

    # ---------- LOHHLA ----------
    preds = glob.glob(f"{LOH}/{sid}/out_*/*HLAlossPrediction*.txt")
    real = []
    for p in preds:
        if os.path.getsize(p) < 200:
            continue
        if "homozygous" in p or "No_Suitable" in p:
            continue
        real.append(p)
    r["lohhla_runs"] = len(real)
    if real:
        best = None
        for p in real:
            try:
                t = pd.read_csv(p, sep="\t")
                if len(t) and "PVal" in t:
                    row = t.iloc[0]
                    if best is None or row["PVal"] < best["PVal"]:
                        best = row
            except Exception:
                pass
        if best is not None:
            r["lohhla_pval"] = round(float(best["PVal"]), 4)
            cn1 = best.get("HLA_type1copyNum_withBAFBin")
            cn2 = best.get("HLA_type2copyNum_withBAFBin")
            r["lohhla_CN"] = f"{cn1:.2f}/{cn2:.2f}" if pd.notna(cn1) else None
            r["lohhla_sites"] = best.get("numMisMatchSitesCov")
    else:
        note(sid, "lohhla", "no usable prediction")

    # ---------- pVACseq ----------
    allep = glob.glob(f"{PVAC}/{sid}/pvacseq_out/MHC_Class_I/*.all_epitopes.tsv")
    filt = glob.glob(f"{PVAC}/{sid}/pvacseq_out/MHC_Class_I/*.filtered.tsv")
    if allep:
        r["pvac_all"] = count_lines(allep[0]) - 1
        r["pvac_filtered"] = (count_lines(filt[0]) - 1) if filt else 0
        if r["pvac_all"] == 0:
            note(sid, "pvacseq", "no epitopes")
        elif r["pvac_filtered"] == 0:
            note(sid, "pvacseq", "all epitopes filtered out")
    else:
        r["pvac_all"] = 0

    rows.append(r)

df = pd.DataFrame(rows)
pd.set_option("display.width", 260)
pd.set_option("display.max_columns", 40)

# ---------------- report ----------------
print("=" * 100)
print(" STAGE COMPLETENESS")
print("=" * 100)

stages = [
    ("normal BAM", lambda d: (d.normal_bam == "ok").sum()),
    ("tumour BAM", lambda d: (d.tumor_bam == "ok").sum()),
    ("truth set", lambda d: (d.injected > 0).sum()),
    ("Mutect2 VCF", lambda d: (d.vcf_records > 0).sum()),
    ("HLA type", lambda d: (d.get("n_alleles", pd.Series(0, index=d.index)) >= 4).sum()),
    ("mhcflurry", lambda d: (d.get("mhcflurry_mut", pd.Series(0, index=d.index)) > 0).sum()),
    ("LOH BAM", lambda d: (d.loh_bam == "ok").sum()),
    ("LOHHLA", lambda d: (d.get("lohhla_runs", pd.Series(0, index=d.index)) > 0).sum()),
    ("pVACseq", lambda d: (d.get("pvac_all", pd.Series(0, index=d.index)) > 0).sum()),
]
print(f" {'stage':<16}{'BRCA':>8}{'OV':>8}{'total':>10}")
for name, fn in stages:
    b = fn(df[df.cohort == "brca"])
    o = fn(df[df.cohort == "ov"])
    flag = "" if b + o == len(df) else "   <-- incomplete"
    print(f" {name:<16}{b:>4}/20{o:>6}/20{b+o:>7}/{len(df)}{flag}")

print()
print("=" * 100)
print(" PER SAMPLE")
print("=" * 100)
cols = [c for c in ["sample", "cohort", "expected_mut", "injected", "usable_%",
                    "vcf_PASS", "recall_%", "het_loci", "hla_reads",
                    "mhcflurry_strong", "lohhla_runs", "lohhla_pval",
                    "pvac_all", "pvac_filtered"] if c in df.columns]
print(df[cols].to_string(index=False))

# ---------------- genotype sanity ----------------
print()
print("=" * 100)
print(" HLA GENOTYPE SANITY")
print("=" * 100)
if genotypes:
    uniq = len({g for _, g in genotypes})
    print(f"  distinct genotypes: {uniq} of {len(genotypes)} samples")
    if uniq < len(genotypes) * 0.8:
        print("  WARNING: many identical genotypes - check that OptiType really ran")
    dupes = {}
    for sid, g in genotypes:
        dupes.setdefault(g, []).append(sid)
    for g, ids in dupes.items():
        if len(ids) > 1:
            print(f"  shared by {len(ids)}: {' '.join(ids)}")
            print(f"    {' '.join(a for a in g if a)}")

# ---------------- LOH design ----------------
design = f"{LOH}/loh_design.tsv"
if os.path.exists(design):
    print()
    print("=" * 100)
    print(" LOH DESIGN")
    print("=" * 100)
    dd = pd.read_csv(design, sep="\t")
    print(dd.target_locus.value_counts().to_string())
    if "depth_A_before" in dd.columns:
        print()
        print("  target depth before -> after:")
        for _, row in dd.iterrows():
            loc = row.target_locus.replace("HLA-", "")
            b = row.get(f"depth_{loc}_before")
            a = row.get(f"depth_{loc}_after")
            if pd.notna(b):
                ratio = a / b if b else 0
                flag = "" if 0.25 < ratio < 0.45 else "   <-- unexpected ratio"
                print(f"    {row['sample']:<6} {row.target_locus:<7} "
                      f"{b:>6.1f} -> {a:>6.1f}  ({ratio:.2f}){flag}")

# ---------------- problems ----------------
print()
print("=" * 100)
print(f" PROBLEMS: {len(problems)}")
print("=" * 100)
if problems:
    by_stage = {}
    for sid, stage, msg in problems:
        by_stage.setdefault(stage, []).append((sid, msg))
    for stage in sorted(by_stage):
        print(f"\n  {stage} ({len(by_stage[stage])}):")
        for sid, msg in by_stage[stage]:
            print(f"    {sid:<6} {msg}")
else:
    print("  none")

out = os.path.expanduser("~/immune_escape_project/results/pipeline_audit.tsv")
os.makedirs(os.path.dirname(out), exist_ok=True)
df.to_csv(out, sep="\t", index=False)
print(f"\nwritten to {out}")
