#!/bin/bash
# =====================================================================
# STEP 14b - call variants at each depth and score against the truth.
#
# Same Mutect2 invocation as step 3, same panel, same reference. Only the
# input depth changes, so recall differences between levels are coverage
# and nothing else.
#
# Read groups are rewritten before calling. Both BAMs descend from the
# same 1000 Genomes individual and carry an identical sample tag, which
# Mutect2 refuses to work with; addreplacerg gives them distinct ones
# without touching the reads.
#
# Only mutations verified present in the full-depth tumour are scored.
# A mutation BAMSurgeon could not place is not a caller failure at any
# depth, and counting it as one would drag every level down equally.
#
# Usage:
#   sbatch step14b_call_coverage.sh <level> <n_brca> <n_ov> [off_b] [off_o]
#
#   sbatch step14b_call_coverage.sh 10 1 0
# =====================================================================
#SBATCH --job-name=step14b_call
#SBATCH --partition=cpu-single
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=16:00:00
#SBATCH --output=/gpfs/bwfor/work/ws/fr_os136-immune_escape/simulation/step14_coverage_sweep/logs/step14b_%j.log
#SBATCH --error=/gpfs/bwfor/work/ws/fr_os136-immune_escape/simulation/step14_coverage_sweep/logs/step14b_%j.err
# =====================================================================

LEVEL="${1:?give a level, e.g. 10}"
N_BRCA="${2:-20}"
N_OV="${3:-20}"
OFF_BRCA="${4:-0}"
OFF_OV="${5:-0}"

WS=/gpfs/bwfor/work/ws/fr_os136-immune_escape
COHORT=${WS}/simulation/cohort
SWEEP=${WS}/simulation/step14_coverage_sweep
LEVEL_DIR=${SWEEP}/cov${LEVEL}x
MANIFEST=${COHORT}/manifest.tsv
REF=${WS}/ref/Homo_sapiens_assembly38.fasta
PANEL=${WS}/simulation/panel/panel.bed
CONDA=/home/fr/fr_fr/fr_os136/miniconda3/etc/profile.d/conda.sh
THREADS="${SLURM_CPUS_PER_TASK:-4}"

source "${CONDA}"
conda activate gatk_env 2>/dev/null || conda activate bio_work
command -v gatk >/dev/null || { echo "ERROR: gatk not on PATH"; exit 1; }

echo "======================================================================"
echo " STEP 14b - Mutect2 at ${LEVEL}x"
echo " batch: ${N_BRCA} BRCA (off ${OFF_BRCA}), ${N_OV} OV (off ${OFF_OV})"
echo " started $(date '+%F %T') on $(hostname)"
echo "======================================================================"

[ -d "${LEVEL_DIR}" ] || { echo "ERROR: ${LEVEL_DIR} missing - run 14a"; exit 1; }

SAMPLES=$(awk -F'\t' -v nb="${N_BRCA}" -v no="${N_OV}" \
                     -v ob="${OFF_BRCA}" -v oo="${OFF_OV}" '
    NR>1 && $2=="brca" { b++; if (b > ob && b <= ob+nb) print $1 }
    NR>1 && $2=="ov"   { o++; if (o > oo && o <= oo+no) print $1 }
' "${MANIFEST}")

OK=0; SKIP=0; FAIL=0; WAIT=0

for SID in ${SAMPLES}; do
    DIR=${LEVEL_DIR}/${SID}
    T_BAM=${DIR}/${SID}_tumor.bam
    N_BAM=${DIR}/${SID}_normal.bam
    TRUTH=${DIR}/truth_set.tsv
    [ -s "${TRUTH}" ] || TRUTH=${COHORT}/${SID}/truth_set.tsv
    M=${DIR}/mutect2

    echo "----------------------------------------------------------------------"
    echo " ${SID}"

    if [ ! -s "${T_BAM}" ] || [ ! -s "${N_BAM}" ]; then
        echo "   not built at this level"
        WAIT=$((WAIT+1)); continue
    fi

    if [ -s "${M}/${SID}.filtered.vcf.gz.tbi" ]; then
        echo "   already called"
        SKIP=$((SKIP+1)); continue
    fi

    mkdir -p "${M}"
    T0=$(date +%s)

    # distinct read groups, or Mutect2 cannot tell the two apart
    for KIND in normal tumor; do
        RG=${DIR}/${SID}_${KIND}.rg.bam
        if [ ! -s "${RG}.bai" ]; then
            samtools addreplacerg -@ "${THREADS}" \
                -r "ID:${SID}_${KIND}" -r "SM:${SID}_${KIND^^}" \
                -r "PL:ILLUMINA" -r "LB:${SID}" -r "PU:unit1" \
                -o "${RG}" "${DIR}/${SID}_${KIND}.bam" 2>/dev/null
            samtools index -@ "${THREADS}" "${RG}"
        fi
    done

    echo "   $(date '+%T') Mutect2 ..."
    gatk --java-options "-Xmx20g" Mutect2 \
        -R "${REF}" \
        -I "${DIR}/${SID}_tumor.rg.bam" \
        -I "${DIR}/${SID}_normal.rg.bam" \
        -normal "${SID}_NORMAL" \
        -L "${PANEL}" \
        -O "${M}/${SID}.unfiltered.vcf.gz" \
        > "${M}/mutect2.log" 2>&1
    RC=$?

    if [ "${RC}" -ne 0 ] || [ ! -s "${M}/${SID}.unfiltered.vcf.gz" ]; then
        echo "   Mutect2 FAILED (exit ${RC}):"
        tail -5 "${M}/mutect2.log" | sed 's/^/     /'
        printf "%s\t%sx\tmutect2_rc%s\n" "${SID}" "${LEVEL}" "${RC}" \
            >> "${SWEEP}/failures.log"
        FAIL=$((FAIL+1)); continue
    fi

    gatk FilterMutectCalls \
        -R "${REF}" \
        -V "${M}/${SID}.unfiltered.vcf.gz" \
        -O "${M}/${SID}.filtered.vcf.gz" \
        >> "${M}/mutect2.log" 2>&1

    T1=$(date +%s)
    N_ALL=$(zcat "${M}/${SID}.unfiltered.vcf.gz" | grep -vc '^#')
    N_PASS=$(zcat "${M}/${SID}.filtered.vcf.gz" \
             | awk -F'\t' '!/^#/ && $7=="PASS"' | wc -l)
    printf "   %dm%02ds   %d records, %d PASS\n" \
           $(( (T1-T0)/60 )) $(( (T1-T0)%60 )) "${N_ALL}" "${N_PASS}"

    # ---- score against the truth ----
    conda activate cptac_env 2>/dev/null
    mkdir -p "${DIR}/comparison"
    python3 - "${TRUTH}" "${M}/${SID}.filtered.vcf.gz" \
              "${DIR}/comparison/truth_vs_calls.tsv" "${SID}" << 'PYEOF'
import sys, gzip
import pandas as pd

truth_p, vcf_p, out_p, sid = sys.argv[1:5]

t = pd.read_csv(truth_p, sep="\t")
if "Variant_Type" in t.columns:
    t = t[t.Variant_Type == "SNP"]      # indels go through step 8

calls = {}
with gzip.open(vcf_p, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        key = (f[0], int(f[1]))
        calls[key] = {"filter": f[6], "pass": f[6] == "PASS"}

rows = []
for _, r in t.iterrows():
    key = (r.Chromosome_hg38, int(r.Start_Position_hg38))
    c = calls.get(key)
    rows.append({
        "Hugo_Symbol": r.get("Hugo_Symbol"),
        "Chromosome_hg38": key[0], "Start_Position_hg38": key[1],
        "Reference_Allele": r.get("Reference_Allele"),
        "Tumor_Seq_Allele2": r.get("Tumor_Seq_Allele2"),
        "VAF": r.get("VAF"),
        "detected_any": c is not None,
        "detected_pass": bool(c and c["pass"]),
        "filters": c["filter"] if c else "",
    })

d = pd.DataFrame(rows)
d.to_csv(out_p, sep="\t", index=False)
if len(d):
    print("   recall %.1f%% (%d of %d), seen %.1f%%"
          % (100*d.detected_pass.mean(), int(d.detected_pass.sum()), len(d),
             100*d.detected_any.mean()))
PYEOF
    conda activate gatk_env 2>/dev/null || conda activate bio_work

    rm -f "${DIR}/${SID}_normal.rg.bam"* "${DIR}/${SID}_tumor.rg.bam"*
    OK=$((OK+1))
done

echo "======================================================================"
echo " STEP 14b finished $(date '+%F %T')  -  ${LEVEL}x"
echo " called ${OK}   skipped ${SKIP}   waiting ${WAIT}   failed ${FAIL}"
echo "======================================================================"
