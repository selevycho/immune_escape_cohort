#!/bin/bash
# =====================================================================
# STEP 15b - call variants at each allele fraction.
#
# Same Mutect2 invocation as step 3. Only the fraction differs between
# levels, so recall against fraction is measured on one set of positions
# rather than compared across mutations that differ in other ways.
#
# Scoring uses each level's own truth set, which lists only the mutations
# verified present in that level's BAM. A mutation BAMSurgeon could not
# place at 0.05 is not a caller failure there.
#
# Usage:
#   sbatch step15b_call_purity.sh <level> <n_brca> <n_ov> [off_b] [off_o]
#
#   sbatch step15b_call_purity.sh vaf05 5 5
# =====================================================================
#SBATCH --job-name=step15b_call
#SBATCH --partition=cpu-single
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=16:00:00
#SBATCH --output=/gpfs/bwfor/work/ws/fr_os136-immune_escape/simulation/step15_purity_sweep/logs/step15b_%j.log
#SBATCH --error=/gpfs/bwfor/work/ws/fr_os136-immune_escape/simulation/step15_purity_sweep/logs/step15b_%j.err
# =====================================================================

LEVEL="${1:?give a level, e.g. vaf05}"
N_BRCA="${2:-5}"
N_OV="${3:-5}"
OFF_BRCA="${4:-0}"
OFF_OV="${5:-0}"

WS=/gpfs/bwfor/work/ws/fr_os136-immune_escape
COHORT=${WS}/simulation/cohort
SWEEP=${WS}/simulation/step15_purity_sweep
LEVEL_DIR=${SWEEP}/${LEVEL}
MANIFEST=${COHORT}/manifest.tsv
REF=${WS}/ref/Homo_sapiens_assembly38.fasta
PANEL=${WS}/simulation/panel/panel.bed
CONDA=/home/fr/fr_fr/fr_os136/miniconda3/etc/profile.d/conda.sh
THREADS="${SLURM_CPUS_PER_TASK:-4}"

source "${CONDA}"
conda activate gatk_env 2>/dev/null || conda activate bio_work
command -v gatk >/dev/null || { echo "ERROR: gatk not on PATH"; exit 1; }

echo "======================================================================"
echo " STEP 15b - Mutect2 at ${LEVEL}"
echo " started $(date '+%F %T') on $(hostname)"
echo "======================================================================"

[ -d "${LEVEL_DIR}" ] || { echo "ERROR: ${LEVEL_DIR} missing - run 15a"; exit 1; }

SAMPLES=$(awk -F'\t' -v nb="${N_BRCA}" -v no="${N_OV}" \
                     -v ob="${OFF_BRCA}" -v oo="${OFF_OV}" '
    NR>1 && $2=="brca" { b++; if (b > ob && b <= ob+nb) print $1 }
    NR>1 && $2=="ov"   { o++; if (o > oo && o <= oo+no) print $1 }
' "${MANIFEST}")

OK=0; SKIP=0; FAIL=0; WAIT=0

for SID in ${SAMPLES}; do
    DIR=${LEVEL_DIR}/${SID}
    M=${DIR}/mutect2

    echo "----------------------------------------------------------------------"
    echo " ${SID}"

    if [ ! -s "${DIR}/${SID}_tumor.bam" ]; then
        echo "   not built at this level"
        WAIT=$((WAIT+1)); continue
    fi
    if [ -s "${M}/${SID}.filtered.vcf.gz.tbi" ]; then
        echo "   already called"
        SKIP=$((SKIP+1)); continue
    fi

    mkdir -p "${M}"
    T0=$(date +%s)

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
        echo "   Mutect2 FAILED (exit ${RC})"
        tail -4 "${M}/mutect2.log" | sed 's/^/     /'
        printf "%s\t%s\tmutect2\n" "${SID}" "${LEVEL}" >> "${SWEEP}/failures.log"
        FAIL=$((FAIL+1)); continue
    fi

    gatk FilterMutectCalls -R "${REF}" \
        -V "${M}/${SID}.unfiltered.vcf.gz" \
        -O "${M}/${SID}.filtered.vcf.gz" \
        >> "${M}/mutect2.log" 2>&1

    T1=$(date +%s)
    N_PASS=$(zcat "${M}/${SID}.filtered.vcf.gz" \
             | awk -F'\t' '!/^#/ && $7=="PASS"' | wc -l)
    printf "   %dm%02ds   %d PASS\n" \
           $(( (T1-T0)/60 )) $(( (T1-T0)%60 )) "${N_PASS}"

    conda activate cptac_env 2>/dev/null
    mkdir -p "${DIR}/comparison"
    python3 - "${DIR}/truth_set.tsv" "${M}/${SID}.filtered.vcf.gz" \
              "${DIR}/comparison/truth_vs_calls.tsv" << 'PYEOF'
import sys, gzip
import pandas as pd

truth_p, vcf_p, out_p = sys.argv[1:4]
t = pd.read_csv(truth_p, sep="\t")

calls = {}
with gzip.open(vcf_p, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        calls[(f[0], int(f[1]))] = {"filter": f[6], "pass": f[6] == "PASS"}

rows = []
for _, r in t.iterrows():
    key = (r.Chromosome_hg38, int(r.Start_Position_hg38))
    c = calls.get(key)
    rows.append({
        "Hugo_Symbol": r.get("Hugo_Symbol"),
        "Chromosome_hg38": key[0], "Start_Position_hg38": key[1],
        "VAF": r.get("VAF"),
        "detected_any": c is not None,
        "detected_pass": bool(c and c["pass"]),
        "filters": c["filter"] if c else "",
    })

d = pd.DataFrame(rows)
d.to_csv(out_p, sep="\t", index=False)
if len(d):
    print("   recall %.1f%% (%d of %d), seen %.1f%%"
          % (100*d.detected_pass.mean(), int(d.detected_pass.sum()),
             len(d), 100*d.detected_any.mean()))
PYEOF
    conda activate gatk_env 2>/dev/null || conda activate bio_work

    rm -f "${DIR}/${SID}_normal.rg.bam"* "${DIR}/${SID}_tumor.rg.bam"*
    OK=$((OK+1))
done

echo "======================================================================"
echo " STEP 15b finished $(date '+%F %T')  -  ${LEVEL}"
echo " called ${OK}   skipped ${SKIP}   waiting ${WAIT}   failed ${FAIL}"
echo "======================================================================"
