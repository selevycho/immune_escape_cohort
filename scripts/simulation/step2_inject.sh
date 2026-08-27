#!/bin/bash
# =====================================================================
# STEP 2 - build the mutation set and inject it with BAMSurgeon.
#
# For each sample in the batch: pull that patient's mutations out of the
# lifted-over MAF, keep the ones inside the panel, write them in BAMSurgeon
# format, inject, sort, index, then verify the injection actually landed.
#
# SNVs only. Indels need a second BAMSurgeon pass through addindel.py,
# which outputs an unsorted BAM and is a known source of failures; they are
# handled separately in test_indels.sh rather than risking the whole cohort.
#
# The verification step is not optional. BAMSurgeon reports success even
# when the realignment has driven reads at a position down to zero mapping
# quality, which happens systematically across the MHC region. Without
# checking, those positions would later be scored as caller failures rather
# than injection failures.
#
# Memory: every parallel BWA process loads the full hg38 index, about 6 GB.
# Two threads with 16 GB is the safe combination - eight threads reached
# 33 GB and were OOM-killed.
#
# Runtime scales with mutation count, not panel size: roughly 3 seconds per
# injected SNV. A sample with 340 mutations takes around 20 minutes, one
# with 15 takes under two.
#
# Usage:
#   sbatch step2_inject.sh <n_brca> <n_ov> [offset_brca] [offset_ov]
#
#   sbatch step2_inject.sh 5 5          first 5 of each cohort
#   sbatch step2_inject.sh 5 5 5 5      the next 5
# =====================================================================
#SBATCH --job-name=step2_inject
#SBATCH --partition=cpu-single
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/bwfor/work/ws/fr_os136-immune_escape/logs/step2_%j.log
#SBATCH --error=/gpfs/bwfor/work/ws/fr_os136-immune_escape/logs/step2_%j.err
# =====================================================================
set -o pipefail

N_BRCA="${1:-5}"
N_OV="${2:-5}"
OFF_BRCA="${3:-0}"
OFF_OV="${4:-0}"

WS=/gpfs/bwfor/work/ws/fr_os136-immune_escape
MANIFEST=${WS}/simulation/cohort/manifest.tsv
PANEL=${WS}/simulation/panel/panel.bed
REF=${WS}/ref/Homo_sapiens_assembly38.fasta
COHORT_DIR=${WS}/simulation/cohort
SCRIPTS=/home/fr/fr_fr/fr_os136/immune_escape_project/scripts/01_cohort_tcga
CONDA=/home/fr/fr_fr/fr_os136/miniconda3/etc/profile.d/conda.sh

MIN_SNV=5
MIN_VAF=0.05
THREADS="${SLURM_CPUS_PER_TASK:-2}"

source "${CONDA}"

echo "======================================================================"
echo " STEP 2 - mutation set and BAMSurgeon injection (SNV only)"
echo " batch: ${N_BRCA} BRCA (offset ${OFF_BRCA}), ${N_OV} OV (offset ${OFF_OV})"
echo " started $(date '+%F %T') on $(hostname), ${THREADS} threads"
echo "======================================================================"

for f in "${MANIFEST}" "${PANEL}" "${REF}"; do
    [ -s "$f" ] || { echo "ERROR: missing $f"; exit 1; }
done
for s in make_mutations.py verify_injection.py; do
    [ -s "${SCRIPTS}/$s" ] || { echo "ERROR: missing ${SCRIPTS}/$s"; exit 1; }
done

SAMPLES=$(awk -F'\t' -v nb="${N_BRCA}" -v no="${N_OV}" \
                     -v ob="${OFF_BRCA}" -v oo="${OFF_OV}" '
    NR>1 && $2=="brca" { b++; if (b > ob && b <= ob+nb) print $1 }
    NR>1 && $2=="ov"   { o++; if (o > oo && o <= oo+no) print $1 }
' "${MANIFEST}")

echo
echo "samples in this batch: $(echo ${SAMPLES} | tr '\n' ' ')"
echo

OK=0; SKIP=0; FAIL=0; WAIT=0
START_ALL=$(date +%s)

for SID in ${SAMPLES}; do
    LINE=$(awk -F'\t' -v s="${SID}" 'NR>1 && $1==s' "${MANIFEST}")
    COHORT=$(echo "${LINE}"  | cut -f2)
    BARCODE=$(echo "${LINE}" | cut -f3)
    NMUT=$(echo "${LINE}"    | cut -f4)

    OUT=${COHORT_DIR}/${SID}
    NORMAL=${OUT}/${SID}_normal.bam
    TUMOR=${OUT}/${SID}_tumor.bam
    MAF=${WS}/liftover/out/${COHORT}.hg38.maf.tsv

    echo "----------------------------------------------------------------------"
    printf " %-6s %-5s %-16s  expecting ~%s mutations\n" \
           "${SID}" "${COHORT}" "${BARCODE}" "${NMUT}"

    if [ ! -s "${NORMAL}" ] || [ ! -s "${NORMAL}.bai" ]; then
        echo "   normal BAM not ready - step 1 has not reached this sample"
        WAIT=$((WAIT+1)); continue
    fi

    if [ -s "${TUMOR}" ] && [ -s "${TUMOR}.bai" ]; then
        echo "   tumour BAM already present - skipping"
        SKIP=$((SKIP+1)); continue
    fi

    # ---------------- mutation set ----------------
    conda activate cptac_env
    if ! python "${SCRIPTS}/make_mutations.py" \
            "${MAF}" "${PANEL}" "${BARCODE}" "${OUT}" "${MIN_VAF}" \
            > "${OUT}/make_mutations.log" 2>&1; then
        echo "   FAILED building the mutation set:"
        tail -4 "${OUT}/make_mutations.log" | sed 's/^/     /'
        printf "%s\tmake_mutations\n" "${SID}" >> "${COHORT_DIR}/step2_failures.log"
        FAIL=$((FAIL+1)); continue
    fi

    N_SNV=$(wc -l < "${OUT}/snv_mutations.txt")
    N_IND=$(wc -l < "${OUT}/indel_mutations.txt" 2>/dev/null || echo 0)
    VAF_MED=$(awk '{print $4}' "${OUT}/snv_mutations.txt" | sort -n \
              | awk '{a[NR]=$1} END {if(NR>0) printf "%.3f", a[int(NR/2)+1]; else print "NA"}')
    echo "   ${N_SNV} SNVs (median VAF ${VAF_MED}), ${N_IND} indels set aside"

    if [ "${N_SNV}" -lt "${MIN_SNV}" ]; then
        echo "   below the ${MIN_SNV}-mutation floor - skipping"
        printf "%s\ttoo_few_%s\n" "${SID}" "${N_SNV}" >> "${COHORT_DIR}/step2_failures.log"
        FAIL=$((FAIL+1)); continue
    fi

    # ---------------- injection ----------------
    conda activate bamsurgeon_env
    PICARD=$(find "${CONDA_PREFIX}" -name "picard.jar" 2>/dev/null | head -1)
    if [ -z "${PICARD}" ]; then
        echo "   picard.jar not found in bamsurgeon_env"
        FAIL=$((FAIL+1)); continue
    fi

    T0=$(date +%s)
    echo "   $(date '+%T') running addsnv.py ..."

    # addsnv.py writes its scratch directory into the working directory
    cd "${OUT}"
    rm -rf addsnv.tmp

    if ! addsnv.py \
            -v "${OUT}/snv_mutations.txt" \
            -f "${NORMAL}" \
            -r "${REF}" \
            -o "${OUT}/tumor_raw.bam" \
            -p "${THREADS}" \
            --aligner mem \
            --picardjar "${PICARD}" \
            --force --insane \
            > "${OUT}/bamsurgeon.log" 2>&1; then
        echo "   FAILED in addsnv.py:"
        tail -6 "${OUT}/bamsurgeon.log" | sed 's/^/     /'
        printf "%s\taddsnv\n" "${SID}" >> "${COHORT_DIR}/step2_failures.log"
        rm -rf "${OUT}/addsnv.tmp" "${OUT}/tumor_raw.bam"*
        FAIL=$((FAIL+1)); continue
    fi

    conda activate bio_work
    if ! samtools sort -@ "${THREADS}" -m 1G -o "${TUMOR}" "${OUT}/tumor_raw.bam" \
            2>> "${OUT}/bamsurgeon.log"; then
        echo "   FAILED sorting the injected BAM"
        printf "%s\tsort\n" "${SID}" >> "${COHORT_DIR}/step2_failures.log"
        rm -f "${TUMOR}"
        FAIL=$((FAIL+1)); continue
    fi
    samtools index -@ "${THREADS}" "${TUMOR}"

    # BAMSurgeon leaves a per-mutation log directory behind; at 40 samples
    # that is thousands of small files for no further use
    rm -rf "${OUT}/tumor_raw.bam"* "${OUT}/addsnv.tmp" \
           "${OUT}"/addsnv_logs_*

    T1=$(date +%s)
    printf "   injected in %dm%02ds  |  %s\n" \
           $(( (T1-T0)/60 )) $(( (T1-T0)%60 )) \
           "$(ls -lh ${TUMOR} | awk '{print $5}')"

    # ---------------- verification ----------------
    conda activate cptac_env
    if python "${SCRIPTS}/verify_injection.py" \
            "${OUT}/truth_set.tsv" "${NORMAL}" "${TUMOR}" \
            "${OUT}/verification" > "${OUT}/verification.log" 2>&1; then
        grep -E "usable as ground truth" "${OUT}/verification.log" | sed 's/^/   /'
        grep -E "Pearson r|mean error" "${OUT}/verification.log" | head -2 | sed 's/^/   /'
    else
        echo "   WARNING: verification did not complete"
        tail -3 "${OUT}/verification.log" | sed 's/^/     /'
    fi

    OK=$((OK+1))
done

END_ALL=$(date +%s)
echo "======================================================================"
echo " STEP 2 finished $(date '+%F %T')"
printf " elapsed %dh%02dm   injected %d   skipped %d   waiting on step 1: %d   failed %d\n" \
       $(( (END_ALL-START_ALL)/3600 )) $(( ((END_ALL-START_ALL)%3600)/60 )) \
       "${OK}" "${SKIP}" "${WAIT}" "${FAIL}"
echo "======================================================================"
echo
echo "cohort progress:"
NB=$(ls ${COHORT_DIR}/*/*_normal.bam 2>/dev/null | wc -l)
NT=$(ls ${COHORT_DIR}/*/*_tumor.bam 2>/dev/null | wc -l)
NM=$(awk 'NR>1' "${MANIFEST}" | wc -l)
echo "  normal BAMs: ${NB} / ${NM}"
echo "  tumour BAMs: ${NT} / ${NM}"
