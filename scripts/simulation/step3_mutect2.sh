#!/bin/bash
# =====================================================================
# STEP 3 - somatic variant calling and comparison against the truth set.
#
# For each sample in the batch: rewrite the read groups so Mutect2 can tell
# tumour from normal, call variants restricted to the panel, filter, then
# compare what came back against the mutations that were actually injected.
#
# The read group rewrite is not cosmetic. Both BAMs descend from the same
# 1000 Genomes individual and therefore carry an identical SM tag; Mutect2
# refuses to run on a pair it cannot distinguish. samtools addreplacerg
# does this as a stream operation with no realignment.
#
# The comparison uses truth_set.tsv rather than the injected variant file,
# because verification in step 2 already established which positions really
# received their mutation. Positions where the injection failed must not be
# counted against the caller.
#
# Usage:
#   sbatch step3_mutect2.sh <n_brca> <n_ov> [offset_brca] [offset_ov]
# =====================================================================
#SBATCH --job-name=step3_mutect2
#SBATCH --partition=cpu-single
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=08:00:00
#SBATCH --output=/gpfs/bwfor/work/ws/fr_os136-immune_escape/logs/step3_%j.log
#SBATCH --error=/gpfs/bwfor/work/ws/fr_os136-immune_escape/logs/step3_%j.err
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

THREADS="${SLURM_CPUS_PER_TASK:-2}"

source "${CONDA}"

echo "======================================================================"
echo " STEP 3 - Mutect2 and truth comparison"
echo " batch: ${N_BRCA} BRCA (offset ${OFF_BRCA}), ${N_OV} OV (offset ${OFF_OV})"
echo " started $(date '+%F %T') on $(hostname)"
echo "======================================================================"

[ -s "${SCRIPTS}/compare_mutect2_truth.py" ] \
    || { echo "ERROR: compare_mutect2_truth.py missing"; exit 1; }

SAMPLES=$(awk -F'\t' -v nb="${N_BRCA}" -v no="${N_OV}" \
                     -v ob="${OFF_BRCA}" -v oo="${OFF_OV}" '
    NR>1 && $2=="brca" { b++; if (b > ob && b <= ob+nb) print $1 }
    NR>1 && $2=="ov"   { o++; if (o > oo && o <= oo+no) print $1 }
' "${MANIFEST}")

echo
echo "samples: $(echo ${SAMPLES} | tr '\n' ' ')"
echo

OK=0; SKIP=0; FAIL=0; WAIT=0
START_ALL=$(date +%s)

for SID in ${SAMPLES}; do
    LINE=$(awk -F'\t' -v s="${SID}" 'NR>1 && $1==s' "${MANIFEST}")
    COHORT=$(echo "${LINE}" | cut -f2)

    OUT=${COHORT_DIR}/${SID}
    NORMAL=${OUT}/${SID}_normal.bam
    TUMOR=${OUT}/${SID}_tumor.bam
    VCF=${OUT}/mutect2/${SID}.filtered.vcf.gz

    echo "----------------------------------------------------------------------"
    printf " %-6s %-5s\n" "${SID}" "${COHORT}"

    if [ ! -s "${TUMOR}" ]; then
        echo "   tumour BAM not ready - step 2 has not reached this sample"
        WAIT=$((WAIT+1)); continue
    fi
    if [ ! -s "${OUT}/truth_set.tsv" ]; then
        echo "   truth set missing - step 2 did not complete here"
        WAIT=$((WAIT+1)); continue
    fi

    if [ -s "${VCF}" ]; then
        echo "   VCF already present - skipping the call"
        SKIP=$((SKIP+1))
    else
        conda activate bio_work
        mkdir -p "${OUT}/mutect2"
        T0=$(date +%s)

        # both BAMs come from the same individual, so SM must be rewritten
        samtools addreplacerg -@ "${THREADS}" \
            -r "ID:TUMOR\tSM:${SID}_TUMOR\tPL:ILLUMINA\tLB:panel\tPU:sim" \
            -o "${OUT}/t.rg.bam" "${TUMOR}" 2>/dev/null \
            || { echo "   FAILED rewriting tumour read groups"; FAIL=$((FAIL+1)); continue; }
        samtools index -@ "${THREADS}" "${OUT}/t.rg.bam"

        samtools addreplacerg -@ "${THREADS}" \
            -r "ID:NORMAL\tSM:${SID}_NORMAL\tPL:ILLUMINA\tLB:panel\tPU:sim" \
            -o "${OUT}/n.rg.bam" "${NORMAL}" 2>/dev/null \
            || { echo "   FAILED rewriting normal read groups"; FAIL=$((FAIL+1)); continue; }
        samtools index -@ "${THREADS}" "${OUT}/n.rg.bam"

        echo "   $(date '+%T') running Mutect2 ..."
        if ! gatk --java-options "-Xmx8g" Mutect2 \
                -R "${REF}" \
                -I "${OUT}/t.rg.bam" -I "${OUT}/n.rg.bam" \
                -normal "${SID}_NORMAL" \
                -L "${PANEL}" \
                --native-pair-hmm-threads "${THREADS}" \
                -O "${OUT}/mutect2/${SID}.unfiltered.vcf.gz" \
                > "${OUT}/mutect2.log" 2>&1; then
            echo "   FAILED in Mutect2:"
            tail -6 "${OUT}/mutect2.log" | sed 's/^/     /'
            printf "%s\tmutect2\n" "${SID}" >> "${COHORT_DIR}/step3_failures.log"
            rm -f "${OUT}"/[tn].rg.bam*
            FAIL=$((FAIL+1)); continue
        fi

        if ! gatk --java-options "-Xmx6g" FilterMutectCalls \
                -R "${REF}" \
                -V "${OUT}/mutect2/${SID}.unfiltered.vcf.gz" \
                -O "${VCF}" >> "${OUT}/mutect2.log" 2>&1; then
            echo "   FAILED in FilterMutectCalls"
            printf "%s\tfilter\n" "${SID}" >> "${COHORT_DIR}/step3_failures.log"
            rm -f "${OUT}"/[tn].rg.bam*
            FAIL=$((FAIL+1)); continue
        fi

        # the read-group copies are only needed during the call
        rm -f "${OUT}"/[tn].rg.bam*

        T1=$(date +%s)
        N_ALL=$(zcat "${VCF}" | grep -vc "^#")
        N_PASS=$(zcat "${VCF}" | grep -v "^#" | awk -F'\t' '$7=="PASS"' | wc -l)
        printf "   called in %dm%02ds  |  %s records, %s PASS\n" \
               $(( (T1-T0)/60 )) $(( (T1-T0)%60 )) "${N_ALL}" "${N_PASS}"
        OK=$((OK+1))
    fi

    # ---------------- compare against the truth ----------------
    conda activate cptac_env
    if python "${SCRIPTS}/compare_mutect2_truth.py" \
            "${OUT}/truth_set.tsv" "${VCF}" "${OUT}/comparison" \
            > "${OUT}/comparison.log" 2>&1; then
        grep -E "recovered, PASS|false positives|Pearson r" \
             "${OUT}/comparison.log" | sed 's/^/   /'
    else
        echo "   WARNING: comparison failed"
        tail -3 "${OUT}/comparison.log" | sed 's/^/     /'
    fi
done

END_ALL=$(date +%s)
echo "======================================================================"
echo " STEP 3 finished $(date '+%F %T')"
printf " elapsed %dh%02dm   called %d   skipped %d   waiting %d   failed %d\n" \
       $(( (END_ALL-START_ALL)/3600 )) $(( ((END_ALL-START_ALL)%3600)/60 )) \
       "${OK}" "${SKIP}" "${WAIT}" "${FAIL}"
echo "======================================================================"
echo
echo "cohort progress:"
NM=$(awk 'NR>1' "${MANIFEST}" | wc -l)
echo "  normal BAMs : $(ls ${COHORT_DIR}/*/*_normal.bam 2>/dev/null | wc -l) / ${NM}"
echo "  tumour BAMs : $(ls ${COHORT_DIR}/*/*_tumor.bam 2>/dev/null | wc -l) / ${NM}"
echo "  VCFs        : $(ls ${COHORT_DIR}/*/mutect2/*.filtered.vcf.gz 2>/dev/null | wc -l) / ${NM}"
