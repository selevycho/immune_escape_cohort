#!/bin/bash
#SBATCH --job-name=step8b_indel_call
#SBATCH --partition=cpu-single
#SBATCH --time=06:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --output=%x_%j.log
#
# Step 8b - call variants on the indel BAMs.
#
# Step 8 placed indels into a copy of the tumour BAM and nothing was ever
# called against it, so the pipeline has been shown to place indels and
# not to find them. This closes that gap.
#
# The BAM already carries the substitutions from step 2, so a single run
# recovers both kinds and the comparison separates them by variant type
# rather than by file. That is deliberate: it is also how a real sample
# behaves, with both kinds present at once.
#
# Read groups are rewritten exactly as step 3 does. Tumour and normal are
# sliced from the same 1000 Genomes individual and carry an identical SM
# tag; Mutect2 refuses a pair it cannot tell apart.
#
# Usage:
#   sbatch step8b_call_indels.sh              all samples with indels
#   sbatch step8b_call_indels.sh B001 B002    named samples
set -uo pipefail

WS=$(ws_find immune_escape)
COHORT="${WS}/simulation/cohort"
INDELS="${WS}/simulation/indels"
PANEL="${WS}/simulation/panel/panel.bed"
REF="${WS}/ref/Homo_sapiens_assembly38.fasta"
TRUTH="${INDELS}/indel_truth_all.tsv"
THREADS=4

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate bio_work

if [ $# -gt 0 ]; then
    SAMPLES="$*"
else
    SAMPLES=$(awk -F'\t' 'NR>1 {print $1}' "${TRUTH}" | sort -u)
fi

echo "======================================================================"
echo " STEP 8b - Mutect2 on the indel BAMs"
echo "======================================================================"
echo "samples: $(echo ${SAMPLES} | tr '\n' ' ')"
echo "started: $(date '+%H:%M:%S')"
echo "----------------------------------------------------------------------"

DONE=0; SKIP=0; FAIL=0

for SID in ${SAMPLES}; do
    TBAM="${INDELS}/${SID}/${SID}_tumor_snv_indel.bam"
    NBAM="${COHORT}/${SID}/${SID}_normal.bam"
    OUT="${INDELS}/${SID}/calls"

    echo ""
    echo " ${SID}"

    if [ ! -s "${TBAM}" ] || [ ! -s "${NBAM}" ]; then
        echo "   BAM missing - skipped"
        SKIP=$((SKIP+1)); continue
    fi
    if [ -s "${OUT}/${SID}.filtered.vcf.gz" ]; then
        echo "   already called - skipped"
        SKIP=$((SKIP+1)); continue
    fi

    mkdir -p "${OUT}"
    T0=$(date +%s)

    # both BAMs come from the same individual, so SM must be rewritten
    echo -n "   read groups ... "
    samtools addreplacerg -@ "${THREADS}" \
        -r "ID:TUMOR\tSM:${SID}_TUMOR\tPL:ILLUMINA\tLB:panel\tPU:sim" \
        -o "${OUT}/tumour.rg.bam" "${TBAM}" 2>"${OUT}/rg.err" \
        && samtools index -@ "${THREADS}" "${OUT}/tumour.rg.bam" \
        || { echo "failed"; FAIL=$((FAIL+1)); continue; }

    samtools addreplacerg -@ "${THREADS}" \
        -r "ID:NORMAL\tSM:${SID}_NORMAL\tPL:ILLUMINA\tLB:panel\tPU:sim" \
        -o "${OUT}/normal.rg.bam" "${NBAM}" 2>>"${OUT}/rg.err" \
        && samtools index -@ "${THREADS}" "${OUT}/normal.rg.bam" \
        || { echo "failed"; FAIL=$((FAIL+1)); continue; }
    echo "ok"

    echo -n "   Mutect2 ... "
    if gatk --java-options "-Xmx12g" Mutect2 \
            -R "${REF}" \
            -I "${OUT}/tumour.rg.bam" \
            -I "${OUT}/normal.rg.bam" \
            -normal "${SID}_NORMAL" \
            -L "${PANEL}" \
            -O "${OUT}/${SID}.raw.vcf.gz" \
            > "${OUT}/mutect2.log" 2>&1; then
        echo "ok"
    else
        echo "failed - see ${OUT}/mutect2.log"
        FAIL=$((FAIL+1)); continue
    fi

    echo -n "   FilterMutectCalls ... "
    if gatk --java-options "-Xmx8g" FilterMutectCalls \
            -R "${REF}" \
            -V "${OUT}/${SID}.raw.vcf.gz" \
            -O "${OUT}/${SID}.filtered.vcf.gz" \
            > "${OUT}/filter.log" 2>&1; then
        echo "ok"
    else
        echo "failed - see ${OUT}/filter.log"
        FAIL=$((FAIL+1)); continue
    fi

    # bcftools in this environment is missing libcrypto and exits before
    # printing anything, so the VCF is read with zcat instead
    VCF="${OUT}/${SID}.filtered.vcf.gz"
    N_ALL=$(zcat "${VCF}" | grep -vc "^#")
    N_PASS=$(zcat "${VCF}" | grep -v "^#" | awk -F'\t' '$7=="PASS"' | wc -l)
    # an indel is a record whose ref and alt differ in length
    N_IND=$(zcat "${VCF}" | grep -v "^#" \
            | awk -F'\t' '$7=="PASS" && length($4)!=length($5)' | wc -l)

    rm -f "${OUT}/tumour.rg.bam" "${OUT}/tumour.rg.bam.bai" \
          "${OUT}/normal.rg.bam" "${OUT}/normal.rg.bam.bai"

    T1=$(date +%s)
    printf "   %d calls, %d PASS, %d of those indels   %dm%02ds\n" \
        "${N_ALL}" "${N_PASS}" "${N_IND}" $(( (T1-T0)/60 )) $(( (T1-T0)%60 ))
    DONE=$((DONE+1))
done

echo ""
echo "======================================================================"
echo " called ${DONE}   skipped ${SKIP}   failed ${FAIL}"
echo " finished: $(date '+%H:%M:%S')"
echo "======================================================================"
