#!/bin/bash
#SBATCH --job-name=step17_negctrl
#SBATCH --partition=cpu-single
#SBATCH --time=08:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --output=%x_%j.log
#
# Step 17 - the negative control.
#
# Everything measured so far asks what the caller recovers from files that
# carry mutations we placed. None of it asks what the caller reports from
# files that carry none. A pipeline that produced calls regardless of
# input would score identically on recall and be worthless, and nothing
# done so far would have caught that.
#
# So Mutect2 is run on pairs where both sides are the same untouched normal
# slice. There is no tumour: the same file is passed twice, under two read
# group names, so the caller has a somatic comparison to make and nothing
# somatic to find. Every PASS call is a false positive by construction.
#
# This is the control the specificity figure has been resting on
# implicitly. Measured, it either confirms the figure or invalidates it.
#
# Usage:
#   sbatch step17_negative_control.sh              all forty
#   sbatch step17_negative_control.sh B001 B002    named samples
set -uo pipefail

WS=$(ws_find immune_escape)
COHORT="${WS}/simulation/cohort"
OUTDIR="${WS}/simulation/step17_negative_control"
PANEL="${WS}/simulation/panel/panel.bed"
REF="${WS}/ref/Homo_sapiens_assembly38.fasta"
THREADS=4

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate bio_work

if [ $# -gt 0 ]; then
    SAMPLES="$*"
else
    SAMPLES=$(awk -F'\t' 'NR>1 {print $1}' "${COHORT}/manifest.tsv")
fi

mkdir -p "${OUTDIR}"

echo "======================================================================"
echo " STEP 17 - negative control: normal against itself"
echo "======================================================================"
echo "samples: $(echo ${SAMPLES} | tr '\n' ' ')"
echo "started: $(date '+%H:%M:%S')"
echo "----------------------------------------------------------------------"

DONE=0; SKIP=0; FAIL=0; TOTAL_PASS=0

for SID in ${SAMPLES}; do
    NBAM="${COHORT}/${SID}/${SID}_normal.bam"
    OUT="${OUTDIR}/${SID}"

    echo ""
    echo " ${SID}"

    if [ ! -s "${NBAM}" ]; then
        echo "   normal BAM missing - skipped"
        SKIP=$((SKIP+1)); continue
    fi
    if [ -s "${OUT}/${SID}.filtered.vcf.gz" ]; then
        N=$(zcat "${OUT}/${SID}.filtered.vcf.gz" | grep -v "^#" \
            | awk -F'\t' '$7=="PASS"' | wc -l)
        echo "   already called - ${N} PASS"
        TOTAL_PASS=$((TOTAL_PASS+N))
        SKIP=$((SKIP+1)); continue
    fi

    mkdir -p "${OUT}"
    T0=$(date +%s)

    # the same file twice, under two sample names, so that Mutect2 has a
    # pair to compare and nothing somatic separates them
    echo -n "   read groups ... "
    samtools addreplacerg -@ "${THREADS}" \
        -r "ID:TUMOR\tSM:${SID}_TUMOR\tPL:ILLUMINA\tLB:panel\tPU:sim" \
        -o "${OUT}/as_tumour.bam" "${NBAM}" 2>"${OUT}/rg.err" \
        && samtools index -@ "${THREADS}" "${OUT}/as_tumour.bam" \
        || { echo "failed"; FAIL=$((FAIL+1)); continue; }

    samtools addreplacerg -@ "${THREADS}" \
        -r "ID:NORMAL\tSM:${SID}_NORMAL\tPL:ILLUMINA\tLB:panel\tPU:sim" \
        -o "${OUT}/as_normal.bam" "${NBAM}" 2>>"${OUT}/rg.err" \
        && samtools index -@ "${THREADS}" "${OUT}/as_normal.bam" \
        || { echo "failed"; FAIL=$((FAIL+1)); continue; }
    echo "ok"

    echo -n "   Mutect2 ... "
    if gatk --java-options "-Xmx12g" Mutect2 \
            -R "${REF}" \
            -I "${OUT}/as_tumour.bam" \
            -I "${OUT}/as_normal.bam" \
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
        echo "failed"
        FAIL=$((FAIL+1)); continue
    fi

    N_RAW=$(zcat "${OUT}/${SID}.raw.vcf.gz" | grep -vc "^#")
    N_PASS=$(zcat "${OUT}/${SID}.filtered.vcf.gz" | grep -v "^#" \
             | awk -F'\t' '$7=="PASS"' | wc -l)
    TOTAL_PASS=$((TOTAL_PASS+N_PASS))

    rm -f "${OUT}/as_tumour.bam" "${OUT}/as_tumour.bam.bai" \
          "${OUT}/as_normal.bam" "${OUT}/as_normal.bam.bai"

    T1=$(date +%s)
    printf "   %d raw, %d PASS   <- every PASS is a false positive   %dm%02ds\n" \
        "${N_RAW}" "${N_PASS}" $(( (T1-T0)/60 )) $(( (T1-T0)%60 ))
    DONE=$((DONE+1))
done

echo ""
echo "======================================================================"
echo " called ${DONE}   skipped ${SKIP}   failed ${FAIL}"
echo " false positives across all samples: ${TOTAL_PASS}"
echo " finished: $(date '+%H:%M:%S')"
echo "======================================================================"
