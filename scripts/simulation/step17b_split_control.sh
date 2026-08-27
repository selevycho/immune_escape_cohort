#!/bin/bash
#SBATCH --job-name=step17b_split
#SBATCH --partition=cpu-single
#SBATCH --time=08:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --output=%x_%j.log
#
# Step 17b - the negative control, done properly.
#
# The first attempt passed one normal BAM to Mutect2 twice, under two
# sample names. That produced 376 PASS calls across 33 samples, none of
# which appears in the real run, and the records show why: tumour and
# normal carry byte-identical read counts and Mutect2 still assigns them
# different genotypes. With the same read names on both sides the caller
# has no independent evidence to compare and its statistics break.
#
# Splitting the reads fixes it. Each read goes to one side or the other by
# a hash of its name, so both halves are real independent sequencing of a
# person with no somatic variants. Any PASS call is then a genuine false
# positive, at half the depth - which makes this a conservative estimate
# rather than an optimistic one.
#
# Usage:
#   sbatch step17b_split_control.sh              all forty
#   sbatch step17b_split_control.sh B001 B002    named samples
set -o pipefail

WS=$(ws_find immune_escape)
COHORT="${WS}/simulation/cohort"
OUTDIR="${WS}/simulation/step17b_split_control"
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
echo " STEP 17b - negative control by splitting the reads"
echo "======================================================================"
echo "samples: $(echo ${SAMPLES} | tr '\n' ' ')"
echo "started: $(date '+%H:%M:%S')"
echo "----------------------------------------------------------------------"

DONE=0; SKIP=0; FAIL=0; TOTAL=0

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
        TOTAL=$((TOTAL+N)); SKIP=$((SKIP+1)); continue
    fi

    mkdir -p "${OUT}"
    T0=$(date +%s)

    # a read and its mate must stay together, so the side is chosen from
    # the read name rather than at random per record
    echo -n "   splitting ... "
    samtools view -h "${NBAM}" | awk -v out="${OUT}" '
        /^@/ { print > (out "/half_a.sam"); print > (out "/half_b.sam"); next }
        {
            h = 0
            n = split($1, ch, "")
            for (i = 1; i <= n; i++) h = (h * 31 + index("ACGTacgt0123456789:_-.", ch[i])) % 997
            if (h % 2 == 0) print > (out "/half_a.sam")
            else            print > (out "/half_b.sam")
        }' || { echo "failed"; FAIL=$((FAIL+1)); continue; }

    samtools addreplacerg -@ "${THREADS}" \
        -r "ID:TUMOR\tSM:${SID}_TUMOR\tPL:ILLUMINA\tLB:panel\tPU:sim" \
        -o "${OUT}/as_tumour.bam" "${OUT}/half_a.sam" 2>"${OUT}/rg.err" \
        && samtools index -@ "${THREADS}" "${OUT}/as_tumour.bam"
    samtools addreplacerg -@ "${THREADS}" \
        -r "ID:NORMAL\tSM:${SID}_NORMAL\tPL:ILLUMINA\tLB:panel\tPU:sim" \
        -o "${OUT}/as_normal.bam" "${OUT}/half_b.sam" 2>>"${OUT}/rg.err" \
        && samtools index -@ "${THREADS}" "${OUT}/as_normal.bam"
    rm -f "${OUT}/half_a.sam" "${OUT}/half_b.sam"

    DA=$(samtools view -c "${OUT}/as_tumour.bam")
    DB=$(samtools view -c "${OUT}/as_normal.bam")
    echo "ok  (${DA} / ${DB} reads)"

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
        echo "failed"; FAIL=$((FAIL+1)); continue
    fi

    echo -n "   FilterMutectCalls ... "
    if gatk --java-options "-Xmx8g" FilterMutectCalls \
            -R "${REF}" -V "${OUT}/${SID}.raw.vcf.gz" \
            -O "${OUT}/${SID}.filtered.vcf.gz" \
            > "${OUT}/filter.log" 2>&1; then
        echo "ok"
    else
        echo "failed"; FAIL=$((FAIL+1)); continue
    fi

    N_RAW=$(zcat "${OUT}/${SID}.raw.vcf.gz" | grep -vc "^#")
    N_PASS=$(zcat "${OUT}/${SID}.filtered.vcf.gz" | grep -v "^#" \
             | awk -F'\t' '$7=="PASS"' | wc -l)
    TOTAL=$((TOTAL+N_PASS))

    rm -f "${OUT}/as_tumour.bam"* "${OUT}/as_normal.bam"*

    T1=$(date +%s)
    printf "   %d raw, %d PASS   %dm%02ds\n" "${N_RAW}" "${N_PASS}" \
        $(( (T1-T0)/60 )) $(( (T1-T0)%60 ))
    DONE=$((DONE+1))
done

echo ""
echo "======================================================================"
echo " called ${DONE}   skipped ${SKIP}   failed ${FAIL}"
echo " false positives in total: ${TOTAL}"
echo "======================================================================"
