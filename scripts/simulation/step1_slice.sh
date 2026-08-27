#!/bin/bash
# =====================================================================
# STEP 1 - fetch the backbone reads for a batch of samples.
#
# For each sample in the batch this slices the panel out of its 1000
# Genomes CRAM, sorts and indexes the result, and reports coverage.
#
# The CRAM is read remotely rather than downloaded: the full file is about
# 20 GB, the panel is 24 Mb, and samtools uses the .crai index to fetch only
# the byte ranges that matter. A whole-file download would take an hour per
# sample; this takes minutes.
#
# Samples already sliced are skipped, so the script is safe to re-run after
# a partial failure.
#
# Usage:
#   step1_slice.sh <n_brca> <n_ov> [offset_brca] [offset_ov]
#
#   step1_slice.sh 5 5        first 5 of each cohort
#   step1_slice.sh 5 5 5 5    the next 5 of each
#
# Submit with:
#   sbatch step1_slice.sh 5 5
# =====================================================================
#SBATCH --job-name=step1_slice
#SBATCH --partition=cpu-single
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=/gpfs/bwfor/work/ws/fr_os136-immune_escape/logs/step1_%j.log
#SBATCH --error=/gpfs/bwfor/work/ws/fr_os136-immune_escape/logs/step1_%j.err
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
CONDA=/home/fr/fr_fr/fr_os136/miniconda3/etc/profile.d/conda.sh

MIN_DEPTH=10
THREADS="${SLURM_CPUS_PER_TASK:-4}"

source "${CONDA}"
conda activate bio_work

echo "======================================================================"
echo " STEP 1 - slice panel from remote CRAM"
echo " batch: ${N_BRCA} BRCA (offset ${OFF_BRCA}), ${N_OV} OV (offset ${OFF_OV})"
echo " started $(date '+%F %T') on $(hostname)"
echo " samtools $(samtools --version | head -1 | awk '{print $2}')"
echo "======================================================================"

for f in "${MANIFEST}" "${PANEL}" "${REF}"; do
    [ -s "$f" ] || { echo "ERROR: missing $f"; exit 1; }
done

# ---- pick the samples for this batch ----
SAMPLES=$(awk -F'\t' -v nb="${N_BRCA}" -v no="${N_OV}" \
                     -v ob="${OFF_BRCA}" -v oo="${OFF_OV}" '
    NR>1 && $2=="brca" { b++; if (b > ob && b <= ob+nb) print $1 }
    NR>1 && $2=="ov"   { o++; if (o > oo && o <= oo+no) print $1 }
' "${MANIFEST}")

N_TOTAL=$(echo "${SAMPLES}" | grep -c . )
echo
echo "samples in this batch (${N_TOTAL}):"
echo "${SAMPLES}" | tr '\n' ' '
echo
echo

OK=0
SKIP=0
FAIL=0
START_ALL=$(date +%s)

for SID in ${SAMPLES}; do
    LINE=$(awk -F'\t' -v s="${SID}" 'NR>1 && $1==s' "${MANIFEST}")
    COHORT=$(echo "${LINE}"   | cut -f2)
    BARCODE=$(echo "${LINE}"  | cut -f3)
    NMUT=$(echo "${LINE}"     | cut -f4)
    BACKBONE=$(echo "${LINE}" | cut -f6)
    SUPERPOP=$(echo "${LINE}" | cut -f8)
    CRAM=$(echo "${LINE}"     | cut -f10)

    OUT=${COHORT_DIR}/${SID}
    BAM=${OUT}/${SID}_normal.bam
    mkdir -p "${OUT}"

    echo "----------------------------------------------------------------------"
    printf " %-6s %-5s %-16s backbone %-9s %-4s  ~%s mutations\n" \
           "${SID}" "${COHORT}" "${BARCODE}" "${BACKBONE}" "${SUPERPOP}" "${NMUT}"

    if [ -s "${BAM}" ] && [ -s "${BAM}.bai" ]; then
        D=$(samtools depth -a -b "${PANEL}" "${BAM}" 2>/dev/null \
            | awk '{s+=$3; n++} END {if(n>0) printf "%.1f", s/n; else print 0}')
        echo "   already done, ${D}x - skipping"
        SKIP=$((SKIP+1))
        continue
    fi

    T0=$(date +%s)
    echo "   $(date '+%T') slicing ..."

    # EBI stalls under parallel load and samtools has no timeout of its
    # own, so without this a hung fetch holds the slot indefinitely
    if ! timeout 45m samtools view -@ "${THREADS}" -b -T "${REF}" -L "${PANEL}" \
            -o "${OUT}/tmp_unsorted.bam" "${CRAM}" 2>"${OUT}/slice.err"; then
        echo "   FAILED at the slice step:"
        tail -3 "${OUT}/slice.err" | sed 's/^/     /'
        echo "${SID}\tslice" >> "${COHORT_DIR}/step1_failures.log"
        FAIL=$((FAIL+1))
        rm -f "${OUT}/tmp_unsorted.bam"
        continue
    fi

    if ! samtools sort -@ "${THREADS}" -m 1G -o "${BAM}" \
            "${OUT}/tmp_unsorted.bam" 2>>"${OUT}/slice.err"; then
        echo "   FAILED at sort"
        echo "${SID}\tsort" >> "${COHORT_DIR}/step1_failures.log"
        FAIL=$((FAIL+1))
        rm -f "${OUT}/tmp_unsorted.bam" "${BAM}"
        continue
    fi
    samtools index -@ "${THREADS}" "${BAM}"
    rm -f "${OUT}/tmp_unsorted.bam"

    T1=$(date +%s)
    SIZE=$(ls -lh "${BAM}" | awk '{print $5}')
    READS=$(samtools view -c -@ "${THREADS}" "${BAM}")
    D=$(samtools depth -a -b "${PANEL}" "${BAM}" \
        | awk '{s+=$3; n++} END {if(n>0) printf "%.1f", s/n; else print 0}')

    printf "   done in %dm%02ds  |  %s  |  %s reads  |  %sx mean depth\n" \
           $(( (T1-T0)/60 )) $(( (T1-T0)%60 )) "${SIZE}" "${READS}" "${D}"

    if [ "${D%.*}" -lt "${MIN_DEPTH}" ]; then
        echo "   WARNING: coverage below ${MIN_DEPTH}x - this sample may not be usable"
        echo "${SID}\tlow_coverage_${D}" >> "${COHORT_DIR}/step1_failures.log"
    fi
    OK=$((OK+1))
done

END_ALL=$(date +%s)
echo "======================================================================"
echo " STEP 1 finished $(date '+%F %T')"
printf " elapsed: %dh%02dm    sliced: %d   skipped: %d   failed: %d\n" \
       $(( (END_ALL-START_ALL)/3600 )) $(( ((END_ALL-START_ALL)%3600)/60 )) \
       "${OK}" "${SKIP}" "${FAIL}"
echo "======================================================================"
echo
echo "current state of the cohort:"
printf " %-6s %-6s %10s %10s\n" "sample" "cohort" "size" "depth"
awk -F'\t' 'NR>1 {print $1"\t"$2}' "${MANIFEST}" | while read -r s c; do
    b=${COHORT_DIR}/${s}/${s}_normal.bam
    if [ -s "$b" ]; then
        sz=$(ls -lh "$b" | awk '{print $5}')
        printf " %-6s %-6s %10s %10s\n" "$s" "$c" "$sz" "ready"
    fi
done
