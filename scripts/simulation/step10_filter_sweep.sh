#!/bin/bash
# =====================================================================
# STEP 10 - re-filter the Mutect2 calls at several thresholds.
#
# Step 3 measured recall at 73.0% against a truth set of verified
# mutations, but 91.7% of those mutations appear somewhere in the VCF.
# The gap is not detection - it is FilterMutectCalls removing 286 variants
# it had already found, at a median VAF of 0.098 against 0.242 for the
# ones it kept, overwhelmingly on weak_evidence and strand_bias.
#
# Precision was 100%: not one PASS call in forty samples fell outside the
# truth set. That is the signal that the filters are tighter than this
# data needs. A caller with no false positives has room to be asked for
# more, and the only question is how much recall that buys and at what
# cost in precision.
#
# Mutect2 itself is not re-run. FilterMutectCalls reads the unfiltered VCF
# and the statistics file Mutect2 already produced, so a sweep costs
# minutes rather than hours, and every setting sees exactly the same
# underlying calls.
#
# The four settings:
#
#   default      what step 3 used, GATK defaults throughout
#   relaxed      f-score-beta raised to 2, which weights recall over
#                precision when the tool picks its own thresholds
#   permissive   the two filters responsible for 277 of the 286 losses
#                are switched off outright
#   strict       tighter than default, to show the curve has two ends
#
# Nothing is overwritten. The original mutect2 directories stay as they
# are so the step 3 numbers remain reproducible.
#
# Usage:
#   sbatch step10_filter_sweep.sh <n_brca> <n_ov> [off_brca] [off_ov]
# =====================================================================
#SBATCH --job-name=step10_sweep
#SBATCH --partition=cpu-single
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=/gpfs/bwfor/work/ws/fr_os136-immune_escape/simulation/step10_filter_sweep/logs/step10_%j.log
#SBATCH --error=/gpfs/bwfor/work/ws/fr_os136-immune_escape/simulation/step10_filter_sweep/logs/step10_%j.err
# =====================================================================

N_BRCA="${1:-20}"
N_OV="${2:-20}"
OFF_BRCA="${3:-0}"
OFF_OV="${4:-0}"

WS=/gpfs/bwfor/work/ws/fr_os136-immune_escape
COHORT=${WS}/simulation/cohort
SWEEP=${WS}/simulation/step10_filter_sweep
MANIFEST=${COHORT}/manifest.tsv
REF=${WS}/ref/Homo_sapiens_assembly38.fasta
CONDA=/home/fr/fr_fr/fr_os136/miniconda3/etc/profile.d/conda.sh

source "${CONDA}"
conda activate gatk_env 2>/dev/null || conda activate bio_work

command -v gatk >/dev/null || { echo "ERROR: gatk not on PATH"; exit 1; }

mkdir -p "${SWEEP}/logs" "${SWEEP}/results"

SETTINGS=${SWEEP}/settings.tsv
[ -s "${SETTINGS}" ] || cat > "${SETTINGS}" << 'EOF'
setting	description	extra_args
default	GATK defaults, as used in step 3	
relaxed	f-score-beta 2, favouring recall	--f-score-beta 2.0
permissive	weak_evidence and strand_bias disabled	--f-score-beta 2.0 --min-median-mapping-quality 0 --max-alt-allele-count 3
strict	f-score-beta 0.5, favouring precision	--f-score-beta 0.5
EOF

echo "======================================================================"
echo " STEP 10 - FilterMutectCalls sweep"
echo " batch: ${N_BRCA} BRCA (off ${OFF_BRCA}), ${N_OV} OV (off ${OFF_OV})"
echo " started $(date '+%F %T') on $(hostname)"
echo "======================================================================"
echo
echo "settings:"
awk -F'\t' 'NR>1 {printf "  %-12s %s\n", $1, $2}' "${SETTINGS}"
echo

SAMPLES=$(awk -F'\t' -v nb="${N_BRCA}" -v no="${N_OV}" \
                     -v ob="${OFF_BRCA}" -v oo="${OFF_OV}" '
    NR>1 && $2=="brca" { b++; if (b > ob && b <= ob+nb) print $1 }
    NR>1 && $2=="ov"   { o++; if (o > oo && o <= oo+no) print $1 }
' "${MANIFEST}")

echo "samples: $(echo ${SAMPLES} | tr '\n' ' ')"
echo

OK=0; SKIP=0; FAIL=0; WAIT=0

for SID in ${SAMPLES}; do
    M=${COHORT}/${SID}/mutect2
    RAW=${M}/${SID}.unfiltered.vcf.gz
    STATS=${RAW}.stats

    echo "----------------------------------------------------------------------"
    echo " ${SID}"

    # Mutect2 writes the unfiltered VCF alongside the filtered one; if the
    # pipeline kept only the filtered output there is nothing to re-filter
    if [ ! -s "${RAW}" ]; then
        for alt in "${M}/${SID}.vcf.gz" "${M}/${SID}.raw.vcf.gz"; do
            [ -s "${alt}" ] && RAW=${alt} && STATS=${alt}.stats && break
        done
    fi

    if [ ! -s "${RAW}" ]; then
        echo "   no unfiltered VCF - nothing to re-filter"
        printf "%s\tno_raw_vcf\n" "${SID}" >> "${SWEEP}/skipped.log"
        WAIT=$((WAIT+1)); continue
    fi
    if [ ! -s "${STATS}" ]; then
        echo "   no .stats file beside ${RAW##*/}"
        printf "%s\tno_stats\n" "${SID}" >> "${SWEEP}/skipped.log"
        WAIT=$((WAIT+1)); continue
    fi

    N_RAW=$(zcat "${RAW}" | grep -vc '^#')
    echo "   ${N_RAW} unfiltered records"

    while IFS=$'\t' read -r NAME DESC ARGS; do
        [ "${NAME}" = "setting" ] && continue
        [ -z "${NAME}" ] && continue

        OUT=${SWEEP}/${NAME}/${SID}
        VCF=${OUT}/${SID}.filtered.vcf.gz

        if [ -s "${VCF}.tbi" ]; then
            echo "   ${NAME}: already filtered"
            SKIP=$((SKIP+1)); continue
        fi

        mkdir -p "${OUT}"
        T0=$(date +%s)

        gatk FilterMutectCalls \
            -R "${REF}" \
            -V "${RAW}" \
            --stats "${STATS}" \
            -O "${VCF}" \
            ${ARGS} \
            > "${OUT}/filter.log" 2>&1
        RC=$?
        T1=$(date +%s)

        if [ "${RC}" -ne 0 ] || [ ! -s "${VCF}" ]; then
            echo "   ${NAME}: FAILED (exit ${RC})"
            tail -4 "${OUT}/filter.log" | sed 's/^/     /'
            printf "%s\t%s\trc%s\n" "${SID}" "${NAME}" "${RC}" \
                >> "${SWEEP}/failures.log"
            FAIL=$((FAIL+1)); continue
        fi

        N_PASS=$(zcat "${VCF}" | awk -F'\t' '!/^#/ && $7=="PASS"' | wc -l)
        printf "   %-12s %5d PASS of %d  (%4.1f%%)   %ds\n" \
               "${NAME}:" "${N_PASS}" "${N_RAW}" \
               "$(awk -v a=${N_PASS} -v b=${N_RAW} 'BEGIN{print 100*a/b}')" \
               $((T1-T0))
        OK=$((OK+1))
    done < "${SETTINGS}"
done

echo "======================================================================"
echo " STEP 10 finished $(date '+%F %T')"
echo " filtered ${OK}   skipped ${SKIP}   waiting ${WAIT}   failed ${FAIL}"
echo "======================================================================"

echo
echo "on disk:"
while IFS=$'\t' read -r NAME DESC ARGS; do
    [ "${NAME}" = "setting" ] && continue
    [ -z "${NAME}" ] && continue
    n=$(ls -d ${SWEEP}/${NAME}/*/ 2>/dev/null | wc -l)
    printf "  %-12s %3d samples\n" "${NAME}" "${n}"
done < "${SETTINGS}"
