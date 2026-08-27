#!/bin/bash
# =====================================================================
# STEP 14a - build the same samples at several sequencing depths.
#
# Step 3 measured recall at one depth, the 30 to 46x the 1000 Genomes
# slices happened to provide, and showed that what a caller finds is set
# by how many reads carry the mutation rather than by anything about the
# caller. That relationship is stated there but only demonstrated at a
# single point.
#
# Downsampling the same BAMs to a range of depths turns it into a curve.
# The mutations are identical at every level - the same file, the same
# injected positions, the same allele fractions - so recall against depth
# is measured on one set of variants rather than compared across samples
# that differ in other ways.
#
# Both tumour and normal are thinned by the same factor. Mutect2 uses the
# normal to reject germline variants, so leaving it at full depth would
# measure something that does not happen in practice.
#
# The subsampling seed is fixed, which makes the levels nested: a read
# present at 10x is present at every higher level. Recall differences
# between levels are then coverage differences and nothing else.
#
# Usage:
#   sbatch step14a_build_coverage.sh <n_brca> <n_ov> [off_b] [off_o] [levels]
#
#   sbatch step14a_build_coverage.sh 1 0
#   sbatch step14a_build_coverage.sh 20 20 0 0 "10 20 30"
# =====================================================================
#SBATCH --job-name=step14a_cov
#SBATCH --partition=cpu-single
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=16:00:00
#SBATCH --output=/gpfs/bwfor/work/ws/fr_os136-immune_escape/simulation/step14_coverage_sweep/logs/step14a_%j.log
#SBATCH --error=/gpfs/bwfor/work/ws/fr_os136-immune_escape/simulation/step14_coverage_sweep/logs/step14a_%j.err
# =====================================================================

N_BRCA="${1:-20}"
N_OV="${2:-20}"
OFF_BRCA="${3:-0}"
OFF_OV="${4:-0}"
LEVELS="${5:-10 15 20 30 40}"

WS=/gpfs/bwfor/work/ws/fr_os136-immune_escape
COHORT=${WS}/simulation/cohort
SWEEP=${WS}/simulation/step14_coverage_sweep
MANIFEST=${COHORT}/manifest.tsv
PANEL=${WS}/simulation/panel/panel.bed
CONDA=/home/fr/fr_fr/fr_os136/miniconda3/etc/profile.d/conda.sh
THREADS="${SLURM_CPUS_PER_TASK:-4}"
SEED=42

source "${CONDA}"
conda activate bio_work

mkdir -p "${SWEEP}/logs" "${SWEEP}/results"
DESIGN=${SWEEP}/levels.tsv
[ -s "${DESIGN}" ] || printf "level\tsample\tcohort\tdepth_target\tdepth_source\tfraction\tdepth_normal\tdepth_tumour\treads_tumour\n" > "${DESIGN}"

echo "======================================================================"
echo " STEP 14a - downsample to fixed depths"
echo " levels: ${LEVELS}x"
echo " batch: ${N_BRCA} BRCA (off ${OFF_BRCA}), ${N_OV} OV (off ${OFF_OV})"
echo " seed ${SEED}, so the levels are nested"
echo " started $(date '+%F %T') on $(hostname)"
echo "======================================================================"

SAMPLES=$(awk -F'\t' -v nb="${N_BRCA}" -v no="${N_OV}" \
                     -v ob="${OFF_BRCA}" -v oo="${OFF_OV}" '
    NR>1 && $2=="brca" { b++; if (b > ob && b <= ob+nb) print $1 }
    NR>1 && $2=="ov"   { o++; if (o > oo && o <= oo+no) print $1 }
' "${MANIFEST}")

echo
echo "samples: $(echo ${SAMPLES} | tr '\n' ' ')"
echo

panel_depth () {   # bam
    samtools depth -a -b "${PANEL}" "$1" 2>/dev/null \
      | awk '{s+=$3;n++} END {printf "%.2f", (n?s/n:0)}'
}

OK=0; SKIP=0; FAIL=0; TOOLOW=0

for SID in ${SAMPLES}; do
    COH=$(awk -F'\t' -v s="${SID}" 'NR>1 && $1==s {print $2}' "${MANIFEST}")
    SRC_T=${COHORT}/${SID}/${SID}_tumor.bam
    SRC_N=${COHORT}/${SID}/${SID}_normal.bam
    TRUTH=${COHORT}/${SID}/truth_set.tsv

    echo "----------------------------------------------------------------------"
    echo " ${SID}  (${COH})"

    if [ ! -s "${SRC_T}" ] || [ ! -s "${SRC_N}" ]; then
        echo "   source BAMs missing"
        FAIL=$((FAIL+1)); continue
    fi

    echo -n "   measuring source depth ... "
    D_SRC=$(panel_depth "${SRC_T}")
    echo "${D_SRC}x"

    for LEVEL in ${LEVELS}; do
        OUT=${SWEEP}/cov${LEVEL}x/${SID}

        if [ -s "${OUT}/${SID}_tumor.bam.bai" ]; then
            echo "   ${LEVEL}x: already built"
            SKIP=$((SKIP+1)); continue
        fi

        # a level above the source depth cannot be produced by
        # downsampling, and padding it would be inventing reads
        ABOVE=$(awk -v l="${LEVEL}" -v d="${D_SRC}" 'BEGIN {print (l > d) ? 1 : 0}')
        if [ "${ABOVE}" = "1" ]; then
            echo "   ${LEVEL}x: above the source depth of ${D_SRC}x - skipping"
            printf "%s\t%sx\tabove_source_%s\n" "${SID}" "${LEVEL}" "${D_SRC}" \
                >> "${SWEEP}/skipped.log"
            TOOLOW=$((TOOLOW+1)); continue
        fi

        FRAC=$(awk -v l="${LEVEL}" -v d="${D_SRC}" 'BEGIN {printf "%.4f", l/d}')
        mkdir -p "${OUT}"
        T0=$(date +%s)

        # samtools -s takes SEED.FRACTION as one number, so the same seed
        # across levels means the retained reads nest
        for KIND in normal tumor; do
            SRC=${COHORT}/${SID}/${SID}_${KIND}.bam
            samtools view -b -@ "${THREADS}" \
                -s "${SEED}${FRAC#0}" \
                -o "${OUT}/${SID}_${KIND}.bam" "${SRC}" 2>/dev/null
            samtools index -@ "${THREADS}" "${OUT}/${SID}_${KIND}.bam"
        done

        if [ ! -s "${OUT}/${SID}_tumor.bam" ]; then
            echo "   ${LEVEL}x: FAILED"
            printf "%s\t%sx\tbuild\n" "${SID}" "${LEVEL}" >> "${SWEEP}/failures.log"
            FAIL=$((FAIL+1)); continue
        fi

        # the truth set is copied rather than referenced, so each level is
        # self-contained and can be handed to step 14b on its own
        [ -s "${TRUTH}" ] && cp "${TRUTH}" "${OUT}/truth_set.tsv"

        D_N=$(panel_depth "${OUT}/${SID}_normal.bam")
        D_T=$(panel_depth "${OUT}/${SID}_tumor.bam")
        N_READS=$(samtools view -c -@ "${THREADS}" "${OUT}/${SID}_tumor.bam")

        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "cov${LEVEL}x" "${SID}" "${COH}" "${LEVEL}" "${D_SRC}" \
            "${FRAC}" "${D_N}" "${D_T}" "${N_READS}" >> "${DESIGN}"

        T1=$(date +%s)
        OFF=$(awk -v a="${D_T}" -v l="${LEVEL}" \
              'BEGIN {d = a - l; if (d < 0) d = -d; print (d > 2) ? "  <-- off" : ""}')
        printf "   %-6s fraction %-8s normal %6sx  tumour %6sx  %ds%s\n" \
               "${LEVEL}x:" "${FRAC}" "${D_N}" "${D_T}" $((T1-T0)) "${OFF}"
        OK=$((OK+1))
    done
done

echo "======================================================================"
echo " STEP 14a finished $(date '+%F %T')"
echo " built ${OK}   skipped ${SKIP}   above source ${TOOLOW}   failed ${FAIL}"
echo "======================================================================"

echo
echo "on disk:"
for LEVEL in ${LEVELS}; do
    n=$(ls -d ${SWEEP}/cov${LEVEL}x/*/ 2>/dev/null | wc -l)
    sz=$(du -sh ${SWEEP}/cov${LEVEL}x 2>/dev/null | cut -f1)
    printf "  %-8s %3d samples  %s\n" "cov${LEVEL}x" "${n}" "${sz:-0}"
done

echo
echo "realised depths:"
awk -F'\t' 'NR>1 {n[$1]++; s[$1]+=$8
    if (min[$1]=="" || $8<min[$1]) min[$1]=$8
    if ($8>max[$1]) max[$1]=$8}
  END {for (k in n) printf "  %-8s target %sx  realised %.1fx  (%.1f - %.1f, n=%d)\n",
       k, substr(k,4,length(k)-4), s[k]/n[k], min[k], max[k], n[k]}' "${DESIGN}" | sort
