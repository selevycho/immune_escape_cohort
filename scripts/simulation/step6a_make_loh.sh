#!/bin/bash
# =====================================================================
# STEP 6a - simulate HLA loss, exactly as it was done for S001.
#
# The S001 run worked: HLA-A coverage was thinned to 35% while B and C
# were left alone, and LOHHLA reported copy numbers 0.689 against 0.201
# with PVal 0.008. Nothing about the method is changed here - only the
# choice of which locus to hit, which now follows each sample's own
# genotype rather than always being HLA-A.
#
# Thinning is locus-wide: samtools view -s drops a random fraction of the
# reads across the whole locus window. An allele-specific version was
# tried and abandoned - on panel data the two alleles of a locus cannot be
# separated cleanly, reads split 152 against 20 where 50/50 was needed.
#
# The locus to hit is the first heterozygous one in the order A, B, C.
# A homozygous locus has no second allele to lose, so there is nothing to
# simulate there; the untouched loci serve as within-sample controls.
#
# Every sample gets its own directory with private copies of both BAMs.
# LOHHLA writes many intermediates beside its input and builds working
# directory names from the BAM filename, so directories cannot be shared.
#
# Usage:
#   sbatch step6a_make_loh.sh <n_brca> <n_ov> [off_brca] [off_ov] [keep]
#
#   sbatch step6a_make_loh.sh 20 20            keep 35%, as for S001
#   sbatch step6a_make_loh.sh 20 20 0 0 0.5    keep half instead
# =====================================================================
#SBATCH --job-name=step6a_loh
#SBATCH --partition=cpu-single
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/bwfor/work/ws/fr_os136-immune_escape/logs/step6a_%j.log
#SBATCH --error=/gpfs/bwfor/work/ws/fr_os136-immune_escape/logs/step6a_%j.err
# =====================================================================

N_BRCA="${1:-20}"
N_OV="${2:-20}"
OFF_BRCA="${3:-0}"
OFF_OV="${4:-0}"
KEEP="${5:-0.35}"

WS=/gpfs/bwfor/work/ws/fr_os136-immune_escape
MANIFEST=${WS}/simulation/cohort/manifest.tsv
COHORT_DIR=${WS}/simulation/cohort
LOH_DIR=${WS}/simulation/lohhla_allelic
CONDA=/home/fr/fr_fr/fr_os136/miniconda3/etc/profile.d/conda.sh
THREADS="${SLURM_CPUS_PER_TASK:-4}"

# hg38 windows, the same ones used for typing and for S001
A_S=29932000; A_E=29956000
B_S=31343000; B_E=31367000
C_S=31258000; C_E=31282000

source "${CONDA}"
conda activate bio_work

echo "======================================================================"
echo " STEP 6a - simulate loss of one class I locus"
echo " batch: ${N_BRCA} BRCA (off ${OFF_BRCA}), ${N_OV} OV (off ${OFF_OV})"
echo " keeping ${KEEP} of the target locus, other loci untouched"
echo " started $(date '+%F %T') on $(hostname)"
echo "======================================================================"

mkdir -p "${LOH_DIR}"
DESIGN=${LOH_DIR}/loh_design.tsv
[ -s "${DESIGN}" ] || printf "sample\tcohort\ttarget_locus\ttarget_alleles\tcontrols\tkeep_fraction\tdepth_A_before\tdepth_A_after\tdepth_B_before\tdepth_B_after\tdepth_C_before\tdepth_C_after\n" > "${DESIGN}"

SAMPLES=$(awk -F'\t' -v nb="${N_BRCA}" -v no="${N_OV}" \
                     -v ob="${OFF_BRCA}" -v oo="${OFF_OV}" '
    NR>1 && $2=="brca" { b++; if (b > ob && b <= ob+nb) print $1 }
    NR>1 && $2=="ov"   { o++; if (o > oo && o <= oo+no) print $1 }
' "${MANIFEST}")

echo
echo "samples: $(echo ${SAMPLES} | tr '\n' ' ')"
echo

OK=0; SKIP=0; FAIL=0; WAIT=0; HOMO=0

depth_of () {   # $1 bam  $2 start  $3 end
    samtools depth -a -r "chr6:${2}-${3}" "$1" 2>/dev/null \
      | awk '{s+=$3;n++} END {printf "%.1f", (n?s/n:0)}'
}

for SID in ${SAMPLES}; do
    COHORT=$(awk -F'\t' -v s="${SID}" 'NR>1 && $1==s {print $2}' "${MANIFEST}")
    SRC=${COHORT_DIR}/${SID}
    NORMAL=${SRC}/${SID}_normal.bam
    TUMOR=${SRC}/${SID}_tumor.bam
    HLA=${SRC}/optitype/${SID}_result.tsv
    OUT=${LOH_DIR}/${SID}

    echo "----------------------------------------------------------------------"
    echo " ${SID}  (${COHORT})"

    if [ ! -s "${TUMOR}" ] || [ ! -s "${HLA}" ]; then
        echo "   tumour BAM or HLA type not ready"
        WAIT=$((WAIT+1)); continue
    fi

    if [ -s "${OUT}/${SID}_tumor_LOH.bam.bai" ]; then
        echo "   already built - skipping"
        SKIP=$((SKIP+1)); continue
    fi

    read -r A1 A2 B1 B2 C1 C2 < <(awk 'NR==2 {print $2,$3,$4,$5,$6,$7}' "${HLA}")
    printf "   A: %-9s %-9s   B: %-9s %-9s   C: %-9s %-9s\n" \
           "${A1}" "${A2}" "${B1}" "${B2}" "${C1}" "${C2}"

    # ---- first heterozygous locus wins, in the order A, B, C ----
    TARGET=""; T_S=""; T_E=""; T_ALLELES=""; CONTROLS=""
    if   [ "${A1}" != "${A2}" ] && [ "${A1#*\*}" != "${A1}" ] && [ "${A2#*\*}" != "${A2}" ]; then
        TARGET=A; T_S=${A_S}; T_E=${A_E}; T_ALLELES="${A1}/${A2}"; CONTROLS="B,C"
    elif [ "${B1}" != "${B2}" ] && [ "${B1#*\*}" != "${B1}" ] && [ "${B2#*\*}" != "${B2}" ]; then
        TARGET=B; T_S=${B_S}; T_E=${B_E}; T_ALLELES="${B1}/${B2}"; CONTROLS="A,C"
    elif [ "${C1}" != "${C2}" ] && [ "${C1#*\*}" != "${C1}" ] && [ "${C2#*\*}" != "${C2}" ]; then
        TARGET=C; T_S=${C_S}; T_E=${C_E}; T_ALLELES="${C1}/${C2}"; CONTROLS="A,B"
    else
        echo "   homozygous at every typed locus - nothing to lose"
        printf "%s\thomozygous_all\n" "${SID}" >> "${LOH_DIR}/skipped.log"
        HOMO=$((HOMO+1)); continue
    fi
    echo "   target: HLA-${TARGET} (${T_ALLELES})   controls: ${CONTROLS}"

    mkdir -p "${OUT}"
    T0=$(date +%s)

    # ---- private copies; the cohort files are never modified ----
    if [ ! -s "${OUT}/${SID}_normal.bam.bai" ]; then
        echo "   copying normal BAM ..."
        cp "${NORMAL}" "${OUT}/${SID}_normal.bam"
        cp "${NORMAL}.bai" "${OUT}/${SID}_normal.bam.bai" 2>/dev/null \
            || samtools index -@ "${THREADS}" "${OUT}/${SID}_normal.bam"
    fi
    if [ ! -s "${OUT}/${SID}_tumor_original.bam.bai" ]; then
        echo "   copying tumour BAM (pre-loss reference) ..."
        cp "${TUMOR}" "${OUT}/${SID}_tumor_original.bam"
        cp "${TUMOR}.bai" "${OUT}/${SID}_tumor_original.bam.bai" 2>/dev/null \
            || samtools index -@ "${THREADS}" "${OUT}/${SID}_tumor_original.bam"
    fi
    SRC_BAM=${OUT}/${SID}_tumor_original.bam

    DA0=$(depth_of "${SRC_BAM}" ${A_S} ${A_E})
    DB0=$(depth_of "${SRC_BAM}" ${B_S} ${B_E})
    DC0=$(depth_of "${SRC_BAM}" ${C_S} ${C_E})

    # ---- everything outside the target, plus a thinned target, merged ----
    echo "   thinning HLA-${TARGET} to ${KEEP} ..."
    samtools view -h "${SRC_BAM}" \
      | awk -v s="${T_S}" -v e="${T_E}" \
            'BEGIN{OFS="\t"} /^@/ {print; next} !($3=="chr6" && $4>=s && $4<=e) {print}' \
      | samtools view -b -o "${OUT}/rest.bam" - 2>/dev/null

    samtools view -b -s "${KEEP}" "${SRC_BAM}" "chr6:${T_S}-${T_E}" \
        > "${OUT}/target_thin.bam" 2>/dev/null

    samtools merge -f -@ "${THREADS}" -o "${OUT}/merged.bam" \
        "${OUT}/rest.bam" "${OUT}/target_thin.bam" 2>/dev/null
    samtools sort -@ "${THREADS}" -m 1G \
        -o "${OUT}/${SID}_tumor_LOH.bam" "${OUT}/merged.bam" 2>/dev/null
    samtools index -@ "${THREADS}" "${OUT}/${SID}_tumor_LOH.bam"
    rm -f "${OUT}/rest.bam" "${OUT}/target_thin.bam" "${OUT}/merged.bam"

    if [ ! -s "${OUT}/${SID}_tumor_LOH.bam" ]; then
        echo "   FAILED building the LOH BAM"
        printf "%s\tbuild\n" "${SID}" >> "${LOH_DIR}/failures.log"
        FAIL=$((FAIL+1)); continue
    fi

    LOHBAM=${OUT}/${SID}_tumor_LOH.bam
    DA1=$(depth_of "${LOHBAM}" ${A_S} ${A_E})
    DB1=$(depth_of "${LOHBAM}" ${B_S} ${B_E})
    DC1=$(depth_of "${LOHBAM}" ${C_S} ${C_E})

    # ---- LOHHLA inputs: all six alleles, as in the S001 run ----
    awk 'NR==2 {
      for (i = 2; i <= 7; i++) {
        if ($i ~ /\*/) {
          split($i, x, "*"); split(x[2], p, ":")
          printf "hla_%s_%s_%s\n", tolower(x[1]), p[1], p[2]
        }
      }
    }' "${HLA}" > "${OUT}/${SID}_alleles.txt"

    # purity and ploidy are exact by construction in a simulation
    printf "Ploidy\ttumorPurity\ttumorPloidy\t\n%s_tumor_LOH\t2\t1.0\t2\t\n" \
        "${SID}" > "${OUT}/${SID}_solutions.txt"

    printf "%s\t%s\tHLA-%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "${SID}" "${COHORT}" "${TARGET}" "${T_ALLELES}" "${CONTROLS}" "${KEEP}" \
        "${DA0}" "${DA1}" "${DB0}" "${DB1}" "${DC0}" "${DC1}" >> "${DESIGN}"

    T1=$(date +%s)
    printf "   built in %dm%02ds  |  %s\n" $(( (T1-T0)/60 )) $(( (T1-T0)%60 )) \
           "$(ls -lh ${LOHBAM} | awk '{print $5}')"
    printf "   HLA-A  %7s -> %7s x%s\n" "${DA0}" "${DA1}" \
           "$([ "${TARGET}" = "A" ] && echo '   <- target' || echo '   control')"
    printf "   HLA-B  %7s -> %7s x%s\n" "${DB0}" "${DB1}" \
           "$([ "${TARGET}" = "B" ] && echo '   <- target' || echo '   control')"
    printf "   HLA-C  %7s -> %7s x%s\n" "${DC0}" "${DC1}" \
           "$([ "${TARGET}" = "C" ] && echo '   <- target' || echo '   control')"

    OK=$((OK+1))
done

echo "======================================================================"
echo " STEP 6a finished $(date '+%F %T')"
echo " built ${OK}   skipped ${SKIP}   homozygous ${HOMO}   waiting ${WAIT}   failed ${FAIL}"
echo "======================================================================"
echo
echo "design:"
column -t -s $'\t' "${DESIGN}" 2>/dev/null | sed 's/^/  /' || cat "${DESIGN}"
echo
echo "target loci:"
awk -F'\t' 'NR>1 {c[$3]++} END {for (k in c) printf "  %-8s %d\n", k, c[k]}' "${DESIGN}"
