#!/bin/bash
# =====================================================================
# STEP 9a - build four parallel collections for measuring LOHHLA.
#
# The design in step 6 answered one question badly. Each sample gave one
# target locus, so sensitivity rested on 20 observations, and the untouched
# loci meant to serve as controls mostly failed to reach the statistics
# themselves - only 6 produced a result. Specificity on six observations is
# barely measured.
#
# Four collections separate the questions:
#
#   A_untouched   nothing removed     every heterozygous locus is a control
#   B_loss20      target kept at 0.20  strong loss
#   C_loss35      target kept at 0.35  the fraction used in step 6
#   D_loss50      target kept at 0.50  weak loss
#
# A gives specificity on data where no loss exists anywhere, so any
# detection is unambiguously false. B, C and D give sensitivity against the
# depth of the loss - the same shape of answer step 3 produced for Mutect2
# against VAF.
#
# Two things are done differently from step 6a.
#
# Reads are removed by name rather than by coordinate. The earlier version
# dropped reads whose start position fell inside the window and separately
# subsampled reads overlapping it, so a read beginning before the window
# and extending into it passed through both paths. Selecting names and
# excluding them everywhere removes the whole fragment, which is also
# closer to what losing a chromosome copy does.
#
# The random seed is fixed, so the three loss levels are nested: a read
# dropped at 0.50 is also dropped at 0.35 and 0.20. The only difference
# between the collections is then how much was removed, with no sampling
# noise on top.
#
# The step 6 directories are left alone so those results stay reproducible.
#
# Usage:
#   sbatch step9a_build_collections.sh <n_brca> <n_ov> [off_brca] [off_ov]
# =====================================================================
#SBATCH --job-name=step9a_build
#SBATCH --partition=cpu-single
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/bwfor/work/ws/fr_os136-immune_escape/simulation/step9_lohhla_panel/logs/step9a_%j.log
#SBATCH --error=/gpfs/bwfor/work/ws/fr_os136-immune_escape/simulation/step9_lohhla_panel/logs/step9a_%j.err
# =====================================================================

N_BRCA="${1:-20}"
N_OV="${2:-20}"
OFF_BRCA="${3:-0}"
OFF_OV="${4:-0}"

WS=/gpfs/bwfor/work/ws/fr_os136-immune_escape
COHORT=${WS}/simulation/cohort
PANEL=${WS}/simulation/step9_lohhla_panel
MANIFEST=${COHORT}/manifest.tsv
CONDA=/home/fr/fr_fr/fr_os136/miniconda3/etc/profile.d/conda.sh
THREADS="${SLURM_CPUS_PER_TASK:-4}"
SEED=42

# hg38 windows, the same ones used for typing and in step 6
declare -A WIN_S=( [A]=29932000 [B]=31343000 [C]=31258000 )
declare -A WIN_E=( [A]=29956000 [B]=31367000 [C]=31282000 )

COLLECTIONS="A_untouched:1.00 B_loss20:0.20 C_loss35:0.35 D_loss50:0.50"

source "${CONDA}"
conda activate bio_work

mkdir -p "${PANEL}/logs" "${PANEL}/results"
DESIGN=${PANEL}/collections.tsv
[ -s "${DESIGN}" ] || printf "collection\tsample\tcohort\ttarget_locus\ttarget_alleles\thet_loci\tkeep_requested\tkeep_realised\treads_before\treads_after\tdepth_A\tdepth_B\tdepth_C\n" > "${DESIGN}"

echo "======================================================================"
echo " STEP 9a - four collections for LOHHLA evaluation"
echo " batch: ${N_BRCA} BRCA (off ${OFF_BRCA}), ${N_OV} OV (off ${OFF_OV})"
echo " seed ${SEED}, so the loss levels are nested subsamples"
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

depth_of () {   # bam locus
    samtools depth -a -r "chr6:${WIN_S[$2]}-${WIN_E[$2]}" "$1" 2>/dev/null \
      | awk '{s+=$3;n++} END {printf "%.1f", (n?s/n:0)}'
}

OK=0; SKIP=0; HOMO=0; FAIL=0

for SID in ${SAMPLES}; do
    COH=$(awk -F'\t' -v s="${SID}" 'NR>1 && $1==s {print $2}' "${MANIFEST}")
    SRC_T=${COHORT}/${SID}/${SID}_tumor.bam
    SRC_N=${COHORT}/${SID}/${SID}_normal.bam
    HLA=${COHORT}/${SID}/optitype/${SID}_result.tsv

    echo "----------------------------------------------------------------------"
    echo " ${SID}  (${COH})"

    if [ ! -s "${SRC_T}" ] || [ ! -s "${HLA}" ]; then
        echo "   tumour BAM or HLA type missing"
        printf "%s\tinput\n" "${SID}" >> "${PANEL}/failures.log"
        FAIL=$((FAIL+1)); continue
    fi

    read -r A1 A2 B1 B2 C1 C2 < <(awk 'NR==2 {print $2,$3,$4,$5,$6,$7}' "${HLA}")

    HET=""
    for spec in "A ${A1} ${A2}" "B ${B1} ${B2}" "C ${C1} ${C2}"; do
        set -- ${spec}
        if [ "$2" != "$3" ] && [ "${2#*\*}" != "$2" ] && [ "${3#*\*}" != "$3" ]; then
            HET="${HET}$1"
        fi
    done

    if [ -z "${HET}" ]; then
        echo "   homozygous at every typed locus - contributes nothing"
        printf "%s\thomozygous_all\n" "${SID}" >> "${PANEL}/skipped.log"
        HOMO=$((HOMO+1)); continue
    fi

    TARGET=${HET:0:1}
    case "${TARGET}" in
        A) T_ALL="${A1}/${A2}" ;;
        B) T_ALL="${B1}/${B2}" ;;
        C) T_ALL="${C1}/${C2}" ;;
    esac
    CONTROLS=$(echo "${HET}" | sed "s/${TARGET}//")
    echo "   heterozygous: ${HET}   target: HLA-${TARGET} (${T_ALL})"
    echo "   controls in the loss collections: ${CONTROLS:--}"
    echo "   controls in A_untouched: ${HET}"

    ALLELES=$(awk 'NR==2 {
      for (i = 2; i <= 7; i++)
        if ($i ~ /\*/) {
          split($i, x, "*"); split(x[2], p, ":")
          printf "hla_%s_%s_%s\n", tolower(x[1]), p[1], p[2]
        }
    }' "${HLA}")

    # names of every read overlapping the target window, gathered once and
    # reused for all three loss levels so the subsamples nest
    NAMES=/tmp/${SID}_names_$$.txt
    samtools view "${SRC_T}" \
        "chr6:${WIN_S[$TARGET]}-${WIN_E[$TARGET]}" 2>/dev/null \
      | cut -f1 | sort -u > "${NAMES}"
    N_NAMES=$(wc -l < "${NAMES}")
    echo "   ${N_NAMES} read names overlap the target window"

    for entry in ${COLLECTIONS}; do
        COL=${entry%%:*}
        KEEP=${entry##*:}
        OUT=${PANEL}/${COL}/${SID}

        if [ -s "${OUT}/${SID}_tumor.bam.bai" ]; then
            echo "   ${COL}: already built"
            SKIP=$((SKIP+1)); continue
        fi

        mkdir -p "${OUT}"
        T0=$(date +%s)

        # the normal is identical everywhere, so it is linked rather than
        # copied four times - 320 MB per sample per collection otherwise
        [ -e "${OUT}/${SID}_normal.bam" ] || ln -s "${SRC_N}" "${OUT}/${SID}_normal.bam"
        [ -e "${OUT}/${SID}_normal.bam.bai" ] || ln -s "${SRC_N}.bai" "${OUT}/${SID}_normal.bam.bai"

        N_BEFORE=$(samtools view -c "${SRC_T}" \
                   "chr6:${WIN_S[$TARGET]}-${WIN_E[$TARGET]}" 2>/dev/null)

        if [ "${COL}" = "A_untouched" ]; then
            cp "${SRC_T}" "${OUT}/${SID}_tumor.bam"
            cp "${SRC_T}.bai" "${OUT}/${SID}_tumor.bam.bai" 2>/dev/null \
                || samtools index -@ "${THREADS}" "${OUT}/${SID}_tumor.bam"
            N_DROP=0
        else
            DROP=/tmp/${SID}_${COL}_drop_$$.txt
            awk -v k="${KEEP}" -v seed="${SEED}" \
                'BEGIN {srand(seed)} {if (rand() > k) print}' \
                "${NAMES}" > "${DROP}"
            N_DROP=$(wc -l < "${DROP}")

            samtools view -h "${SRC_T}" 2>/dev/null \
              | awk -v dropfile="${DROP}" '
                  BEGIN { while ((getline line < dropfile) > 0) d[line] = 1 }
                  /^@/ { print; next }
                  !($1 in d) { print }' \
              | samtools sort -@ "${THREADS}" -m 1G \
                  -o "${OUT}/${SID}_tumor.bam" - 2>/dev/null
            samtools index -@ "${THREADS}" "${OUT}/${SID}_tumor.bam"
            rm -f "${DROP}"
        fi

        if [ ! -s "${OUT}/${SID}_tumor.bam" ]; then
            echo "   ${COL}: FAILED to build"
            printf "%s\t%s\tbuild\n" "${SID}" "${COL}" >> "${PANEL}/failures.log"
            FAIL=$((FAIL+1)); continue
        fi

        echo "${ALLELES}" > "${OUT}/${SID}_alleles.txt"
        printf "Ploidy\ttumorPurity\ttumorPloidy\t\n%s_tumor\t2\t1.0\t2\t\n" \
            "${SID}" > "${OUT}/${SID}_solutions.txt"

        DA=$(depth_of "${OUT}/${SID}_tumor.bam" A)
        DB=$(depth_of "${OUT}/${SID}_tumor.bam" B)
        DC=$(depth_of "${OUT}/${SID}_tumor.bam" C)
        N_AFTER=$(samtools view -c "${OUT}/${SID}_tumor.bam" \
                  "chr6:${WIN_S[$TARGET]}-${WIN_E[$TARGET]}" 2>/dev/null)

        REALISED=$(awk -v a="${N_AFTER}" -v b="${N_BEFORE}" \
                   'BEGIN {printf "%.3f", (b ? a/b : 0)}')

        if [ "${COL}" = "A_untouched" ]; then
            printf "%s\t%s\t%s\tNONE\t-\t%s\t1.000\t%s\t%s\t%s\t%s\t%s\t%s\n" \
                "${COL}" "${SID}" "${COH}" "${HET}" "${REALISED}" \
                "${N_BEFORE}" "${N_AFTER}" "${DA}" "${DB}" "${DC}" >> "${DESIGN}"
        else
            printf "%s\t%s\t%s\tHLA-%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
                "${COL}" "${SID}" "${COH}" "${TARGET}" "${T_ALL}" "${HET}" \
                "${KEEP}" "${REALISED}" "${N_BEFORE}" "${N_AFTER}" \
                "${DA}" "${DB}" "${DC}" >> "${DESIGN}"
        fi

        T1=$(date +%s)
        FLAG=""
        if [ "${COL}" != "A_untouched" ]; then
            OFF=$(awk -v r="${REALISED}" -v k="${KEEP}" \
                  'BEGIN {d = r - k; if (d < 0) d = -d; print (d > 0.05) ? 1 : 0}')
            [ "${OFF}" = "1" ] && FLAG="   <-- off target" \
                && printf "%s\t%s\tratio_%s_vs_%s\n" "${SID}" "${COL}" \
                   "${REALISED}" "${KEEP}" >> "${PANEL}/failures.log"
        fi
        printf "   %-13s kept %s of %s reads (%s)   A %6s B %6s C %6s   %ds%s\n" \
               "${COL}:" "${N_AFTER}" "${N_BEFORE}" "${REALISED}" \
               "${DA}" "${DB}" "${DC}" $((T1-T0)) "${FLAG}"
        OK=$((OK+1))
    done

    rm -f "${NAMES}"
done

echo "======================================================================"
echo " STEP 9a finished $(date '+%F %T')"
echo " built ${OK}   skipped ${SKIP}   homozygous ${HOMO}   failed ${FAIL}"
echo "======================================================================"

echo
echo "collections on disk:"
for entry in ${COLLECTIONS}; do
    COL=${entry%%:*}
    n=$(ls -d ${PANEL}/${COL}/*/ 2>/dev/null | wc -l)
    sz=$(du -sh ${PANEL}/${COL} 2>/dev/null | cut -f1)
    printf "  %-14s %3d samples   %s\n" "${COL}" "${n}" "${sz}"
done

echo
echo "expected observations once 9b has run:"
awk -F'\t' 'NR>1 {
    if ($1 == "A_untouched") { ctrl += length($6) }
    else { tgt++; ctrl += length($6) - 1 }
} END {
    printf "  targets  %d   (loss collections, one locus each)\n", tgt
    printf "  controls %d   (A plus the untouched loci elsewhere)\n", ctrl
}' "${DESIGN}"

echo
echo "realised fractions:"
awk -F'\t' 'NR>1 && $1 != "A_untouched" {
    n[$1]++; s[$1] += $8
    if (min[$1] == "" || $8 < min[$1]) min[$1] = $8
    if ($8 > max[$1]) max[$1] = $8
} END {
    for (k in n) printf "  %-13s median %.3f, range %.3f - %.3f  (n=%d)\n",
        k, s[k]/n[k], min[k], max[k], n[k]
}' "${DESIGN}"
