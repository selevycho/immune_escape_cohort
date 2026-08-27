#!/bin/bash
# =====================================================================
# STEP 5 - neoantigen prediction with mhcflurry.
#
# Takes the mutations that were injected into each sample, translates them
# into mutant peptides, and predicts binding against that sample's own HLA
# genotype as called by OptiType.
#
# One caveat worth stating plainly: the HLA genotype belongs to the 1000
# Genomes backbone, while the mutations belong to a TCGA patient. The
# predictions are therefore "what these mutations would produce in someone
# with this genotype", not that patient's real neoantigens. For measuring
# what the pipeline can detect this is fine; for biological claims about
# the donor it is not.
#
# mhcflurry lives in its own environment because tensorflow does not
# coexist with the rest of the stack.
#
# Usage:
#   sbatch step5_neoantigens.sh <n_brca> <n_ov> [offset_brca] [offset_ov] [lengths]
#
#   sbatch step5_neoantigens.sh 20 20              everything that is ready
#   sbatch step5_neoantigens.sh 5 5 0 0 9,10       9-mers and 10-mers only
# =====================================================================
#SBATCH --job-name=step5_neo
#SBATCH --partition=cpu-single
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=08:00:00
#SBATCH --output=/gpfs/bwfor/work/ws/fr_os136-immune_escape/logs/step5_%j.log
#SBATCH --error=/gpfs/bwfor/work/ws/fr_os136-immune_escape/logs/step5_%j.err
# =====================================================================
set -o pipefail

N_BRCA="${1:-20}"
N_OV="${2:-20}"
OFF_BRCA="${3:-0}"
OFF_OV="${4:-0}"
LENGTHS="${5:-8,9,10,11}"

WS=/gpfs/bwfor/work/ws/fr_os136-immune_escape
MANIFEST=${WS}/simulation/cohort/manifest.tsv
COHORT_DIR=${WS}/simulation/cohort
PROT=${WS}/ref/gencode.v46.pc_translations.fa.gz
SCRIPTS=/home/fr/fr_fr/fr_os136/immune_escape_project/scripts/01_cohort_tcga
CONDA=/home/fr/fr_fr/fr_os136/miniconda3/etc/profile.d/conda.sh

source "${CONDA}"
conda activate mhc_env

echo "======================================================================"
echo " STEP 5 - neoantigen prediction (mhcflurry)"
echo " batch: ${N_BRCA} BRCA (offset ${OFF_BRCA}), ${N_OV} OV (offset ${OFF_OV})"
echo " peptide lengths: ${LENGTHS}"
echo " started $(date '+%F %T') on $(hostname)"
echo "======================================================================"

[ -s "${PROT}" ] || { echo "ERROR: ${PROT} missing"; exit 1; }
[ -s "${SCRIPTS}/make_neoantigens.py" ] \
    || { echo "ERROR: make_neoantigens.py missing"; exit 1; }
command -v mhcflurry-predict >/dev/null \
    || { echo "ERROR: mhcflurry not available in mhc_env"; exit 1; }

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
    COHORT=$(echo "${LINE}"   | cut -f2)
    BACKBONE=$(echo "${LINE}" | cut -f6)

    OUT=${COHORT_DIR}/${SID}
    TRUTH=${OUT}/truth_set.tsv
    HLA=${OUT}/optitype/${SID}_result.tsv
    NEO=${OUT}/neoantigens

    echo "----------------------------------------------------------------------"
    printf " %-6s %-5s backbone %s\n" "${SID}" "${COHORT}" "${BACKBONE}"

    if [ ! -s "${TRUTH}" ]; then
        echo "   truth set missing - step 2 has not reached this sample"
        WAIT=$((WAIT+1)); continue
    fi
    if [ ! -s "${HLA}" ]; then
        echo "   HLA type missing - step 4 has not reached this sample"
        WAIT=$((WAIT+1)); continue
    fi
    if [ -s "${NEO}/neoantigens_per_mutation.tsv" ]; then
        n=$(awk 'NR>1' "${NEO}/neoantigens_per_mutation.tsv" | wc -l)
        s=$(awk -F'\t' 'NR>1 {t+=$6} END {print t+0}' "${NEO}/neoantigens_per_mutation.tsv")
        echo "   already predicted: ${n} mutations, ${s} strong binders"
        SKIP=$((SKIP+1)); continue
    fi

    T0=$(date +%s)
    echo "   $(date '+%T') predicting ..."

    if ! python "${SCRIPTS}/make_neoantigens.py" \
            "${TRUTH}" "${HLA}" "${PROT}" "${NEO}" "${LENGTHS}" \
            > "${OUT}/neoantigens.log" 2>&1; then
        echo "   FAILED:"
        tail -6 "${OUT}/neoantigens.log" | sed 's/^/     /'
        printf "%s\tneoantigens\n" "${SID}" >> "${COHORT_DIR}/step5_failures.log"
        FAIL=$((FAIL+1)); continue
    fi

    T1=$(date +%s)
    grep -E "missense inputs|peptides tested|strong binders|weak binders|at least one binder" \
         "${OUT}/neoantigens.log" | sed 's/^/   /'
    printf "   done in %dm%02ds\n" $(( (T1-T0)/60 )) $(( (T1-T0)%60 ))
    OK=$((OK+1))
done

END_ALL=$(date +%s)
echo "======================================================================"
echo " STEP 5 finished $(date '+%F %T')"
printf " elapsed %dh%02dm   predicted %d   skipped %d   waiting %d   failed %d\n" \
       $(( (END_ALL-START_ALL)/3600 )) $(( ((END_ALL-START_ALL)%3600)/60 )) \
       "${OK}" "${SKIP}" "${WAIT}" "${FAIL}"
echo "======================================================================"

# ---------------- cohort summary ----------------
echo
echo "neoantigens across the cohort:"
printf " %-6s %-5s %8s %8s %8s %8s %10s\n" \
       "sample" "coh" "mutations" "peptides" "strong" "weak" "yield"

TOT_MUT=0; TOT_STRONG=0; TOT_WEAK=0; N=0
while IFS=$'\t' read -r sid cohort rest; do
    [ "${sid}" = "sample_id" ] && continue
    F=${COHORT_DIR}/${sid}/neoantigens/neoantigens_per_mutation.tsv
    [ -s "${F}" ] || continue
    read -r nm st wk pep < <(awk -F'\t' '
        NR>1 {n++; s+=$6; w+=$7; p+=$4}
        END {print n+0, s+0, w+0, p+0}' "${F}")
    yielded=$(awk -F'\t' 'NR>1 && ($6+$7)>0 {n++} END {print n+0}' "${F}")
    printf " %-6s %-5s %8s %8s %8s %8s %6s/%s\n" \
           "${sid}" "${cohort}" "${nm}" "${pep}" "${st}" "${wk}" "${yielded}" "${nm}"
    TOT_MUT=$((TOT_MUT+nm)); TOT_STRONG=$((TOT_STRONG+st)); TOT_WEAK=$((TOT_WEAK+wk))
    N=$((N+1))
done < "${MANIFEST}"

echo
if [ "${N}" -gt 0 ]; then
    echo "  samples with predictions : ${N}"
    echo "  missense mutations used  : ${TOT_MUT}"
    echo "  strong binders           : ${TOT_STRONG}"
    echo "  weak binders             : ${TOT_WEAK}"
    printf "  strong binders per mutation: %.2f\n" \
           "$(echo "${TOT_STRONG} ${TOT_MUT}" | awk '{print ($2>0 ? $1/$2 : 0)}')"
else
    echo "  nothing predicted yet"
fi
