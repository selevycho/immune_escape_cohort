#!/bin/bash
# =====================================================================
# STEP 9d - the same collections at a different coverage threshold.
#
# Step 9b ran everything at minCoverageFilter 3, chosen because a single
# test sample gained usable positions at that setting: 1 position at the
# LOHHLA default of 5, 8 at 3, 21 at 2, with the loss detected only at 3.
#
# Across the full panel that choice looks less clear. Two results argue
# against it. Detection stopped tracking the number of discriminating
# positions - the association was p = 0.006 at threshold 5 in step 6b and
# p = 0.23 here - and loci with 10 or more positions detected less often
# than loci with 5 to 9, which cannot be right if more positions mean more
# information. Both are what would happen if the extra positions admitted
# at 3 carry mostly noise: a position covered by three reads gives an
# allele ratio that swings on one read.
#
# Rather than argue about it, the same collections are run again at 5 and
# the two are compared directly. Nothing else changes - same BAMs, same
# allele references, same script - so any difference is the threshold.
#
# Results go to out_<locus>_mc<threshold>, leaving the step 9b output in
# place.
#
# Usage:
#   sbatch step9d_threshold_sweep.sh <collection> <mincov> <n_brca> <n_ov> [ob] [oo]
#
#   sbatch step9d_threshold_sweep.sh A_untouched 5 10 0 0 0
# =====================================================================
#SBATCH --job-name=step9d_mc
#SBATCH --partition=cpu-single
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=16:00:00
#SBATCH --output=/gpfs/bwfor/work/ws/fr_os136-immune_escape/simulation/step9_lohhla_panel/logs/step9d_%j.log
#SBATCH --error=/gpfs/bwfor/work/ws/fr_os136-immune_escape/simulation/step9_lohhla_panel/logs/step9d_%j.err
# =====================================================================

COLLECTION="${1:?collection: A_untouched B_loss20 C_loss35 D_loss50}"
MIN_COV="${2:?minCoverageFilter value, e.g. 5}"
N_BRCA="${3:-20}"
N_OV="${4:-20}"
OFF_BRCA="${5:-0}"
OFF_OV="${6:-0}"

WS=/gpfs/bwfor/work/ws/fr_os136-immune_escape
PANEL=${WS}/simulation/step9_lohhla_panel
COL_DIR=${PANEL}/${COLLECTION}
MANIFEST=${WS}/simulation/cohort/manifest.tsv
LOHHLA_HOME=/home/fr/fr_fr/fr_os136/immune_escape_project/soft/lohhla
HLA_FASTA_ALL=${LOHHLA_HOME}/data/hla_all_lohhla.fasta
HLA_EXON=${LOHHLA_HOME}/data/hla.dat
CONDA=/home/fr/fr_fr/fr_os136/miniconda3/etc/profile.d/conda.sh

source "${CONDA}"
conda activate lohhla_env
export PATH=/home/fr/fr_fr/fr_os136/miniconda3/envs/bio_work/bin:$PATH

echo "======================================================================"
echo " STEP 9d - ${COLLECTION} at minCoverageFilter ${MIN_COV}"
echo " batch: ${N_BRCA} BRCA (off ${OFF_BRCA}), ${N_OV} OV (off ${OFF_OV})"
echo " started $(date '+%F %T') on $(hostname)"
echo "======================================================================"

[ -d "${COL_DIR}" ] || { echo "ERROR: ${COL_DIR} missing"; exit 1; }

SAMPLES=$(awk -F'\t' -v nb="${N_BRCA}" -v no="${N_OV}" \
                     -v ob="${OFF_BRCA}" -v oo="${OFF_OV}" '
    NR>1 && $2=="brca" { b++; if (b > ob && b <= ob+nb) print $1 }
    NR>1 && $2=="ov"   { o++; if (o > oo && o <= oo+no) print $1 }
' "${MANIFEST}")

N_RUN=0; N_OK=0; N_FAIL=0; N_SKIP=0

for SID in ${SAMPLES}; do
    DIR=${COL_DIR}/${SID}
    ALLELES=${DIR}/${SID}_alleles.txt
    [ -s "${DIR}/${SID}_tumor.bam" ] && [ -s "${ALLELES}" ] || continue

    echo "----------------------------------------------------------------------"
    echo " ${SID}"
    cd "${DIR}"

    for LOCUS in a b c; do
        grep "^hla_${LOCUS}_" "${ALLELES}" | sort -u > "${DIR}/al_${LOCUS}.txt"
        [ "$(wc -l < ${DIR}/al_${LOCUS}.txt)" -ge 2 ] || continue

        OUTDIR=${DIR}/out_${LOCUS}_mc${MIN_COV}
        PRED=$(ls "${OUTDIR}"/*HLAlossPrediction*.txt 2>/dev/null | head -1)
        if [ -n "${PRED}" ] && [ "$(wc -l < ${PRED})" -gt 1 ]; then
            echo "   HLA-${LOCUS^^}: already done at mc${MIN_COV}"
            N_SKIP=$((N_SKIP+1)); continue
        fi

        SUB=${DIR}/sub_${LOCUS}.fasta
        if [ ! -s "${SUB}" ]; then
            python3 - "${DIR}/al_${LOCUS}.txt" "${HLA_FASTA_ALL}" "${SUB}" << 'PYEOF'
import sys
keep = tuple({l.strip() for l in open(sys.argv[1]) if l.strip()})
out, w, n = [], False, 0
for line in open(sys.argv[2]):
    if line.startswith(">"):
        w = any(line[1:].strip().startswith(k) for k in keep)
        if w:
            n += 1
    if w:
        out.append(line)
open(sys.argv[3], "w").writelines(out)
PYEOF
        fi
        [ "$(grep -c '^>' ${SUB} 2>/dev/null)" -ge 2 ] || continue

        rm -rf "${OUTDIR}"
        mkdir -p "${OUTDIR}/${SID}_normal" "${OUTDIR}/${SID}_tumor"

        echo "   HLA-${LOCUS^^}: running at mc${MIN_COV} ..."
        N_RUN=$((N_RUN+1))
        T0=$(date +%s)

        Rscript "${LOHHLA_HOME}/LOHHLAscript.R" \
            --patientId "${SID}" \
            --outputDir "${OUTDIR}" \
            --normalBAMfile "${SID}_normal.bam" \
            --tumorBAMfile "${SID}_tumor.bam" \
            --hlaPath "${DIR}/al_${LOCUS}.txt" \
            --HLAfastaLoc "${SUB}" \
            --HLAexonLoc "${HLA_EXON}" \
            --CopyNumLoc "${DIR}/${SID}_solutions.txt" \
            --mappingStep TRUE --fishingStep FALSE --coverageStep TRUE \
            --plottingStep FALSE --cleanUp FALSE \
            --minCoverageFilter "${MIN_COV}" \
            --numMisMatch 1 --ignoreWarnings TRUE \
            --novoDir "" --gatkDir "" \
            > "${DIR}/lohhla_${LOCUS}_mc${MIN_COV}.log" \
            2> "${DIR}/lohhla_${LOCUS}_mc${MIN_COV}.err"
        RC=$?
        T1=$(date +%s)

        PRED=$(ls "${OUTDIR}"/*HLAlossPrediction*.txt 2>/dev/null | head -1)
        if [ -n "${PRED}" ] && [ "$(wc -l < ${PRED})" -gt 1 ]; then
            printf "       %dm%02ds  " $(( (T1-T0)/60 )) $(( (T1-T0)%60 ))
            awk -F'\t' 'NR==1 {for (i=1;i<=NF;i++) h[$i]=i; next}
              { printf "CN %.3f/%.3f  P=%.4g  sites=%s\n",
                $h["HLA_type1copyNum_withBAFBin"],
                $h["HLA_type2copyNum_withBAFBin"],
                $h["PVal"], $h["numMisMatchSitesCov"] }' "${PRED}" 2>/dev/null
            N_OK=$((N_OK+1))
        else
            R=$(grep -m1 "^Error" "${DIR}/lohhla_${LOCUS}_mc${MIN_COV}.err" \
                2>/dev/null | cut -c1-50)
            echo "       no prediction: ${R}"
            printf "%s\t%s\t%s\thla_%s\n" "${COLLECTION}" "mc${MIN_COV}" \
                "${SID}" "${LOCUS}" >> "${PANEL}/step9d_failures.log"
            N_FAIL=$((N_FAIL+1))
        fi
    done
done

echo "======================================================================"
echo " STEP 9d finished $(date '+%F %T')"
echo " ${COLLECTION} at mc${MIN_COV}"
echo " runs ${N_RUN}   succeeded ${N_OK}   failed ${N_FAIL}   skipped ${N_SKIP}"
echo "======================================================================"
