#!/bin/bash
# =====================================================================
# STEP 9b - run LOHHLA across the four collections.
#
# Identical to step 6b in how it invokes LOHHLA: one locus per run, all
# IMGT subtypes of each allele in the reference FASTA, the same patched
# script. Only the input directory changes, so any difference in the
# results comes from the collections rather than from the tool.
#
# Running one locus at a time matters here more than before. LOHHLA
# processes every locus in the allele file within a single invocation and
# aborts the whole run if one of them fails, so a locus whose alleles
# differ at too few covered positions would take the others down with it.
# With four collections that would compound.
#
# Usage:
#   sbatch step9b_lohhla_panel.sh <collection> <n_brca> <n_ov> [off_b] [off_o]
#
#   sbatch step9b_lohhla_panel.sh A_untouched 10 0 0 0
# =====================================================================
#SBATCH --job-name=step9b_lohhla
#SBATCH --partition=cpu-single
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=16:00:00
#SBATCH --output=/gpfs/bwfor/work/ws/fr_os136-immune_escape/simulation/step9_lohhla_panel/logs/step9b_%j.log
#SBATCH --error=/gpfs/bwfor/work/ws/fr_os136-immune_escape/simulation/step9_lohhla_panel/logs/step9b_%j.err
# =====================================================================

COLLECTION="${1:?give a collection: A_untouched B_loss20 C_loss35 D_loss50}"
N_BRCA="${2:-20}"
N_OV="${3:-20}"
OFF_BRCA="${4:-0}"
OFF_OV="${5:-0}"

WS=/gpfs/bwfor/work/ws/fr_os136-immune_escape
PANEL=${WS}/simulation/step9_lohhla_panel
COL_DIR=${PANEL}/${COLLECTION}
MANIFEST=${WS}/simulation/cohort/manifest.tsv
LOHHLA_HOME=/home/fr/fr_fr/fr_os136/immune_escape_project/soft/lohhla
HLA_FASTA_ALL=${LOHHLA_HOME}/data/hla_all_lohhla.fasta
HLA_EXON=${LOHHLA_HOME}/data/hla.dat
CONDA=/home/fr/fr_fr/fr_os136/miniconda3/etc/profile.d/conda.sh
# Positions are used where both BAMs cover them at least
# this deeply. The LOHHLA default of 5 leaves too few
# positions on panel data; at 3 a test sample gave 8 usable
# positions instead of 1 and the loss was detected.
MIN_COV=3

source "${CONDA}"
conda activate lohhla_env
# samtools 1.23 from bio_work: the patched script needs view -N, and
# lohhla_env is pinned to 1.6 by its R dependencies
export PATH=/home/fr/fr_fr/fr_os136/miniconda3/envs/bio_work/bin:$PATH

echo "======================================================================"
echo " STEP 9b - LOHHLA on ${COLLECTION}"
echo " batch: ${N_BRCA} BRCA (off ${OFF_BRCA}), ${N_OV} OV (off ${OFF_OV})"
echo " started $(date '+%F %T') on $(hostname)"
echo "======================================================================"

[ -d "${COL_DIR}" ] || { echo "ERROR: ${COL_DIR} does not exist"; exit 1; }
for f in "${LOHHLA_HOME}/LOHHLAscript.R" "${HLA_FASTA_ALL}" "${HLA_EXON}"; do
    [ -s "$f" ] || { echo "ERROR: missing $f"; exit 1; }
done

SAMPLES=$(awk -F'\t' -v nb="${N_BRCA}" -v no="${N_OV}" \
                     -v ob="${OFF_BRCA}" -v oo="${OFF_OV}" '
    NR>1 && $2=="brca" { b++; if (b > ob && b <= ob+nb) print $1 }
    NR>1 && $2=="ov"   { o++; if (o > oo && o <= oo+no) print $1 }
' "${MANIFEST}")

N_RUN=0; N_OK=0; N_FAIL=0; N_SKIP=0; N_WAIT=0

for SID in ${SAMPLES}; do
    DIR=${COL_DIR}/${SID}
    NORMAL=${SID}_normal.bam
    TUMOR=${SID}_tumor.bam
    ALLELES=${DIR}/${SID}_alleles.txt
    SOLUTIONS=${DIR}/${SID}_solutions.txt

    echo "----------------------------------------------------------------------"
    echo " ${SID}"

    if [ ! -s "${DIR}/${TUMOR}" ] || [ ! -s "${ALLELES}" ]; then
        echo "   not built yet - run step 9a first"
        N_WAIT=$((N_WAIT+1)); continue
    fi

    cd "${DIR}"

    for LOCUS in a b c; do
        grep "^hla_${LOCUS}_" "${ALLELES}" | sort -u > "${DIR}/al_${LOCUS}.txt"
        N_AL=$(wc -l < "${DIR}/al_${LOCUS}.txt")
        if [ "${N_AL}" -lt 2 ]; then
            rm -f "${DIR}/al_${LOCUS}.txt"
            continue
        fi

        OUTDIR=${DIR}/out_${LOCUS}
        PRED=$(ls "${OUTDIR}"/*HLAlossPrediction*.txt 2>/dev/null | head -1)
        if [ -n "${PRED}" ] && [ "$(wc -l < ${PRED})" -gt 1 ]; then
            echo "   HLA-${LOCUS^^}: already done"
            N_SKIP=$((N_SKIP+1)); continue
        fi

        # every IMGT subtype of these two alleles: reducing the reference to
        # one sequence per allele was tried in step 6 and cut the usable
        # mismatch positions from 95 to 3
        SUB=${DIR}/sub_${LOCUS}.fasta
        rm -f "${SUB}"*
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
print("       FASTA: %d sequences" % n)
PYEOF

        if [ "$(grep -c '^>' ${SUB})" -lt 2 ]; then
            echo "   HLA-${LOCUS^^}: fewer than two sequences matched"
            N_FAIL=$((N_FAIL+1)); continue
        fi

        rm -rf "${OUTDIR}"
        mkdir -p "${OUTDIR}/${NORMAL%.bam}" "${OUTDIR}/${TUMOR%.bam}"

        echo "   HLA-${LOCUS^^}: $(date '+%T') running ..."
        N_RUN=$((N_RUN+1))
        T0=$(date +%s)

        Rscript "${LOHHLA_HOME}/LOHHLAscript.R" \
            --patientId "${SID}" \
            --outputDir "${OUTDIR}" \
            --normalBAMfile "${NORMAL}" \
            --tumorBAMfile "${TUMOR}" \
            --hlaPath "${DIR}/al_${LOCUS}.txt" \
            --HLAfastaLoc "${SUB}" \
            --HLAexonLoc "${HLA_EXON}" \
            --CopyNumLoc "${SOLUTIONS}" \
            --mappingStep TRUE \
            --fishingStep FALSE \
            --coverageStep TRUE \
            --plottingStep TRUE \
            --cleanUp FALSE \
            --minCoverageFilter "${MIN_COV}" \
            --numMisMatch 1 \
            --ignoreWarnings TRUE \
            --novoDir "" --gatkDir "" \
            > "${DIR}/lohhla_${LOCUS}.log" 2> "${DIR}/lohhla_${LOCUS}.err"
        RC=$?
        T1=$(date +%s)

        PRED=$(ls "${OUTDIR}"/*HLAlossPrediction*.txt 2>/dev/null | head -1)
        if [ -n "${PRED}" ] && [ "$(wc -l < ${PRED})" -gt 1 ]; then
            printf "       done in %dm%02ds\n" $(( (T1-T0)/60 )) $(( (T1-T0)%60 ))
            awk -F'\t' 'NR==1 {for (i=1;i<=NF;i++) h[$i]=i; next}
              { printf "       %s vs %s   CN %.3f / %.3f   P=%.4g   sites=%s\n",
                $h["HLA_A_type1"], $h["HLA_A_type2"],
                $h["HLA_type1copyNum_withBAFBin"], $h["HLA_type2copyNum_withBAFBin"],
                $h["PVal"], $h["numMisMatchSitesCov"] }' "${PRED}" 2>/dev/null
            N_OK=$((N_OK+1))
        else
            REASON=$(grep -m1 "^Error" "${DIR}/lohhla_${LOCUS}.err" 2>/dev/null | cut -c1-60)
            echo "       no prediction (exit ${RC}): ${REASON}"
            printf "%s\t%s\thla_%s\trc%s\n" "${COLLECTION}" "${SID}" \
                "${LOCUS}" "${RC}" >> "${PANEL}/step9b_failures.log"
            N_FAIL=$((N_FAIL+1))
        fi
    done
done

echo "======================================================================"
echo " STEP 9b finished $(date '+%F %T')  -  ${COLLECTION}"
echo " locus runs ${N_RUN}   succeeded ${N_OK}   failed ${N_FAIL}"
echo " skipped ${N_SKIP}   waiting ${N_WAIT}"
echo "======================================================================"
