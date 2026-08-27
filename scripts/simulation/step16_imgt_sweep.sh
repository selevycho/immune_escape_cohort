#!/bin/bash
# =====================================================================
# STEP 16 - how many IMGT subtypes the reference needs.
#
# LOHHLA aligns reads against the two alleles it is given and counts which
# reads prefer which. What it is given is a FASTA, and the question is how
# many sequences that FASTA should hold.
#
# An allele call like A*02:01 corresponds to dozens of IMGT sequences that
# differ outside the typed exons - A*02:01:01:01, A*02:01:01:02L and so
# on. During the earlier debugging the reference was reduced to one
# sequence per allele on the reasoning that the rest were redundant, and
# detection collapsed: a sandbox sample went from 95 usable mismatch
# positions to 3. The full set was restored and the matter left there.
#
# This measures it properly. The same sample and the same BAMs are run
# against references holding 1, 5, 20 and all matching subtypes, and the
# number of discriminating positions and the resulting p-value are
# recorded at each. If the relationship is as steep as the sandbox
# suggested, it is worth stating as a finding rather than a footnote:
# anyone reducing an IMGT reference for speed would lose most of their
# detection power without any error to warn them.
#
# Subtypes are taken in IMGT order rather than at random, which is what
# someone trimming a reference by hand would do.
#
# Usage:
#   sbatch step16_imgt_sweep.sh <collection> <sample> [locus] [counts]
#
#   sbatch step16_imgt_sweep.sh C_loss35 O006 b "1 5 20 0"
#   (0 means all available)
# =====================================================================
#SBATCH --job-name=step16_imgt
#SBATCH --partition=cpu-single
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=/gpfs/bwfor/work/ws/fr_os136-immune_escape/simulation/step16_imgt_sweep/logs/step16_%j.log
#SBATCH --error=/gpfs/bwfor/work/ws/fr_os136-immune_escape/simulation/step16_imgt_sweep/logs/step16_%j.err
# =====================================================================

COLLECTION="${1:?collection, e.g. C_loss35}"
SID="${2:?sample, e.g. O006}"
LOCUS="${3:-}"
COUNTS="${4:-1 5 20 0}"

WS=/gpfs/bwfor/work/ws/fr_os136-immune_escape
PANEL=${WS}/simulation/step9_lohhla_panel
SWEEP=${WS}/simulation/step16_imgt_sweep
SRC=${PANEL}/${COLLECTION}/${SID}
LOHHLA_HOME=/home/fr/fr_fr/fr_os136/immune_escape_project/soft/lohhla
HLA_FASTA_ALL=${LOHHLA_HOME}/data/hla_all_lohhla.fasta
HLA_EXON=${LOHHLA_HOME}/data/hla.dat
DESIGN=${PANEL}/collections.tsv
CONDA=/home/fr/fr_fr/fr_os136/miniconda3/etc/profile.d/conda.sh
MIN_COV=5

source "${CONDA}"
conda activate lohhla_env
export PATH=/home/fr/fr_fr/fr_os136/miniconda3/envs/bio_work/bin:$PATH

mkdir -p "${SWEEP}/logs" "${SWEEP}/results"
RESULTS=${SWEEP}/results/imgt_sweep.tsv
[ -s "${RESULTS}" ] || printf "collection\tsample\tlocus\tn_requested\tn_sequences\tallele1\tallele2\tsites\tcn1\tcn2\tpval\tdetected\truntime_s\n" > "${RESULTS}"

echo "======================================================================"
echo " STEP 16 - IMGT subtype count"
echo " ${COLLECTION} / ${SID}   counts: ${COUNTS}"
echo " started $(date '+%F %T') on $(hostname)"
echo "======================================================================"

[ -d "${SRC}" ] || { echo "ERROR: ${SRC} missing"; exit 1; }

# without a locus argument, use the one that was thinned - it is the only
# one where a detection is expected and a change in it means something
if [ -z "${LOCUS}" ]; then
    TGT=$(awk -F'\t' -v c="${COLLECTION}" -v s="${SID}" \
          'NR>1 && $1==c && $2==s {print $4}' "${DESIGN}")
    case "${TGT}" in
        HLA-A) LOCUS=a ;;
        HLA-B) LOCUS=b ;;
        HLA-C) LOCUS=c ;;
        *) echo "ERROR: no target locus recorded for ${SID} in ${COLLECTION}"
           echo "       pass one explicitly as the third argument"
           exit 1 ;;
    esac
    echo " target locus from the design: HLA-${LOCUS^^}"
fi

ALLELES=${SRC}/al_${LOCUS}.txt
if [ ! -s "${ALLELES}" ]; then
    grep "^hla_${LOCUS}_" "${SRC}/${SID}_alleles.txt" | sort -u > "${ALLELES}"
fi
N_AL=$(wc -l < "${ALLELES}")
[ "${N_AL}" -ge 2 ] || { echo "ERROR: locus is homozygous or untyped"; exit 1; }

A1=$(sed -n 1p "${ALLELES}")
A2=$(sed -n 2p "${ALLELES}")
echo " alleles: ${A1} / ${A2}"

TOTAL=$(grep -c "^>\(${A1}\|${A2}\)" "${HLA_FASTA_ALL}" 2>/dev/null || echo 0)
echo " subtypes available in the IMGT reference: ${TOTAL}"
echo

for N in ${COUNTS}; do
    LABEL=$([ "${N}" = "0" ] && echo "all" || echo "${N}")
    OUTDIR=${SWEEP}/${SID}_${LOCUS}_n${LABEL}
    SUB=${SWEEP}/${SID}_${LOCUS}_n${LABEL}.fasta

    echo "----------------------------------------------------------------------"
    echo " ${LABEL} subtype(s) per allele"

    # take the first N sequences of each allele in the order IMGT lists
    # them, which is what trimming a reference by hand produces
    python3 - "${ALLELES}" "${HLA_FASTA_ALL}" "${SUB}" "${N}" << 'PYEOF'
import sys
from collections import defaultdict

allele_file, fasta, out, n = sys.argv[1:5]
n = int(n)
keep = [l.strip() for l in open(allele_file) if l.strip()]

seqs = defaultdict(list)
name, buf = None, []
for line in open(fasta):
    if line.startswith(">"):
        if name:
            for k in keep:
                if name.startswith(k):
                    seqs[k].append((name, "".join(buf)))
                    break
        name, buf = line[1:].strip(), []
    else:
        buf.append(line)
if name:
    for k in keep:
        if name.startswith(k):
            seqs[k].append((name, "".join(buf)))
            break

written = 0
with open(out, "w") as fh:
    for k in keep:
        chosen = seqs[k] if n == 0 else seqs[k][:n]
        for nm, sq in chosen:
            fh.write(">%s\n%s" % (nm, sq))
            written += 1
        print("       %s: %d of %d subtypes" % (k, len(chosen), len(seqs[k])))
print("       %d sequences written" % written)
PYEOF

    N_SEQ=$(grep -c '^>' "${SUB}" 2>/dev/null || echo 0)
    if [ "${N_SEQ}" -lt 2 ]; then
        echo "       fewer than two sequences - cannot run"
        continue
    fi

    rm -rf "${OUTDIR}"
    mkdir -p "${OUTDIR}/${SID}_normal" "${OUTDIR}/${SID}_tumor"

    cd "${SRC}"
    T0=$(date +%s)

    Rscript "${LOHHLA_HOME}/LOHHLAscript.R" \
        --patientId "${SID}" \
        --outputDir "${OUTDIR}" \
        --normalBAMfile "${SID}_normal.bam" \
        --tumorBAMfile "${SID}_tumor.bam" \
        --hlaPath "${ALLELES}" \
        --HLAfastaLoc "${SUB}" \
        --HLAexonLoc "${HLA_EXON}" \
        --CopyNumLoc "${SRC}/${SID}_solutions.txt" \
        --mappingStep TRUE --fishingStep FALSE --coverageStep TRUE \
        --plottingStep FALSE --cleanUp FALSE \
        --minCoverageFilter "${MIN_COV}" \
        --numMisMatch 1 --ignoreWarnings TRUE \
        --novoDir "" --gatkDir "" \
        > "${SWEEP}/logs/${SID}_${LOCUS}_n${LABEL}.log" \
        2> "${SWEEP}/logs/${SID}_${LOCUS}_n${LABEL}.err"
    RC=$?
    T1=$(date +%s)

    PRED=$(ls "${OUTDIR}"/*HLAlossPrediction*.txt 2>/dev/null | head -1)
    if [ -n "${PRED}" ] && [ "$(wc -l < ${PRED})" -gt 1 ]; then
        read -r SITES CN1 CN2 PV < <(awk -F'\t' '
            NR==1 {for (i=1;i<=NF;i++) h[$i]=i; next}
            NR==2 {print $h["numMisMatchSitesCov"],
                         $h["HLA_type1copyNum_withBAFBin"],
                         $h["HLA_type2copyNum_withBAFBin"],
                         $h["PVal"]}' "${PRED}")
        DET=$(awk -v p="${PV}" 'BEGIN {print (p != "NA" && p+0 < 0.05) ? "yes" : "no"}')
        printf "       sites %-5s CN %.3f/%.3f  P=%-12s %s   %ds\n" \
               "${SITES}" "${CN1}" "${CN2}" "${PV}" "${DET}" $((T1-T0))
        printf "%s\t%s\thla_%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "${COLLECTION}" "${SID}" "${LOCUS}" "${LABEL}" "${N_SEQ}" \
            "${A1}" "${A2}" "${SITES}" "${CN1}" "${CN2}" "${PV}" \
            "${DET}" $((T1-T0)) >> "${RESULTS}"
    else
        R=$(grep -m1 "^Error" "${SWEEP}/logs/${SID}_${LOCUS}_n${LABEL}.err" \
            2>/dev/null | cut -c1-50)
        echo "       no prediction (exit ${RC}): ${R}"
        printf "%s\t%s\thla_%s\t%s\t%s\t%s\t%s\tNA\tNA\tNA\tNA\tfailed\t%d\n" \
            "${COLLECTION}" "${SID}" "${LOCUS}" "${LABEL}" "${N_SEQ}" \
            "${A1}" "${A2}" $((T1-T0)) >> "${RESULTS}"
    fi
done

echo "======================================================================"
echo " STEP 16 finished $(date '+%F %T')"
echo "======================================================================"
echo
echo "results so far:"
awk -F'\t' 'NR==1 {printf "  %-8s %-6s %-8s %-6s %-8s %s\n",
                   "sample","locus","subtypes","seqs","sites","P"; next}
  {printf "  %-8s %-6s %-8s %-6s %-8s %s\n", $2, $3, $4, $5, $8, $11}' \
  "${RESULTS}"
