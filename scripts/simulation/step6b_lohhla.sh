#!/bin/bash
# =====================================================================
# STEP 6b - run the patched LOHHLA, one locus per invocation.
#
# LOHHLA processes every locus in the allele file within a single run and
# aborts the whole run if any one of them fails - a locus whose two
# alleles differ at only a handful of covered positions gives its t-test
# constant data, and R halts, discarding the loci that had worked. Running
# each locus separately means a failure costs that locus only.
#
# The allele FASTA keeps every IMGT subtype of each called allele, not one
# representative. Reducing it to a single sequence per allele was tried
# and made things worse: the reference run that succeeded used 834
# sequences for six alleles, and cutting that to six dropped the number of
# usable mismatch positions from 95 to 3.
#
# LOHHLA is a modified copy: novoalign replaced by bwa mem -a, GATK 3 jar
# calls by samtools equivalents, hg19 HLA coordinates by hg38 ones with
# chr prefixes, invalid multi-letter optparse short flags removed, and a
# name-sort inserted before FASTQ conversion so read pairs survive. The
# statistics are untouched.
#
# It also insists on bare BAM filenames rather than paths, and does not
# create the per-BAM working directories it expects - both handled below.
#
# Usage:
#   sbatch step6b_lohhla.sh <n_brca> <n_ov> [off_brca] [off_ov]
# =====================================================================
#SBATCH --job-name=step6b_lohhla
#SBATCH --partition=cpu-single
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/bwfor/work/ws/fr_os136-immune_escape/logs/step6b_%j.log
#SBATCH --error=/gpfs/bwfor/work/ws/fr_os136-immune_escape/logs/step6b_%j.err
# =====================================================================

N_BRCA="${1:-20}"
N_OV="${2:-20}"
OFF_BRCA="${3:-0}"
OFF_OV="${4:-0}"

WS=/gpfs/bwfor/work/ws/fr_os136-immune_escape
MANIFEST=${WS}/simulation/cohort/manifest.tsv
LOH_DIR=${WS}/simulation/lohhla_allelic
DESIGN=${LOH_DIR}/loh_design.tsv
LOHHLA_HOME=/home/fr/fr_fr/fr_os136/immune_escape_project/soft/lohhla
HLA_FASTA_ALL=${LOHHLA_HOME}/data/hla_all_lohhla.fasta
HLA_EXON=${LOHHLA_HOME}/data/hla.dat
CONDA=/home/fr/fr_fr/fr_os136/miniconda3/etc/profile.d/conda.sh

MIN_COV=5

source "${CONDA}"
conda activate lohhla_env
# samtools 1.23 from bio_work: the patched script needs "view -N" (1.12+)
# and lohhla_env is pinned to 1.6 by its R dependencies
export PATH=/home/fr/fr_fr/fr_os136/miniconda3/envs/bio_work/bin:$PATH

echo "======================================================================"
echo " STEP 6b - LOHHLA, one locus per run"
echo " batch: ${N_BRCA} BRCA (off ${OFF_BRCA}), ${N_OV} OV (off ${OFF_OV})"
echo " started $(date '+%F %T') on $(hostname)"
echo "======================================================================"

for f in "${LOHHLA_HOME}/LOHHLAscript.R" "${HLA_FASTA_ALL}" "${HLA_EXON}"; do
    [ -s "$f" ] || { echo "ERROR: missing $f"; exit 1; }
done

SAMPLES=$(awk -F'\t' -v nb="${N_BRCA}" -v no="${N_OV}" \
                     -v ob="${OFF_BRCA}" -v oo="${OFF_OV}" '
    NR>1 && $2=="brca" { b++; if (b > ob && b <= ob+nb) print $1 }
    NR>1 && $2=="ov"   { o++; if (o > oo && o <= oo+no) print $1 }
' "${MANIFEST}")

echo
echo "samples: $(echo ${SAMPLES} | tr '\n' ' ')"
echo

N_RUN=0; N_OK=0; N_FAIL=0; N_WAIT=0

for SID in ${SAMPLES}; do
    DIR=${LOH_DIR}/${SID}
    NORMAL=${SID}_normal.bam
    TUMOR=${SID}_tumor_LOH.bam
    ALLELES=${DIR}/${SID}_alleles.txt
    SOLUTIONS=${DIR}/${SID}_solutions.txt

    echo "----------------------------------------------------------------------"
    echo " ${SID}"

    if [ ! -s "${DIR}/${TUMOR}" ] || [ ! -s "${ALLELES}" ]; then
        echo "   LOH BAM or allele file missing - run step 6a first"
        N_WAIT=$((N_WAIT+1)); continue
    fi

    TARGETS=$(awk -F'\t' -v s="${SID}" 'NR>1 && $1==s {print $3}' "${DESIGN}" 2>/dev/null | tr '\n' ',' | sed 's/,$//')
    echo "   loci with simulated loss: ${TARGETS:-unknown}"

    cd "${DIR}"

    for LOCUS in a b c; do
        # only heterozygous loci are worth running; LOHHLA refuses the rest
        grep "^hla_${LOCUS}_" "${ALLELES}" | sort -u > "${DIR}/al_${LOCUS}.txt"
        N_AL=$(wc -l < "${DIR}/al_${LOCUS}.txt")
        if [ "${N_AL}" -lt 2 ]; then
            echo "   HLA-${LOCUS^^}: homozygous or untyped - skipping"
            rm -f "${DIR}/al_${LOCUS}.txt"
            continue
        fi

        OUTDIR=${DIR}/out_${LOCUS}
        PRED=$(ls "${OUTDIR}"/*HLAlossPrediction*.txt 2>/dev/null | head -1)
        if [ -n "${PRED}" ] && [ "$(wc -l < ${PRED})" -gt 1 ]; then
            echo "   HLA-${LOCUS^^}: already done"
            N_OK=$((N_OK+1)); continue
        fi

        # every IMGT subtype of these two alleles - see the note above
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
            echo "   HLA-${LOCUS^^}: fewer than two sequences matched - skipping"
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
            echo "       no prediction (exit ${RC}): $(tail -2 ${DIR}/lohhla_${LOCUS}.err | head -1)"
            printf "%s\t%s\trc%s\n" "${SID}" "hla_${LOCUS}" "${RC}" \
                >> "${LOH_DIR}/step6b_failures.log"
            N_FAIL=$((N_FAIL+1))
        fi
    done
done

echo "======================================================================"
echo " STEP 6b finished $(date '+%F %T')"
echo " locus runs ${N_RUN}   succeeded ${N_OK}   failed ${N_FAIL}   waiting ${N_WAIT}"
echo "======================================================================"

echo
echo "all predictions:"
printf " %-6s %-6s %-22s %-22s %8s %8s %10s %6s %s\n" \
       "sample" "locus" "allele 1" "allele 2" "CN1" "CN2" "PVal" "sites" "target"
for d in ${LOH_DIR}/*/out_*/; do
    P=$(ls "${d}"*HLAlossPrediction*.txt 2>/dev/null | head -1)
    [ -n "${P}" ] && [ "$(wc -l < ${P})" -gt 1 ] || continue
    s=$(echo "$d" | sed "s|${LOH_DIR}/||; s|/out_.*||")
    loc=$(echo "$d" | sed 's|.*/out_||; s|/||')
    tgt=$(awk -F'\t' -v x="$s" -v l="HLA-${loc^^}" \
          'NR>1 && $1==x && $3==l {print "yes"}' "${DESIGN}" 2>/dev/null | head -1)
    awk -F'\t' -v s="$s" -v loc="HLA-${loc^^}" -v tgt="${tgt:-control}" '
      NR==1 {for (i=1;i<=NF;i++) h[$i]=i; next}
      { printf " %-6s %-6s %-22s %-22s %8.3f %8.3f %10.4g %6s %s\n",
        s, loc, $h["HLA_A_type1"], $h["HLA_A_type2"],
        $h["HLA_type1copyNum_withBAFBin"], $h["HLA_type2copyNum_withBAFBin"],
        $h["PVal"], $h["numMisMatchSitesCov"],
        (tgt == "yes" ? "TARGET" : "control") }' "${P}" 2>/dev/null
done
