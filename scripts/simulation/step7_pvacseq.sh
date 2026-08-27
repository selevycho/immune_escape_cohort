#!/bin/bash
# =====================================================================
# STEP 7 - neoantigen prediction with pVACseq and NetMHCpan.
#
# This is the second, independent route to the same answer. The mhcflurry
# path in step 5 took mutations from the lifted MAF with HGVSp already
# annotated by TCGA and cut peptides straight out of GENCODE protein
# sequences. Here the Mutect2 VCF goes through VEP, which recomputes the
# protein consequence from the genome, and pVACseq builds the peptides
# itself using the Wildtype and Frameshift plugins before handing them to
# NetMHCpan-4.1.
#
# Running both routes on the same input is the point. They share no code
# and start from different files, so agreement between them is evidence
# that neither carries a systematic error; disagreement localises where
# one of them does.
#
# Three details that cost time to find:
#
#   Sample order in the VCF. Mutect2 writes the normal in column 10 and
#   the tumour in column 11. Passing them the wrong way round produces a
#   valid-looking run with an empty result, since the "tumour" then has no
#   somatic variants at all.
#
#   Duplicate HLA alleles. A homozygous locus appears twice in the
#   OptiType output; pVACseq wants each allele once.
#
#   VEP needs its plugins on an explicit --dir_plugins path, and VEP and
#   pVACseq live in separate conda environments because VEP's perl stack
#   and pvactools' tensorflow cannot be resolved together.
#
# Usage:
#   sbatch step7_pvacseq.sh <n_brca> <n_ov> [off_brca] [off_ov] [algorithms]
#
#   sbatch step7_pvacseq.sh 20 20                      NetMHCpan alone
#   sbatch step7_pvacseq.sh 10 0 0 0 "NetMHCpan MHCflurry"   both engines
# =====================================================================
#SBATCH --job-name=step7_pvac
#SBATCH --partition=cpu-single
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/bwfor/work/ws/fr_os136-immune_escape/logs/step7_%j.log
#SBATCH --error=/gpfs/bwfor/work/ws/fr_os136-immune_escape/logs/step7_%j.err
# =====================================================================

N_BRCA="${1:-20}"
N_OV="${2:-20}"
OFF_BRCA="${3:-0}"
OFF_OV="${4:-0}"
ALGORITHMS="${5:-NetMHCpan}"

WS=/gpfs/bwfor/work/ws/fr_os136-immune_escape
MANIFEST=${WS}/simulation/cohort/manifest.tsv
COHORT=${WS}/simulation/cohort
PVAC=${WS}/simulation/pvacseq/cohort
REF=${WS}/ref/Homo_sapiens_assembly38.fasta
VEP_CACHE=${WS}/ref/vep_cache
VEP_PLUGINS=/home/fr/fr_fr/fr_os136/vep_plugins
NETMHC=${WS}/soft/netmhcpan/netMHCpan-4.1
CONDA=/home/fr/fr_fr/fr_os136/miniconda3/etc/profile.d/conda.sh
CPUS="${SLURM_CPUS_PER_TASK:-4}"

EPITOPE_LENGTHS=8,9,10,11

source "${CONDA}"
mkdir -p "${PVAC}"

echo "======================================================================"
echo " STEP 7 - pVACseq with ${ALGORITHMS}"
echo " batch: ${N_BRCA} BRCA (off ${OFF_BRCA}), ${N_OV} OV (off ${OFF_OV})"
echo " started $(date '+%F %T') on $(hostname)"
echo "======================================================================"

[ -d "${VEP_CACHE}/homo_sapiens" ] || { echo "ERROR: VEP cache missing"; exit 1; }
[ -s "${NETMHC}/netMHCpan" ] || { echo "ERROR: netMHCpan missing"; exit 1; }
[ -s "${VEP_PLUGINS}/Wildtype.pm" ] || { echo "ERROR: VEP plugins missing"; exit 1; }

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
    VCF=${COHORT}/${SID}/mutect2/${SID}.filtered.vcf.gz
    HLA_TSV=${COHORT}/${SID}/optitype/${SID}_result.tsv
    OUT=${PVAC}/${SID}
    ANNOT=${OUT}/${SID}.vep.vcf
    RESULT_DIR=${OUT}/pvacseq_out

    echo "----------------------------------------------------------------------"
    echo " ${SID}"

    if [ ! -s "${VCF}" ] || [ ! -s "${HLA_TSV}" ]; then
        echo "   VCF or HLA type not ready"
        WAIT=$((WAIT+1)); continue
    fi

    FINAL=$(ls "${RESULT_DIR}"/MHC_Class_I/*.filtered.tsv 2>/dev/null | head -1)
    if [ -n "${FINAL}" ] && [ "$(wc -l < ${FINAL})" -gt 1 ]; then
        echo "   already done: $(( $(wc -l < ${FINAL}) - 1 )) epitopes"
        SKIP=$((SKIP+1)); continue
    fi

    mkdir -p "${OUT}"
    T0=$(date +%s)

    # ---- HLA alleles, unique, in pVACseq notation ----
    ALLELES=$(awk 'NR==2 {
        n = 0
        for (i = 2; i <= 7; i++)
            if ($i ~ /\*/ && !seen[$i]++) {
                if (n++) printf ","
                printf "HLA-%s", $i
            }
    }' "${HLA_TSV}")
    echo "   HLA: ${ALLELES}"

    if [ -z "${ALLELES}" ]; then
        echo "   no usable HLA alleles"
        printf "%s\tno_alleles\n" "${SID}" >> "${PVAC}/failures.log"
        FAIL=$((FAIL+1)); continue
    fi

    # ---- 1. VEP ----
    if [ ! -s "${ANNOT}" ]; then
        conda activate vep_env
        echo "   [1/2] $(date '+%T') VEP ..."
        vep \
            --input_file "${VCF}" \
            --output_file "${ANNOT}" \
            --format vcf --vcf \
            --symbol --terms SO --tsl --biotype --hgvs \
            --fasta "${REF}" \
            --offline --cache --dir_cache "${VEP_CACHE}" \
            --plugin Frameshift --plugin Wildtype \
            --dir_plugins "${VEP_PLUGINS}" \
            --pick --transcript_version \
            --fork "${CPUS}" --force_overwrite \
            > "${OUT}/vep.log" 2>&1
        RC=$?
        if [ "${RC}" -ne 0 ] || [ ! -s "${ANNOT}" ]; then
            echo "   VEP FAILED (exit ${RC}):"
            tail -6 "${OUT}/vep.log" | sed 's/^/     /'
            printf "%s\tvep_rc%s\n" "${SID}" "${RC}" >> "${PVAC}/failures.log"
            FAIL=$((FAIL+1)); continue
        fi
        echo "         $(grep -vc '^#' ${ANNOT}) records annotated"
    else
        echo "   [1/2] VEP output already present"
    fi

    # ---- 2. pVACseq ----
    # Mutect2 column order: 10 = normal, 11 = tumour
    NORMAL_NAME=$(grep '^#CHROM' "${ANNOT}" | cut -f10)
    TUMOR_NAME=$(grep '^#CHROM' "${ANNOT}" | cut -f11)
    echo "   [2/2] $(date '+%T') pVACseq: tumour ${TUMOR_NAME}, normal ${NORMAL_NAME}"

    conda activate pvac_env
    export PATH=${NETMHC}:$PATH

    rm -rf "${RESULT_DIR}"
    pvacseq run \
        "${ANNOT}" \
        "${TUMOR_NAME}" \
        "${ALLELES}" \
        ${ALGORITHMS} \
        "${RESULT_DIR}" \
        -e1 "${EPITOPE_LENGTHS}" \
        --normal-sample-name "${NORMAL_NAME}" \
        --pass-only \
        --normal-cov 5 --tdna-cov 5 \
        --normal-vaf 0.02 --tdna-vaf 0.05 \
        --maximum-transcript-support-level 5 \
        --n-threads "${CPUS}" \
        > "${OUT}/pvacseq.log" 2>&1
    RC=$?
    T1=$(date +%s)

    ALL=$(ls "${RESULT_DIR}"/MHC_Class_I/*.all_epitopes.tsv 2>/dev/null | head -1)
    FILT=$(ls "${RESULT_DIR}"/MHC_Class_I/*.filtered.tsv 2>/dev/null | head -1)

    if [ -z "${ALL}" ] || [ "$(wc -l < ${ALL})" -le 1 ]; then
        echo "   no epitopes produced (exit ${RC})"
        grep -iE "error|exception|empty" "${OUT}/pvacseq.log" | tail -4 | sed 's/^/     /'
        printf "%s\tno_epitopes_rc%s\n" "${SID}" "${RC}" >> "${PVAC}/failures.log"
        FAIL=$((FAIL+1)); continue
    fi

    N_ALL=$(( $(wc -l < ${ALL}) - 1 ))
    N_FILT=0
    [ -n "${FILT}" ] && N_FILT=$(( $(wc -l < ${FILT}) - 1 ))

    printf "         %dm%02ds   %d epitopes, %d after filtering\n" \
           $(( (T1-T0)/60 )) $(( (T1-T0)%60 )) "${N_ALL}" "${N_FILT}"

    # Strong and weak binders are counted from all_epitopes with the same
    # percentile thresholds used for mhcflurry in step 5, so the two routes
    # are directly comparable.
    #
    # pVACseq's own filtered.tsv is not used for this. Its filter chain -
    # binding, coverage, transcript support, then one best epitope per
    # mutation - is designed to shortlist vaccine candidates and reduced
    # 49690 epitopes to three on the first sample. That is the right answer
    # to a different question.
    if [ -n "${ALL}" ]; then
        awk -F'\t' 'NR==1 {
            for (i=1;i<=NF;i++) if ($i == "Best MT Percentile") col=i
            next
        }
        col && $col != "NA" {
            if ($col+0 < 0.5) sb++
            else if ($col+0 < 2) wb++
        }
        END { printf "         strong (<0.5%%): %d   weak (0.5-2%%): %d\n", sb+0, wb+0 }' "${ALL}"
    fi

    AGG=$(ls "${RESULT_DIR}"/MHC_Class_I/*.aggregated.tsv 2>/dev/null | head -1)
    if [ -n "${AGG}" ]; then
        echo "         mutations with an epitope: $(( $(wc -l < ${AGG}) - 1 ))"
    fi

    OK=$((OK+1))
done

END_ALL=$(date +%s)
echo "======================================================================"
echo " STEP 7 finished $(date '+%F %T')"
printf " elapsed %dh%02dm   done %d   skipped %d   waiting %d   failed %d\n" \
       $(( (END_ALL-START_ALL)/3600 )) $(( ((END_ALL-START_ALL)%3600)/60 )) \
       "${OK}" "${SKIP}" "${WAIT}" "${FAIL}"
echo "======================================================================"

echo
echo "cohort summary:"
printf " %-6s %-5s %10s %10s %10s\n" "sample" "coh" "epitopes" "filtered" "strong"
TOT_E=0; TOT_F=0; N=0
while IFS=$'\t' read -r sid cohort rest; do
    [ "${sid}" = "sample_id" ] && continue
    A=$(ls ${PVAC}/${sid}/pvacseq_out/MHC_Class_I/*.all_epitopes.tsv 2>/dev/null | head -1)
    [ -n "${A}" ] || continue
    F=$(ls ${PVAC}/${sid}/pvacseq_out/MHC_Class_I/*.filtered.tsv 2>/dev/null | head -1)
    ne=$(( $(wc -l < ${A}) - 1 ))
    nf=0; [ -n "${F}" ] && nf=$(( $(wc -l < ${F}) - 1 ))
    ns=$(awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++) if ($i == "Best MT Percentile") c=i; next}
                      c && $c != "NA" && $c+0 < 0.5 {n++} END {print n+0}' "${A}")
    printf " %-6s %-5s %10s %10s %10s\n" "${sid}" "${cohort}" "${ne}" "${nf}" "${ns}"
    TOT_E=$((TOT_E+ne)); TOT_F=$((TOT_F+nf)); N=$((N+1))
done < "${MANIFEST}"

echo
echo "  samples with predictions : ${N}"
echo "  epitopes total           : ${TOT_E}"
echo "  after filtering          : ${TOT_F}"
