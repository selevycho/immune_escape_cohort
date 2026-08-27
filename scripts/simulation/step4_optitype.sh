#!/bin/bash
# =====================================================================
# STEP 4 - HLA class I typing with OptiType, batched.
#
# Built directly on run_optitype_S001.sh, which worked: same partition,
# same 8 CPUs, same -v flag, same region windows. The only change is the
# loop over a batch of samples instead of one hard-coded path.
#
# An earlier attempt with 4 CPUs hung indefinitely without writing a
# single byte. razers3 is configured for 16 threads in config.ini, and
# under cgroups it appears to stall rather than fall back. Eight CPUs is
# what the working version used, so that is what this uses.
#
# Typing runs on the NORMAL BAM: the HLA genotype is germline and does
# not change with the injected mutations. The answer also determines
# whether LOHHLA can be simulated at all, since a sample homozygous at
# every class I locus has no allele to lose.
#
# Usage:
#   sbatch step4_optitype.sh <n_brca> <n_ov> [offset_brca] [offset_ov]
# =====================================================================
#SBATCH --job-name=step4_hla
#SBATCH --partition=cpu-single
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=/gpfs/bwfor/work/ws/fr_os136-immune_escape/logs/step4_%j.log
#SBATCH --error=/gpfs/bwfor/work/ws/fr_os136-immune_escape/logs/step4_%j.err
# =====================================================================

N_BRCA="${1:-20}"
N_OV="${2:-20}"
OFF_BRCA="${3:-0}"
OFF_OV="${4:-0}"

WS=/gpfs/bwfor/work/ws/fr_os136-immune_escape
MANIFEST=${WS}/simulation/cohort/manifest.tsv
COHORT_DIR=${WS}/simulation/cohort
OPTITYPE_BIN=/home/fr/fr_fr/fr_os136/miniconda3/envs/optitype_stable/bin/OptiTypePipeline.py
CPUS="${SLURM_CPUS_PER_TASK:-8}"

REGIONS="chr6:29932000-29956000 chr6:31258000-31282000 chr6:31343000-31367000"

source /home/fr/fr_fr/fr_os136/miniconda3/etc/profile.d/conda.sh

echo "======================================================================"
echo " STEP 4 - OptiType HLA class I typing"
echo " batch: ${N_BRCA} BRCA (offset ${OFF_BRCA}), ${N_OV} OV (offset ${OFF_OV})"
echo " started $(date '+%F %T') on $(hostname), ${CPUS} CPUs"
echo "======================================================================"

SAMPLES=$(awk -F'\t' -v nb="${N_BRCA}" -v no="${N_OV}" \
                     -v ob="${OFF_BRCA}" -v oo="${OFF_OV}" '
    NR>1 && $2=="brca" { b++; if (b > ob && b <= ob+nb) print $1 }
    NR>1 && $2=="ov"   { o++; if (o > oo && o <= oo+no) print $1 }
' "${MANIFEST}")

echo
echo "samples: $(echo ${SAMPLES} | tr '\n' ' ')"
echo

OK=0; SKIP=0; FAIL=0; WAIT=0

for SID in ${SAMPLES}; do
    OUT=${COHORT_DIR}/${SID}
    BAM=${OUT}/${SID}_normal.bam
    FQDIR=${OUT}/hla_fastq
    OUTDIR=${OUT}/optitype

    echo "----------------------------------------------------------------------"
    echo " ${SID}"

    if [ ! -s "${BAM}" ]; then
        echo "   normal BAM not ready"
        WAIT=$((WAIT+1)); continue
    fi
    if [ -s "${OUTDIR}/${SID}_result.tsv" ]; then
        echo -n "   already typed: "
        awk 'NR==2 {print $2,$3,$4,$5,$6,$7}' "${OUTDIR}/${SID}_result.tsv"
        SKIP=$((SKIP+1)); continue
    fi

    rm -rf "${OUTDIR}"
    mkdir -p "${FQDIR}" "${OUTDIR}"
    T0=$(date +%s)

    conda activate bio_work

    echo "   [$(date '+%T')] extracting HLA class I reads ..."
    samtools view -@ "${CPUS}" -b "${BAM}" ${REGIONS} \
        > "${FQDIR}/hla_region.bam" 2>/dev/null
    N_READS=$(samtools view -c "${FQDIR}/hla_region.bam")
    echo "   reads extracted: ${N_READS}"

    if [ "${N_READS}" -lt 200 ]; then
        echo "   too few reads - skipping"
        printf "%s\tlow_reads_%s\n" "${SID}" "${N_READS}" \
            >> "${COHORT_DIR}/step4_failures.log"
        rm -f "${FQDIR}/hla_region.bam"
        FAIL=$((FAIL+1)); continue
    fi

    echo "   [$(date '+%T')] converting to FASTQ ..."
    samtools sort -n -@ "${CPUS}" -m 2G \
        -o "${FQDIR}/hla_namesorted.bam" "${FQDIR}/hla_region.bam" 2>/dev/null
    samtools fastq -@ "${CPUS}" \
        -1 "${FQDIR}/${SID}_1.fastq.gz" \
        -2 "${FQDIR}/${SID}_2.fastq.gz" \
        -0 /dev/null -s /dev/null -n \
        "${FQDIR}/hla_namesorted.bam" 2>/dev/null
    rm -f "${FQDIR}/hla_region.bam" "${FQDIR}/hla_namesorted.bam"

    if [ ! -s "${FQDIR}/${SID}_1.fastq.gz" ]; then
        echo "   FASTQ conversion produced nothing"
        printf "%s\tfastq\n" "${SID}" >> "${COHORT_DIR}/step4_failures.log"
        FAIL=$((FAIL+1)); continue
    fi

    echo "   [$(date '+%T')] running OptiType ..."
    conda activate optitype_stable

    RC=0
    python "${OPTITYPE_BIN}" \
        -i "${FQDIR}/${SID}_1.fastq.gz" "${FQDIR}/${SID}_2.fastq.gz" \
        --dna -v \
        -o "${OUTDIR}" \
        --prefix "${SID}" \
        > "${OUT}/optitype.log" 2>&1 || RC=$?

    if [ "${RC}" -ne 0 ]; then
        echo "   OptiType exited with code ${RC}"
        tail -8 "${OUT}/optitype.log" | sed 's/^/     /'
        printf "%s\trc%s\n" "${SID}" "${RC}" >> "${COHORT_DIR}/step4_failures.log"
        FAIL=$((FAIL+1)); continue
    fi

    # OptiType may write into a timestamped subdirectory
    if [ ! -s "${OUTDIR}/${SID}_result.tsv" ]; then
        FRESH=$(find "${OUTDIR}" -name "*result.tsv" -newermt "@${T0}" | head -1)
        [ -n "${FRESH}" ] && cp "${FRESH}" "${OUTDIR}/${SID}_result.tsv"
    fi

    if [ ! -s "${OUTDIR}/${SID}_result.tsv" ]; then
        echo "   no result file produced"
        printf "%s\tno_result\n" "${SID}" >> "${COHORT_DIR}/step4_failures.log"
        FAIL=$((FAIL+1)); continue
    fi

    T1=$(date +%s)
    printf "   done in %dm%02ds\n" $(( (T1-T0)/60 )) $(( (T1-T0)%60 ))
    awk 'NR==2 {printf "   A: %-10s %-10s\n   B: %-10s %-10s\n   C: %-10s %-10s\n   reads used: %s\n",
                $2,$3,$4,$5,$6,$7,$8}' "${OUTDIR}/${SID}_result.tsv"

    rm -f "${FQDIR}"/*.fastq.gz
    OK=$((OK+1))
done

echo "======================================================================"
echo " STEP 4 finished $(date '+%F %T')"
echo " typed ${OK}   skipped ${SKIP}   waiting ${WAIT}   failed ${FAIL}"
echo "======================================================================"

echo
echo "genotypes so far:"
printf " %-6s  %-10s %-10s %-10s %-10s %-10s %-10s  %s\n" \
       "sample" "A1" "A2" "B1" "B2" "C1" "C2" "het"
N=0; HETA=0
for f in ${COHORT_DIR}/*/optitype/*_result.tsv; do
    [ -s "$f" ] || continue
    s=$(basename "$f" _result.tsv)
    read -r a1 a2 b1 b2 c1 c2 < <(awk 'NR==2 {print $2,$3,$4,$5,$6,$7}' "$f")
    het=""
    [ "$a1" != "$a2" ] && { het="${het}A"; HETA=$((HETA+1)); }
    [ "$b1" != "$b2" ] && het="${het}B"
    [ "$c1" != "$c2" ] && het="${het}C"
    [ -z "$het" ] && het="-"
    printf " %-6s  %-10s %-10s %-10s %-10s %-10s %-10s  %s\n" \
           "$s" "$a1" "$a2" "$b1" "$b2" "$c1" "$c2" "$het"
    N=$((N+1))
done
echo
echo "  typed: ${N}   heterozygous at HLA-A: ${HETA}"
