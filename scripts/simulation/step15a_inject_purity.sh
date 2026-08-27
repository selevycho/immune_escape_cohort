#!/bin/bash
# =====================================================================
# STEP 15a - inject the same mutations at fixed allele fractions.
#
# Step 3 already gives recall against VAF, but the fractions there come
# from the donors and arrive tangled with everything else: one mutation
# sits at 0.08 in a well-covered gene, another at 0.45 in a poorly covered
# one, and the bins mix genes, positions and coverage together.
#
# Here the same mutations are re-injected into the same backbone at four
# fixed fractions. Position, gene, coverage and read set are held constant
# and only the fraction moves, so the curve measures what the fraction
# does rather than what the donors happened to carry.
#
# The fraction is what a pathologist reports as tumour purity, modulated
# by whether every clone carries the mutation. Fixing it turns the
# question into one an experimenter can act on: at what purity does this
# panel stop finding mutations.
#
# Only substitutions are injected. addsnv.py cannot place an indel, and
# mixing the two would confound the curve with the indel failure rate.
#
# Usage:
#   sbatch step15a_inject_purity.sh <n_brca> <n_ov> [off_b] [off_o] [vafs]
#
#   sbatch step15a_inject_purity.sh 1 0
#   sbatch step15a_inject_purity.sh 5 5 0 0 "0.05 0.10 0.20 0.40"
# =====================================================================
#SBATCH --job-name=step15a_vaf
#SBATCH --partition=cpu-single
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=20:00:00
#SBATCH --output=/gpfs/bwfor/work/ws/fr_os136-immune_escape/simulation/step15_purity_sweep/logs/step15a_%j.log
#SBATCH --error=/gpfs/bwfor/work/ws/fr_os136-immune_escape/simulation/step15_purity_sweep/logs/step15a_%j.err
# =====================================================================

N_BRCA="${1:-5}"
N_OV="${2:-5}"
OFF_BRCA="${3:-0}"
OFF_OV="${4:-0}"
VAFS="${5:-0.05 0.10 0.20 0.40}"

WS=/gpfs/bwfor/work/ws/fr_os136-immune_escape
COHORT=${WS}/simulation/cohort
SWEEP=${WS}/simulation/step15_purity_sweep
MANIFEST=${COHORT}/manifest.tsv
REF=${WS}/ref/Homo_sapiens_assembly38.fasta
CONDA=/home/fr/fr_fr/fr_os136/miniconda3/etc/profile.d/conda.sh
THREADS="${SLURM_CPUS_PER_TASK:-2}"

MHC_START=29600000
MHC_END=33100000

source "${CONDA}"

mkdir -p "${SWEEP}/logs" "${SWEEP}/results"
DESIGN=${SWEEP}/levels.tsv
[ -s "${DESIGN}" ] || printf "level\tsample\tcohort\tvaf\tmutations\tinjected\tlanded\tmedian_observed_vaf\n" > "${DESIGN}"

echo "======================================================================"
echo " STEP 15a - re-inject at fixed allele fractions"
echo " fractions: ${VAFS}"
echo " batch: ${N_BRCA} BRCA (off ${OFF_BRCA}), ${N_OV} OV (off ${OFF_OV})"
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

OK=0; SKIP=0; FAIL=0

for SID in ${SAMPLES}; do
    COH=$(awk -F'\t' -v s="${SID}" 'NR>1 && $1==s {print $2}' "${MANIFEST}")
    SRC_N=${COHORT}/${SID}/${SID}_normal.bam
    TRUTH=${COHORT}/${SID}/truth_set.tsv

    echo "----------------------------------------------------------------------"
    echo " ${SID}  (${COH})"

    if [ ! -s "${SRC_N}" ] || [ ! -s "${TRUTH}" ]; then
        echo "   normal BAM or truth set missing"
        FAIL=$((FAIL+1)); continue
    fi

    # substitutions only, outside the MHC - injection there fails
    # systematically because bwa realigns through alternate haplotypes
    N_SNV=$(awk -F'\t' -v s=${MHC_START} -v e=${MHC_END} \
        'NR>1 && $7=="SNP" && !($2=="chr6" && $3>=s && $3<e)' "${TRUTH}" | wc -l)
    echo "   ${N_SNV} substitutions outside the MHC"

    if [ "${N_SNV}" -lt 5 ]; then
        echo "   too few to measure a rate - skipping"
        printf "%s\ttoo_few_%s\n" "${SID}" "${N_SNV}" >> "${SWEEP}/skipped.log"
        SKIP=$((SKIP+1)); continue
    fi

    for VAF in ${VAFS}; do
        TAG=vaf$(awk -v v="${VAF}" 'BEGIN {printf "%02d", v*100}')
        OUT=${SWEEP}/${TAG}/${SID}

        if [ -s "${OUT}/${SID}_tumor.bam.bai" ]; then
            echo "   ${TAG}: already built"
            SKIP=$((SKIP+1)); continue
        fi

        mkdir -p "${OUT}"
        T0=$(date +%s)

        # the same positions as the original truth set, with the fraction
        # replaced; the file doubles as the answer key for step 15b
        conda activate cptac_env
        python3 - "${TRUTH}" "${VAF}" "${OUT}" "${MHC_START}" "${MHC_END}" << 'PYEOF'
import sys
import pandas as pd

truth_p, vaf, out, mhc_s, mhc_e = sys.argv[1:6]
vaf = float(vaf)
mhc_s, mhc_e = int(mhc_s), int(mhc_e)

t = pd.read_csv(truth_p, sep="\t")
if "Variant_Type" in t.columns:
    t = t[t.Variant_Type == "SNP"]
t = t[~((t.Chromosome_hg38 == "chr6") &
        (t.Start_Position_hg38 >= mhc_s) &
        (t.Start_Position_hg38 < mhc_e))]

# BAMSurgeon wants chrom, start, end, VAF, alt
rows = ["%s\t%d\t%d\t%.4f\t%s" % (r.Chromosome_hg38,
                                  int(r.Start_Position_hg38),
                                  int(r.Start_Position_hg38),
                                  vaf, r.Tumor_Seq_Allele2)
        for _, r in t.iterrows()]
open("%s/snvs.txt" % out, "w").write("\n".join(rows) + "\n")

t = t.copy()
t["VAF"] = vaf
t.to_csv("%s/truth_set.tsv" % out, sep="\t", index=False)
print("   %d positions at VAF %.2f" % (len(rows), vaf))
PYEOF

        N_IN=$(wc -l < "${OUT}/snvs.txt")

        # the normal is the input: injecting into an already-mutated
        # tumour would stack this run's mutations on the earlier ones
        conda activate bio_work
        cp "${SRC_N}" "${OUT}/input.bam"
        cp "${SRC_N}.bai" "${OUT}/input.bam.bai" 2>/dev/null \
            || samtools index -@ "${THREADS}" "${OUT}/input.bam"
        cp "${SRC_N}" "${OUT}/${SID}_normal.bam"
        cp "${SRC_N}.bai" "${OUT}/${SID}_normal.bam.bai" 2>/dev/null \
            || samtools index -@ "${THREADS}" "${OUT}/${SID}_normal.bam"

        conda activate bamsurgeon_env
        PICARD=$(find "${CONDA_PREFIX}" -name "picard.jar" 2>/dev/null | head -1)
        cd "${OUT}"
        rm -rf addsnv.tmp

        echo "   ${TAG}: $(date '+%T') injecting ${N_IN} SNVs ..."
        addsnv.py \
            -v "${OUT}/snvs.txt" \
            -f "${OUT}/input.bam" -r "${REF}" \
            -o "${OUT}/raw.bam" \
            -p "${THREADS}" --aligner mem \
            --picardjar "${PICARD}" --force --ignoresnps --tagreads \
            > "${OUT}/addsnv.log" 2>&1
        RC=$?

        conda activate bio_work
        if [ "${RC}" -ne 0 ] || [ ! -s "${OUT}/raw.bam" ]; then
            echo "   ${TAG}: addsnv FAILED (exit ${RC})"
            tail -4 "${OUT}/addsnv.log" | sed 's/^/     /'
            printf "%s\t%s\taddsnv\n" "${SID}" "${TAG}" >> "${SWEEP}/failures.log"
            rm -rf "${OUT}/raw.bam"* "${OUT}/addsnv.tmp" "${OUT}/input.bam"*
            FAIL=$((FAIL+1)); continue
        fi

        samtools sort -@ "${THREADS}" -m 1G \
            -o "${OUT}/${SID}_tumor.bam" "${OUT}/raw.bam" 2>/dev/null
        samtools index -@ "${THREADS}" "${OUT}/${SID}_tumor.bam"
        rm -rf "${OUT}/raw.bam"* "${OUT}/addsnv.tmp" "${OUT}"/addsnv_logs_* \
               "${OUT}/input.bam"*

        # verify rather than assume: BAMSurgeon reports success at
        # positions where realignment later removed the edited reads
        LANDED=0
        SUM_VAF=0
        while IFS=$'\t' read -r C P E V A; do
            OBS=$(samtools mpileup -B -q 0 -Q 0 \
                    --ff UNMAP,SECONDARY,QCFAIL \
                    -r "${C}:${P}-${P}" -f "${REF}" \
                    "${OUT}/${SID}_tumor.bam" 2>/dev/null \
                  | awk -v alt="${A}" '{
                      # gsub needs an assignable target, so the base
                      # string is copied into a variable first; calling it
                      # on toupper($5) directly makes awk abort and every
                      # position then reports zero
                      b = toupper($5); a = toupper(alt)
                      n = gsub(a, "", b)
                      if ($4 > 0) printf "%d %.4f", n, n/$4
                    }')
            set -- ${OBS}
            if [ -n "$1" ] && [ "$1" -gt 0 ] 2>/dev/null; then
                LANDED=$((LANDED+1))
                SUM_VAF=$(awk -v s="${SUM_VAF}" -v v="$2" 'BEGIN {print s+v}')
            fi
        done < "${OUT}/snvs.txt"

        MED=$(awk -v s="${SUM_VAF}" -v n="${LANDED}" \
              'BEGIN {printf "%.4f", (n ? s/n : 0)}')

        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "${TAG}" "${SID}" "${COH}" "${VAF}" "${N_SNV}" "${N_IN}" \
            "${LANDED}" "${MED}" >> "${DESIGN}"

        T1=$(date +%s)
        printf "       %d of %d landed, mean observed VAF %s   %dm%02ds\n" \
               "${LANDED}" "${N_IN}" "${MED}" \
               $(( (T1-T0)/60 )) $(( (T1-T0)%60 ))
        OK=$((OK+1))
    done
done

echo "======================================================================"
echo " STEP 15a finished $(date '+%F %T')"
echo " built ${OK}   skipped ${SKIP}   failed ${FAIL}"
echo "======================================================================"

echo
echo "landing rate by fraction:"
awk -F'\t' 'NR>1 {n[$4]++; inj[$4]+=$6; land[$4]+=$7; obs[$4]+=$8}
  END {for (v in n) printf "  VAF %-6s %4d of %4d landed (%.0f%%), observed %.3f\n",
       v, land[v], inj[v], 100*land[v]/inj[v], obs[v]/n[v]}' "${DESIGN}" | sort
