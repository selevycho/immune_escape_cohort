#!/bin/bash
# =====================================================================
# STEP 8 - add indels to the simulated tumours.
#
# Steps 1-7 injected substitutions only. Indels are 15.2% of the panel
# mutations in BRCA and 12.5% in OV, so leaving them out understates what
# a caller has to find. This step takes the same donors' indels from the
# lifted MAF and injects them into a copy of each tumour BAM.
#
# The SNV tumours are not modified. Everything here goes to a parallel
# directory so the SNV-only results stay reproducible and the two can be
# compared directly.
#
# What was learned testing this on single samples:
#
#   BAMSurgeon injects indels reliably - 39 of 40 deletions and 37 of 37
#   insertions landed at 33x coverage with VAFs from 8%. An earlier
#   report of 2% success was an artefact of the check, not the injection:
#   mpileup marks a deletion in the column BEFORE the deleted base, as
#   -1G, and shows the base itself as *. Looking only at the deletion's
#   own coordinate finds nothing.
#
#   addindel.py emits an unsorted BAM. Indexing fails unless a sort is
#   inserted between the two.
#
#   Indels inside the MHC are lost, as SNVs are: BAMSurgeon realigns
#   through bwa and the hg38 reference carries alternate MHC haplotypes,
#   so reads return with zero mapping quality. Those positions are
#   excluded up front rather than counted as failures.
#
# Usage:
#   sbatch step8_inject_indels.sh <n_brca> <n_ov> [off_brca] [off_ov]
# =====================================================================
#SBATCH --job-name=step8_indel
#SBATCH --partition=cpu-single
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/bwfor/work/ws/fr_os136-immune_escape/logs/step8_%j.log
#SBATCH --error=/gpfs/bwfor/work/ws/fr_os136-immune_escape/logs/step8_%j.err
# =====================================================================

N_BRCA="${1:-20}"
N_OV="${2:-20}"
OFF_BRCA="${3:-0}"
OFF_OV="${4:-0}"

WS=/gpfs/bwfor/work/ws/fr_os136-immune_escape
MANIFEST=${WS}/simulation/cohort/manifest.tsv
COHORT=${WS}/simulation/cohort
INDEL=${WS}/simulation/indels
PANEL=${WS}/simulation/panel/panel.bed
REF=${WS}/ref/Homo_sapiens_assembly38.fasta
CONDA=/home/fr/fr_fr/fr_os136/miniconda3/etc/profile.d/conda.sh
THREADS="${SLURM_CPUS_PER_TASK:-2}"

MIN_VAF=0.05
MIN_INDELS=1

source "${CONDA}"
mkdir -p "${INDEL}"

echo "======================================================================"
echo " STEP 8 - indel injection"
echo " batch: ${N_BRCA} BRCA (off ${OFF_BRCA}), ${N_OV} OV (off ${OFF_OV})"
echo " started $(date '+%F %T') on $(hostname)"
echo "======================================================================"

TRUTH_ALL=${INDEL}/indel_truth_all.tsv
[ -s "${TRUTH_ALL}" ] || printf "sample\tcohort\tchrom\tpos\tend\tref\talt\ttype\tgene\tVAF\tconsequence\n" > "${TRUTH_ALL}"

SAMPLES=$(awk -F'\t' -v nb="${N_BRCA}" -v no="${N_OV}" \
                     -v ob="${OFF_BRCA}" -v oo="${OFF_OV}" '
    NR>1 && $2=="brca" { b++; if (b > ob && b <= ob+nb) print $1 }
    NR>1 && $2=="ov"   { o++; if (o > oo && o <= oo+no) print $1 }
' "${MANIFEST}")

echo
echo "samples: $(echo ${SAMPLES} | tr '\n' ' ')"
echo

OK=0; SKIP=0; FAIL=0; WAIT=0; NONE=0

for SID in ${SAMPLES}; do
    LINE=$(awk -F'\t' -v s="${SID}" 'NR>1 && $1==s' "${MANIFEST}")
    COH=$(echo "${LINE}"     | cut -f2)
    BARCODE=$(echo "${LINE}" | cut -f3)
    MAF=${WS}/liftover/out/${COH}.hg38.maf.tsv

    SRC=${COHORT}/${SID}/${SID}_tumor.bam
    OUT=${INDEL}/${SID}
    FINAL=${OUT}/${SID}_tumor_snv_indel.bam

    echo "----------------------------------------------------------------------"
    printf " %-6s %-5s %s\n" "${SID}" "${COH}" "${BARCODE}"

    if [ ! -s "${SRC}" ]; then
        echo "   SNV tumour BAM not ready"
        WAIT=$((WAIT+1)); continue
    fi
    if [ -s "${FINAL}.bai" ]; then
        echo "   already built - skipping"
        SKIP=$((SKIP+1)); continue
    fi

    mkdir -p "${OUT}"
    T0=$(date +%s)

    # ---- this donor's indels, inside the panel, outside the MHC ----
    conda activate cptac_env
    python3 - "${MAF}" "${PANEL}" "${BARCODE}" "${SID}" "${COH}" \
              "${OUT}" "${MIN_VAF}" << 'PYEOF'
import sys
import numpy as np
import pandas as pd

maf_p, panel_p, barcode, sid, cohort, out, min_vaf = sys.argv[1:8]
min_vaf = float(min_vaf)

MHC = ("chr6", 29600000, 33100000)
INDEL_CLASS = {"Frame_Shift_Del", "Frame_Shift_Ins",
               "In_Frame_Del", "In_Frame_Ins"}

panel = pd.read_csv(panel_p, sep="\t", header=None,
                    names=["chrom", "start", "end", "name"])
by = {c: g[["start", "end"]].to_numpy() for c, g in panel.groupby("chrom")}

m = pd.read_csv(maf_p, sep="\t", low_memory=False)
m = m[m.Tumor_Sample_Barcode == barcode]
m = m[m.Variant_Classification.isin(INDEL_CLASS)]
n_donor = len(m)

den = m.t_ref_count.fillna(0) + m.t_alt_count.fillna(0)
m = m[den > 0].copy()
m["VAF"] = (m.t_alt_count / den).round(4)
m = m[m.VAF >= min_vaf]

keep = []
for c, p in zip(m.Chromosome_hg38.values, m.Start_Position_hg38.values):
    iv = by.get(c)
    keep.append(False if iv is None
                else bool(((iv[:, 0] <= int(p) - 1) & (int(p) - 1 < iv[:, 1])).any()))
m = m[keep]
n_panel = len(m)

in_mhc = ((m.Chromosome_hg38 == MHC[0]) &
          (m.Start_Position_hg38 >= MHC[1]) &
          (m.Start_Position_hg38 < MHC[2]))
n_mhc = int(in_mhc.sum())
m = m[~in_mhc]

# BAMSurgeon format: chrom, start, end, VAF, type[, sequence]
rows, truth = [], []
for _, r in m.iterrows():
    c = r.Chromosome_hg38
    s = int(r.Start_Position_hg38)
    e = int(r.End_Position_hg38)
    vaf = float(r.VAF)
    ref = str(r.Reference_Allele)
    alt = str(r.Tumor_Seq_Allele2)
    if r.Variant_Type == "DEL" or alt == "-":
        rows.append("%s\t%d\t%d\t%.4f\tDEL" % (c, s, max(e, s + 1), vaf))
        kind = "DEL"
    else:
        seq = alt if alt not in ("-", "") else "A"
        rows.append("%s\t%d\t%d\t%.4f\tINS\t%s" % (c, s, s, vaf, seq))
        kind = "INS"
    truth.append("%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%.4f\t%s" %
                 (sid, cohort, c, s, e, ref, alt, kind,
                  r.Hugo_Symbol, vaf, r.Variant_Classification))

open("%s/indels.txt" % out, "w").write("\n".join(rows) + ("\n" if rows else ""))
open("%s/truth_rows.tsv" % out, "w").write("\n".join(truth) + ("\n" if truth else ""))

print("   donor indels: %d, in panel: %d, dropped from MHC: %d, to inject: %d"
      % (n_donor, n_panel, n_mhc, len(rows)))
if rows:
    v = m.VAF
    print("   VAF %.3f - %.3f (median %.3f)   types: %s"
          % (v.min(), v.max(), v.median(), m.Variant_Type.value_counts().to_dict()))
PYEOF

    N_IND=$(wc -l < "${OUT}/indels.txt" 2>/dev/null || echo 0)
    if [ "${N_IND}" -lt "${MIN_INDELS}" ]; then
        echo "   only ${N_IND} indels - below the ${MIN_INDELS} floor, skipping"
        printf "%s\ttoo_few_%s\n" "${SID}" "${N_IND}" >> "${INDEL}/skipped.log"
        NONE=$((NONE+1)); continue
    fi

    # ---- copy the SNV tumour; the cohort file is never modified ----
    WORK=${OUT}/${SID}_input.bam
    if [ ! -s "${WORK}.bai" ]; then
        conda activate bio_work
        echo "   copying the SNV tumour BAM ..."
        cp "${SRC}" "${WORK}"
        cp "${SRC}.bai" "${WORK}.bai" 2>/dev/null \
            || samtools index -@ "${THREADS}" "${WORK}"
    fi

    # ---- inject ----
    conda activate bamsurgeon_env
    PICARD=$(find "${CONDA_PREFIX}" -name "picard.jar" 2>/dev/null | head -1)
    if [ -z "${PICARD}" ]; then
        echo "   picard.jar not found"
        FAIL=$((FAIL+1)); continue
    fi

    cd "${OUT}"
    rm -rf addindel.tmp
    echo "   $(date '+%T') addindel.py on ${N_IND} indels ..."

    if ! addindel.py \
            -v "${OUT}/indels.txt" \
            -f "${WORK}" -r "${REF}" \
            -o "${OUT}/raw.bam" \
            -p "${THREADS}" --aligner mem \
            --picardjar "${PICARD}" --force --insane \
            > "${OUT}/addindel.log" 2>&1; then
        echo "   FAILED in addindel.py:"
        tail -6 "${OUT}/addindel.log" | sed 's/^/     /'
        printf "%s\taddindel\n" "${SID}" >> "${INDEL}/failures.log"
        rm -rf "${OUT}/raw.bam"* "${OUT}/addindel.tmp"
        FAIL=$((FAIL+1)); continue
    fi

    # addindel never emits a coordinate-sorted BAM
    conda activate bio_work
    samtools sort -@ "${THREADS}" -m 1G -o "${FINAL}" "${OUT}/raw.bam" 2>/dev/null
    samtools index -@ "${THREADS}" "${FINAL}"
    rm -rf "${OUT}/raw.bam"* "${OUT}/addindel.tmp" "${OUT}"/addindel_logs_*

    if [ ! -s "${FINAL}" ]; then
        echo "   FAILED sorting"
        printf "%s\tsort\n" "${SID}" >> "${INDEL}/failures.log"
        FAIL=$((FAIL+1)); continue
    fi

    T1=$(date +%s)
    printf "   injected in %dm%02ds  |  %s\n" \
           $(( (T1-T0)/60 )) $(( (T1-T0)%60 )) \
           "$(ls -lh ${FINAL} | awk '{print $5}')"

    # ---- verify, reading the pileup the way indels are actually written ----
    HIT=0
    while IFS=$'\t' read -r chrom start end vaf typ seq; do
        if [ "${typ}" = "DEL" ]; then
            win="${chrom}:$((start-1))-$((start+1))"
        else
            win="${chrom}:${start}-${start}"
        fi
        n=$(samtools mpileup -B -q 0 -Q 0 -r "${win}" -f "${REF}" "${FINAL}" 2>/dev/null \
            | awk '{print $5}' | grep -o '[+-][0-9]\+\|\*' | wc -l)
        [ "${n}" -gt 0 ] && HIT=$((HIT+1))
    done < "${OUT}/indels.txt"

    printf "   verified: %d / %d landed (%d%%)\n" \
           "${HIT}" "${N_IND}" $(( 100 * HIT / N_IND ))

    if [ "${HIT}" -lt $(( N_IND / 2 )) ]; then
        printf "%s\tlow_yield_%s_of_%s\n" "${SID}" "${HIT}" "${N_IND}" \
            >> "${INDEL}/failures.log"
    fi

    cat "${OUT}/truth_rows.tsv" >> "${TRUTH_ALL}"
    rm -f "${OUT}/truth_rows.tsv" "${WORK}" "${WORK}.bai"

    OK=$((OK+1))
done

echo "======================================================================"
echo " STEP 8 finished $(date '+%F %T')"
echo " injected ${OK}   skipped ${SKIP}   too few indels ${NONE}   waiting ${WAIT}   failed ${FAIL}"
echo "======================================================================"

echo
echo "cohort:"
printf " %-6s %-5s %8s %10s %12s\n" "sample" "coh" "indels" "BAM" "types"
TOT=0; N=0
while IFS=$'\t' read -r sid cohort rest; do
    [ "${sid}" = "sample_id" ] && continue
    f=${INDEL}/${sid}/indels.txt
    b=${INDEL}/${sid}/${sid}_tumor_snv_indel.bam
    [ -s "$f" ] || continue
    n=$(wc -l < "$f")
    types=$(awk -F'\t' '{print $5}' "$f" | sort | uniq -c | tr '\n' ' ' | tr -s ' ')
    size=$([ -s "$b" ] && ls -lh "$b" | awk '{print $5}' || echo "-")
    printf " %-6s %-5s %8s %10s %12s\n" "${sid}" "${cohort}" "${n}" "${size}" "${types}"
    TOT=$((TOT+n)); N=$((N+1))
done < "${MANIFEST}"
echo
echo "  samples with indels : ${N}"
echo "  indels total        : ${TOT}"
