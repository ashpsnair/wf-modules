#!/bin/bash

#PBS -l select=1:ncpus=8:mem=16g
#PBS -l walltime=02:00:00
#PBS -P 11003581
#PBS -N base-filter
#PBS -j oe
#PBS -m ae

# =============================================================================
# base-filter.sh
#
# PURPOSE: Apply quality-based hard filters to the common (intersected) VCFs
#          produced by commonVC_variants.sh.
#
# FILTERS APPLIED (per variant):
#   - FILTER == PASS                (already guaranteed by commonVC, kept as
#                                    a safety check on isec output)
#   - Tumor depth  (FORMAT/DP[1])  >= 20
#   - Normal depth (FORMAT/DP[0])  >= 10
#   - Tumor ALT reads (FORMAT/AD[1:1]) >= 5
#   - Tumor VAF    (AD[1:1]/DP[1]) >= 0.10  (10 %)
#
# INPUT  : $INTERSECTION_DIR/<sample>/0002.vcf.gz
#            (common Mutect2 records from commonVC_variants.sh)
#
# OUTPUT : $BASE_FILTER_DIR/<sample>.base_filtered.vcf.gz
#          $BASE_DIR/base_filter_summary.csv
#
# NOTES  : FORMAT field indexing uses bcftools 0-based sample order:
#            [0] = NORMAL sample (first column after FORMAT in Mutect2 VCF)
#            [1] = TUMOR  sample (second column)
#          Verify the order in your VCFs with:
#            bcftools view -h <sample>.vcf.gz | grep "^#CHROM"
#          Adjust [0]/[1] indices if your column order differs.
# =============================================================================

set -euo pipefail

cd "$PBS_O_WORKDIR"

echo "=== base-filter STARTED: $(date) ==="

module load bcftools/1.15.1

# ---------------------------------------------------------------------------
# PATHS
# ---------------------------------------------------------------------------

BASE_DIR="/home/users/nus/ash.ps/scratch/YS-organoid-downstream/analysis2.1"

INTERSECTION_DIR="${BASE_DIR}/01_commonVC"
BASE_FILTER_DIR="${BASE_DIR}/02_base_filtered"

mkdir -p "$BASE_FILTER_DIR"

# ---------------------------------------------------------------------------
# THRESHOLDS
# ---------------------------------------------------------------------------

MIN_TUMOR_DP=20
MIN_NORMAL_DP=10
MIN_ALT_READS=5
MIN_TUMOR_VAF=0.10

# ---------------------------------------------------------------------------
# SUMMARY CSV
# ---------------------------------------------------------------------------

SUMMARY_CSV="${BASE_DIR}/base_filter_summary.csv"
echo "sample,common_input,base_filtered,removed" > "$SUMMARY_CSV"

# ---------------------------------------------------------------------------
# MAIN LOOP
# ---------------------------------------------------------------------------

echo ""
echo "Thresholds:"
echo "  Tumor DP  >= $MIN_TUMOR_DP"
echo "  Normal DP >= $MIN_NORMAL_DP"
echo "  ALT reads >= $MIN_ALT_READS"
echo "  Tumor VAF >= $MIN_TUMOR_VAF"

for isec_dir in "${INTERSECTION_DIR}"/*; do

    [[ -d "$isec_dir" ]] || continue

    sample=$(basename "$isec_dir")
    input_vcf="${isec_dir}/0002.vcf.gz"
    out_vcf="${BASE_FILTER_DIR}/${sample}.base_filtered.vcf.gz"

    echo ""
    echo "--- $sample ---"

    if [[ ! -f "$input_vcf" ]]; then
        echo "[WARN] Missing input: $input_vcf — skipping"
        continue
    fi

    n_before=$(bcftools view -H "$input_vcf" | wc -l)

    # -----------------------------------------------------------------------
    # Apply filters
    # FORMAT/DP[1]    = tumor total depth
    # FORMAT/DP[0]    = normal total depth
    # FORMAT/AD[1:1]  = tumor ALT allele read count (second allele, tumor sample)
    # -----------------------------------------------------------------------

    bcftools filter \
        --include \
        "FORMAT/DP[1] >= ${MIN_TUMOR_DP} &&
         FORMAT/DP[0] >= ${MIN_NORMAL_DP} &&
         FORMAT/AD[1:1] >= ${MIN_ALT_READS} &&
         (FORMAT/AD[1:1] / FORMAT/DP[1]) >= ${MIN_TUMOR_VAF}" \
        -O z \
        -o "$out_vcf" \
        "$input_vcf"

    bcftools index -f --tbi "$out_vcf"

    n_after=$(bcftools view -H "$out_vcf" | wc -l)
    n_removed=$((n_before - n_after))

    echo "  Input : $n_before"
    echo "  Output: $n_after (removed $n_removed)"

    echo "$sample,$n_before,$n_after,$n_removed" >> "$SUMMARY_CSV"

done

# ---------------------------------------------------------------------------
# DONE
# ---------------------------------------------------------------------------

echo ""
echo "=== base-filter COMPLETE: $(date) ==="
echo ""
echo "Filtered VCFs : $BASE_FILTER_DIR"
echo "Summary CSV   : $SUMMARY_CSV"
echo ""
cat "$SUMMARY_CSV"