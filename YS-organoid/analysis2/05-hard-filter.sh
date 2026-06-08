#!/bin/bash

#PBS -l select=1:ncpus=8:mem=16g
#PBS -l walltime=02:00:00
#PBS -P 11003581
#PBS -N hard-filter-1
#PBS -j oe
#PBS -m ae

# =============================================================================
# hard-filter-1.sh
#
# PURPOSE: Apply stringent hard filters to base-filtered VCFs.
#          This filter branch is independent of cellLRK and is designed
#          to produce a high-confidence somatic variant set.
#
# FILTERS APPLIED:
#   1. Tumor VAF >= 20%                      (from FORMAT fields in VCF)
#   2. ExAC_ALL   <= 0.01  or  "."           (from ANNOVAR annotation in INFO)
#   3. gnomAD_genome_ALL <= 0.01 or "."      (from ANNOVAR annotation in INFO)
#
# INPUT  : $ANNOTATED_DIR/<sample>/<sample>.hg38_multianno.vcf.gz
#            (ANNOVAR-annotated base-filtered VCF)
#
# OUTPUT : $HARD1_DIR/<sample>.hard1_filtered.vcf.gz
#          $BASE_DIR/hard_filter1_summary.csv
#
# NOTES  : Population AF values are embedded in the INFO field by ANNOVAR's
#          -vcfinput mode. The awk filter reads the INFO string and parses
#          ExAC_ALL and gnomAD_genome_ALL tags. Missing/dot values are treated
#          as 0 (i.e. rare / not in database → kept).
#          FORMAT-based VAF filter uses bcftools first, then population AF
#          uses awk on the INFO field in a second pass.
# =============================================================================

set -euo pipefail

cd "$PBS_O_WORKDIR"

echo "=== hard-filter-1 STARTED: $(date) ==="

module load bcftools/1.15.1

# ---------------------------------------------------------------------------
# PATHS
# ---------------------------------------------------------------------------

BASE_DIR="/home/users/nus/ash.ps/scratch/YS-organoid-downstream/analysis2.1"

ANNOTATED_DIR="${BASE_DIR}/03_annotated"
HARD1_DIR="${BASE_DIR}/05_hard_filter1"

TEMP_DIR="/tmp/hard1_$$"

mkdir -p "$HARD1_DIR" "$TEMP_DIR"

# ---------------------------------------------------------------------------
# THRESHOLDS
# ---------------------------------------------------------------------------

MIN_TUMOR_VAF=0.20
MAX_POP_AF=0.01

# ---------------------------------------------------------------------------
# SUMMARY CSV
# ---------------------------------------------------------------------------

SUMMARY_CSV="${BASE_DIR}/hard_filter1_summary.csv"
echo "sample,input,after_vaf20,after_popAF,removed_total" > "$SUMMARY_CSV"

# ---------------------------------------------------------------------------
# MAIN LOOP
# ---------------------------------------------------------------------------

echo ""
echo "Thresholds:"
echo "  Tumor VAF     >= $MIN_TUMOR_VAF (20%)"
echo "  Max pop AF    <= $MAX_POP_AF    (ExAC + gnomAD)"

for sample_dir in "${ANNOTATED_DIR}"/*; do

    [[ -d "$sample_dir" ]] || continue

    sample=$(basename "$sample_dir")
    input_vcf="${sample_dir}/${sample}.hg38_multianno.vcf.gz"
    out_vcf="${HARD1_DIR}/${sample}.hard1_filtered.vcf.gz"

    echo ""
    echo "--- $sample ---"

    if [[ ! -f "$input_vcf" ]]; then
        echo "[WARN] Missing annotated VCF: $input_vcf — skipping"
        continue
    fi

    n_input=$(bcftools view -H "$input_vcf" | wc -l)
    echo "  Input: $n_input"

    # -----------------------------------------------------------------------
    # PASS 1: Tumor VAF >= 20% (FORMAT-based filter via bcftools)
    # FORMAT/AD[1:1] = tumor ALT read count
    # FORMAT/DP[1]   = tumor total depth
    # -----------------------------------------------------------------------

    bcftools filter \
        --include "(FORMAT/AD[1:1] / FORMAT/DP[1]) >= ${MIN_TUMOR_VAF}" \
        -O z \
        -o "$TEMP_DIR/${sample}.vaf20.vcf.gz" \
        "$input_vcf"

    bcftools index -f --tbi "$TEMP_DIR/${sample}.vaf20.vcf.gz"

    n_after_vaf=$(bcftools view -H "$TEMP_DIR/${sample}.vaf20.vcf.gz" | wc -l)
    echo "  After VAF>=20%: $n_after_vaf"

    # -----------------------------------------------------------------------
    # PASS 2: Population AF filter (INFO field from ANNOVAR annotation)
    # Parses ExAC_ALL and gnomAD_genome_ALL tags from the INFO column.
    # Variants with AF > MAX_POP_AF in either database are excluded.
    # Missing values (.) are treated as AF=0 (rare, kept).
    # -----------------------------------------------------------------------

    bcftools view -h "$TEMP_DIR/${sample}.vaf20.vcf.gz" \
        > "$TEMP_DIR/${sample}.header.txt"

    bcftools view -H "$TEMP_DIR/${sample}.vaf20.vcf.gz" | \
    awk -v max_af="$MAX_POP_AF" '
    BEGIN { OFS="\t" }
    {
        split($8, info_fields, ";")

        exac_af   = 0
        gnomad_af = 0

        for (i in info_fields) {
            split(info_fields[i], kv, "=")
            if (kv[1] == "ExAC_ALL" && kv[2] != "." && kv[2] != "")
                exac_af = kv[2] + 0
            if (kv[1] == "gnomAD_genome_ALL" && kv[2] != "." && kv[2] != "")
                gnomad_af = kv[2] + 0
        }

        if (exac_af > max_af || gnomad_af > max_af)
            next

        print
    }
    ' > "$TEMP_DIR/${sample}.popAF_pass.txt"

    # Reconstruct valid VCF from header + filtered records
    cat "$TEMP_DIR/${sample}.header.txt" \
        "$TEMP_DIR/${sample}.popAF_pass.txt" | \
        bcftools view -O z -o "$out_vcf"

    bcftools index -f --tbi "$out_vcf"

    n_after_popAF=$(bcftools view -H "$out_vcf" | wc -l)
    n_removed=$((n_input - n_after_popAF))

    echo "  After pop AF filter: $n_after_popAF"
    echo "  Total removed: $n_removed"

    echo "$sample,$n_input,$n_after_vaf,$n_after_popAF,$n_removed" \
        >> "$SUMMARY_CSV"

    # -----------------------------------------------------------------------
    # Clean temp
    # -----------------------------------------------------------------------

    rm -f \
        "$TEMP_DIR/${sample}.vaf20.vcf.gz" \
        "$TEMP_DIR/${sample}.vaf20.vcf.gz.tbi" \
        "$TEMP_DIR/${sample}.header.txt" \
        "$TEMP_DIR/${sample}.popAF_pass.txt"

    echo "[OK] Done: $sample"

done

# ---------------------------------------------------------------------------
# CLEANUP
# ---------------------------------------------------------------------------

rm -rf "$TEMP_DIR"

# ---------------------------------------------------------------------------
# DONE
# ---------------------------------------------------------------------------

echo ""
echo "=== hard-filter-1 COMPLETE: $(date) ==="
echo ""
echo "Filtered VCFs : $HARD1_DIR"
echo "Summary CSV   : $SUMMARY_CSV"
echo ""
cat "$SUMMARY_CSV"