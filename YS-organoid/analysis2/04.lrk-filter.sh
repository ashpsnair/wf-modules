#!/bin/bash

#PBS -l select=1:ncpus=8:mem=16g
#PBS -l walltime=02:00:00
#PBS -P 11003581
#PBS -N cellLRK-filter
#PBS -j oe
#PBS -m ae

# =============================================================================
# cellLRK-filter.sh
#
# PURPOSE: Reproduce the variant filtering strategy described in the Cell
#          paper (LRK organoid study). Applied to ANNOVAR-annotated files.
#
# FILTER CRITERIA (from methods):
#   REMOVE if annotated as Benign or Likely_benign in ClinVar
#   REMOVE if population AF > 0.01 in ExAC (exac03) or gnomAD genome
#   REMOVE variants on chrX or chrY
#
# INPUT  : $ANNOTATED_DIR/<sample>/<sample>.hg38_multianno.txt
#          $ANNOTATED_DIR/<sample>/<sample>.hg38_multianno.vcf.gz
#
# OUTPUT : $CELLLRK_DIR/<sample>.cellLRK_filtered.vcf.gz
#          $CELLLRK_DIR/<sample>.regions.bed  (intermediate; kept for QC)
#          $BASE_DIR/cellLRK_filter_summary.csv
#
# NOTES  : The awk step reads the multianno TXT to apply all annotation-based
#          filters, then generates a BED file of passing regions. bcftools
#          view -R then extracts those variants from the annotated VCF.
#          Column indices for ExAC_ALL and gnomAD depend on the exact
#          ANNOVAR protocol order; the awk header-parser auto-detects them.
# =============================================================================

set -euo pipefail

cd "$PBS_O_WORKDIR"

echo "=== cellLRK-filter STARTED: $(date) ==="

module load bcftools/1.15.1

# ---------------------------------------------------------------------------
# PATHS
# ---------------------------------------------------------------------------

BASE_DIR="/home/users/nus/ash.ps/scratch/YS-organoid-downstream/analysis2.1"

ANNOTATED_DIR="${BASE_DIR}/03_annotated"
CELLLRK_DIR="${BASE_DIR}/04_cellLRK_filtered"

mkdir -p "$CELLLRK_DIR"

SUMMARY_CSV="${BASE_DIR}/cellLRK_filter_summary.csv"
echo "sample,input_total,after_sex_chr,after_clinvar,after_popAF,final" \
    > "$SUMMARY_CSV"

# ---------------------------------------------------------------------------
# POPULATION AF CUTOFF
# ---------------------------------------------------------------------------

MAX_POP_AF=0.01

# ---------------------------------------------------------------------------
# MAIN LOOP
# ---------------------------------------------------------------------------

for sample_dir in "${ANNOTATED_DIR}"/*; do

    [[ -d "$sample_dir" ]] || continue

    sample=$(basename "$sample_dir")
    multianno="${sample_dir}/${sample}.hg38_multianno.txt"
    anno_vcf_gz="${sample_dir}/${sample}.hg38_multianno.vcf.gz"

    echo ""
    echo "=== Processing: $sample ==="

    if [[ ! -f "$multianno" ]]; then
        echo "[WARN] Missing multianno TXT: $multianno — skipping"
        continue
    fi

    if [[ ! -f "$anno_vcf_gz" ]]; then
        echo "[WARN] Missing annotated VCF: $anno_vcf_gz — skipping"
        continue
    fi

    regions_bed="${CELLLRK_DIR}/${sample}.regions.bed"
    out_vcf="${CELLLRK_DIR}/${sample}.cellLRK_filtered.vcf.gz"

    total=$(($(wc -l < "$multianno") - 1))
    echo "  Total annotated: $total"

    # -----------------------------------------------------------------------
    # AWK: Parse multianno header to locate columns, then apply filters
    #      and emit BED regions for passing variants
    # -----------------------------------------------------------------------

    awk -F '\t' '

    BEGIN { OFS="\t" }

    # ---------- Header: auto-detect column indices -------------------------
    NR == 1 {
        for (i = 1; i <= NF; i++) {
            if ($i == "Chr")                        chr_col    = i
            if ($i == "Start")                      start_col  = i
            if ($i == "End")                        end_col    = i
            if ($i ~ /CLNSIG/ || $i ~ /clnsig/)    clinvar_col = i
            if ($i == "ExAC_ALL")                   exac_col   = i
            if ($i ~ /^gnomAD_genome_ALL$/)         gnomad_col = i
            # Fallback: any column matching gnomad + genome + ALL
            if ($i ~ /gnomad/ && $i ~ /genome/ && $i ~ /ALL/ && gnomad_col == "")
                gnomad_col = i
        }

        print "Column mapping:" > "/dev/stderr"
        print "  Chr="    chr_col    > "/dev/stderr"
        print "  Start="  start_col  > "/dev/stderr"
        print "  End="    end_col    > "/dev/stderr"
        print "  CLNSIG=" clinvar_col > "/dev/stderr"
        print "  ExAC="   exac_col   > "/dev/stderr"
        print "  gnomAD=" gnomad_col > "/dev/stderr"

        # Counters
        n_total   = 0
        n_sex     = 0
        n_clinvar = 0
        n_popAF   = 0

        next
    }

    # ---------- Data rows: apply filters -----------------------------------
    {
        n_total++

        # 1. Remove sex chromosomes
        if ($chr_col == "chrX" || $chr_col == "chrY") {
            n_sex++
            next
        }

        # 2. Remove ClinVar benign / likely benign
        if (clinvar_col != "" && $clinvar_col != "." && $clinvar_col != "") {
            if ($clinvar_col ~ /Benign/ || $clinvar_col ~ /Likely_benign/) {
                n_clinvar++
                next
            }
        }

        # 3. Population AF filter (ExAC and gnomAD)
        exac_af   = 0
        gnomad_af = 0

        if (exac_col != "" && $exac_col != "." && $exac_col != "")
            exac_af = $exac_col + 0

        if (gnomad_col != "" && $gnomad_col != "." && $gnomad_col != "")
            gnomad_af = $gnomad_col + 0

        if (exac_af > '"$MAX_POP_AF"' || gnomad_af > '"$MAX_POP_AF"') {
            n_popAF++
            next
        }

        # 4. Emit BED region (0-based start)
        print $chr_col, ($start_col - 1), $end_col

    }

    # ---------- End: print filter statistics to stderr ---------------------
    END {
        print "Filter stats:" > "/dev/stderr"
        print "  Total:            " n_total   > "/dev/stderr"
        print "  Removed sex chr:  " n_sex     > "/dev/stderr"
        print "  Removed ClinVar:  " n_clinvar > "/dev/stderr"
        print "  Removed pop AF:   " n_popAF   > "/dev/stderr"
        print "  Remaining:        " (n_total - n_sex - n_clinvar - n_popAF) > "/dev/stderr"
    }

    ' "$multianno" \
    2> >(tee /dev/stderr | grep "Filter stats" -A 5 > "${CELLLRK_DIR}/${sample}.filter_stats.txt") \
    > "$regions_bed"

    n_regions=$(wc -l < "$regions_bed")
    echo "  BED regions passing filters: $n_regions"

    # -----------------------------------------------------------------------
    # Extract passing variants from the annotated VCF using the BED file
    # -----------------------------------------------------------------------

    if [[ "$n_regions" -eq 0 ]]; then
        echo "[WARN] No variants passed filters for $sample"
        final=0
    else
        bcftools view \
            -R "$regions_bed" \
            -O z \
            -o "$out_vcf" \
            "$anno_vcf_gz"

        bcftools index -f --tbi "$out_vcf"
        final=$(bcftools view -H "$out_vcf" | wc -l)
    fi

    echo "  Final variants: $final"

    # -----------------------------------------------------------------------
    # Collect per-step counts for the summary CSV
    # Parsing from awk stderr output saved to file
    # -----------------------------------------------------------------------

    n_sex_removed=$(grep "Removed sex chr" "${CELLLRK_DIR}/${sample}.filter_stats.txt" | awk '{print $NF}' || echo 0)
    n_cv_removed=$(grep  "Removed ClinVar" "${CELLLRK_DIR}/${sample}.filter_stats.txt" | awk '{print $NF}' || echo 0)
    n_af_removed=$(grep  "Removed pop AF"  "${CELLLRK_DIR}/${sample}.filter_stats.txt" | awk '{print $NF}' || echo 0)

    after_sex=$((total - n_sex_removed))
    after_clinvar=$((after_sex - n_cv_removed))
    after_popAF=$((after_clinvar - n_af_removed))

    echo "$sample,$total,$after_sex,$after_clinvar,$after_popAF,$final" \
        >> "$SUMMARY_CSV"

    echo "[OK] Done: $sample"

done

# ---------------------------------------------------------------------------
# DONE
# ---------------------------------------------------------------------------

echo ""
echo "=== cellLRK-filter COMPLETE: $(date) ==="
echo ""
echo "Filtered VCFs : $CELLLRK_DIR"
echo "Summary CSV   : $SUMMARY_CSV"
echo ""
cat "$SUMMARY_CSV"