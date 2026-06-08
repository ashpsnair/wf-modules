#!/bin/bash

#PBS -l select=1:ncpus=8:mem=32gb
#PBS -l walltime=02:00:00
#PBS -P 11003581
#PBS -N annovar_signature_filter
#PBS -j oe
#PBS -M your_email@domain.com
#PBS -m ae

set -eo pipefail

cd $PBS_O_WORKDIR

module load bcftools/1.15.1

echo "=== ANNOVAR SIGNATURE FILTER STARTED: $(date) ==="

############################################################
# PATHS
############################################################

BASE_DIR="/home/users/nus/ash.ps/scratch/YS-organoid-downstream/lrk-organoid"

ANNOTATED_DIR="/home/users/nus/ash.ps/scratch/YS-organoid-downstream/analysis1/annotated"

OUT_DIR="${BASE_DIR}/filtered-signature-vcfs"

LOG_DIR="${BASE_DIR}/logs"

BGZIP="/home/project/11003581/Tools/bin/bgzip"

mkdir -p "$OUT_DIR" "$LOG_DIR"

############################################################
# SUMMARY
############################################################

SUMMARY="${OUT_DIR}/signature_filter_summary.csv"

echo "sample,total_variants,final_variants" \
> "$SUMMARY"

############################################################
# PROCESS
############################################################

for sample_dir in ${ANNOTATED_DIR}/*; do

    sample=$(basename "$sample_dir")

    multianno="${sample_dir}/${sample}.hg38_multianno.txt"

    original_vcf="${sample_dir}/${sample}.hg38_multianno.vcf"

    compressed_vcf="${sample_dir}/${sample}.hg38_multianno.vcf.gz"

    echo ""
    echo "========================================="
    echo "Processing: ${sample}"
    echo "========================================="

    ########################################################
    # CHECK FILES
    ########################################################

    if [[ ! -f "$multianno" ]]; then
        echo "[WARN] Missing multianno file"
        continue
    fi

    if [[ ! -f "$original_vcf" ]]; then
        echo "[WARN] Missing VCF file"
        continue
    fi

    ########################################################
    # BGZIP + INDEX
    ########################################################

    if [[ ! -f "$compressed_vcf" ]]; then

        echo "[INFO] Compressing VCF"

        $BGZIP \
            -c "$original_vcf" \
            > "$compressed_vcf"

    fi

    if [[ ! -f "${compressed_vcf}.tbi" ]]; then

        echo "[INFO] Indexing VCF"

        bcftools index -f "$compressed_vcf"

    fi

    ########################################################
    # OUTPUT FILES
    ########################################################

    regions_bed="${OUT_DIR}/${sample}.regions.bed"

    filtered_vcf="${OUT_DIR}/${sample}.signature_filtered.vcf.gz"

    ########################################################
    # COUNTS
    ########################################################

    total=$(($(wc -l < "$multianno") - 1))

    echo "Total variants: $total"

    ########################################################
    # CREATE FILTERED BED REGIONS
    ########################################################

    awk -F '\t' '
    BEGIN{
        OFS="\t"
    }

    ########################################################
    # HEADER
    ########################################################

    NR==1{

        for(i=1;i<=NF;i++){

            if($i=="Chr")
                chr=i

            if($i=="Start")
                start=i

            if($i=="End")
                end=i

            if($i ~ /CLNSIG/)
                clinvar=i

            if($i ~ /ExAC_ALL/)
                exac=i

            if(tolower($i) ~ /gnomad/ && gnomad=="")
                gnomad=i
        }

        print "Detected columns:" > "/dev/stderr"
        print "Chr =", chr > "/dev/stderr"
        print "Start =", start > "/dev/stderr"
        print "End =", end > "/dev/stderr"
        print "CLNSIG =", clinvar > "/dev/stderr"
        print "ExAC =", exac > "/dev/stderr"
        print "gnomad =", gnomad > "/dev/stderr"

        next
    }

    ########################################################
    # FILTERING
    ########################################################

    {

        ####################################################
        # REMOVE chrX / chrY
        ####################################################

        if($chr=="chrX" || $chr=="chrY")
            next

        ####################################################
        # REMOVE BENIGN / LIKELY BENIGN
        ####################################################

        if($clinvar ~ /Benign/ || $clinvar ~ /Likely_benign/)
            next

        ####################################################
        # POPULATION AF FILTER
        ####################################################

        exac_af=0
        gnomad_af=0

        if($exac != "." && $exac != "")
            exac_af=$exac

        if($gnomad != "." && $gnomad != "")
            gnomad_af=$gnomad

        if(exac_af > 0.01 || gnomad_af > 0.01)
            next

        ####################################################
        # KEEP REGION
        ####################################################

        print $chr, $start-1, $end

    }

    ' "$multianno" \
    > "$regions_bed"

    ########################################################
    # FILTER ORIGINAL VCF
    ########################################################

    bcftools view \
        -R "$regions_bed" \
        -Oz \
        -o "$filtered_vcf" \
        "$compressed_vcf"

    ########################################################
    # INDEX FILTERED VCF
    ########################################################

    bcftools index -f "$filtered_vcf"

    ########################################################
    # FINAL COUNT
    ########################################################

    final=$(bcftools view -H "$filtered_vcf" | wc -l)

    echo "Final retained variants: $final"

    echo "$sample,$total,$final" \
    >> "$SUMMARY"

done

############################################################
# FINAL
############################################################

echo ""
echo "========================================="
echo "FILTERING COMPLETE"
echo "========================================="

echo ""
echo "Output directory:"
echo "$OUT_DIR"

echo ""
echo "Summary:"
cat "$SUMMARY"