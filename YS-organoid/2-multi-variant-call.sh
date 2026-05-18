#!/bin/bash

#PBS -l select=2:ncpus=64:mem=128gb
#PBS -l walltime=01:00:00
#PBS -P 11003581
#PBS -N merge-multi-caller
#PBS -j oe

set -euo pipefail
shopt -s nullglob

cd $PBS_O_WORKDIR

echo "=== JOB STARTED: $(date) ==="

module load bcftools/1.15.1

############################################################
# REFERENCE
############################################################

REF="/home/project/11003581/Ref/Homo_sapiens/GATK/hg38/Homo_sapiens_assembly38.fasta"

############################################################
# DIRECTORIES
############################################################

base_dir="/home/users/nus/ash.ps/scratch/YS-organoid-downstream/analysis2"

SNV_DIR="$base_dir/raw-vcfs/strelka/snv"
INDEL_DIR="$base_dir/raw-vcfs/strelka/indel"
STRELKA_MERGED_DIR="$base_dir/raw-vcfs/strelka/merged"

MUTECT2_DIR="$base_dir/raw-vcfs/mutect2"
INTERSECTION_DIR="$base_dir/common-snvs"
FINAL_DIR="$base_dir/final-common"

TEMP_DIR="/tmp/vcf_processing_$$"

mkdir -p "$STRELKA_MERGED_DIR" "$INTERSECTION_DIR" "$FINAL_DIR" "$TEMP_DIR"

CSV_FILE="$FINAL_DIR/variant_counts.csv"
echo "samplename,mutect2-uniq,strelka-uniq,common" > "$CSV_FILE"

############################################################
# SANITY CHECK
############################################################

echo "=== SANITY CHECK ==="

ls "${SNV_DIR}"/*.vcf.gz | head
ls "${INDEL_DIR}"/*.vcf.gz | head
ls "${MUTECT2_DIR}"/*.vcf.gz | head

############################################################
# STEP 1: CONCAT STRELKA SNVs + INDELs
############################################################

echo "=== STEP 1: Concatenating Strelka SNVs + INDELs ==="

for mutect2_file in "${MUTECT2_DIR}"/*.mutect2.filtered.vcf.gz; do
    
    sample_name=$(basename "$mutect2_file" .mutect2.filtered.vcf.gz)

    snv_file="${SNV_DIR}/${sample_name}.strelka.somatic_snvs.vcf.gz"
    indel_file="${INDEL_DIR}/${sample_name}.strelka.somatic_indels.vcf.gz"
    merged_file="${STRELKA_MERGED_DIR}/${sample_name}.merged.vcf.gz"

    if [[ ! -f "$snv_file" ]]; then
        echo "[ERROR] Missing SNV for $sample_name"
        exit 1
    fi

    if [[ ! -f "$indel_file" ]]; then
        echo "[ERROR] Missing INDEL for $sample_name"
        exit 1
    fi

    if [[ ! -f "$merged_file" ]]; then
        
        echo "[INFO] Concatenating: $sample_name"
        
        bcftools concat -a -O z -o "$merged_file" \
            "$snv_file" "$indel_file"
        
        # Sort after concat
        bcftools sort "$merged_file" -O z -o "${merged_file%.vcf.gz}.sorted.vcf.gz"
        mv "${merged_file%.vcf.gz}.sorted.vcf.gz" "$merged_file"
        
        bcftools index "$merged_file"
    fi

    echo "[OK] Ready: $sample_name"
done

############################################################
# STEP 2: FILTER + NORMALIZE + INTERSECT
############################################################

echo "=== STEP 2: Normalize + Intersection ==="

for mutect2_file in "${MUTECT2_DIR}"/*.mutect2.filtered.vcf.gz; do
    
    sample_name=$(basename "$mutect2_file" .mutect2.filtered.vcf.gz)
    strelka_file="${STRELKA_MERGED_DIR}/${sample_name}.merged.vcf.gz"
    out_dir="${INTERSECTION_DIR}/${sample_name}_intersection"
    
    if [[ ! -f "$strelka_file" ]]; then
        echo "[ERROR] Missing merged Strelka for $sample_name"
        exit 1
    fi
    
    echo "Processing: $sample_name"
    
    ########################################################
    # CLEAN OLD OUTPUT (IMPORTANT)
    ########################################################
    
    rm -rf "$out_dir"
    
    ########################################################
    # PASS FILTER
    ########################################################
    
    bcftools view -f PASS "$mutect2_file" -Oz -o "$TEMP_DIR/${sample_name}.mutect2.pass.vcf.gz"
    bcftools view -f PASS "$strelka_file" -Oz -o "$TEMP_DIR/${sample_name}.strelka.pass.vcf.gz"
    
    bcftools index "$TEMP_DIR/${sample_name}.mutect2.pass.vcf.gz"
    bcftools index "$TEMP_DIR/${sample_name}.strelka.pass.vcf.gz"
    
    ########################################################
    # NORMALIZATION
    ########################################################
    
    bcftools norm -m -both -f "$REF" \
        "$TEMP_DIR/${sample_name}.mutect2.pass.vcf.gz" \
        -Oz -o "$TEMP_DIR/${sample_name}.mutect2.norm.vcf.gz"
    
    bcftools norm -m -both -f "$REF" \
        "$TEMP_DIR/${sample_name}.strelka.pass.vcf.gz" \
        -Oz -o "$TEMP_DIR/${sample_name}.strelka.norm.vcf.gz"
    
    bcftools index "$TEMP_DIR/${sample_name}.mutect2.norm.vcf.gz"
    bcftools index "$TEMP_DIR/${sample_name}.strelka.norm.vcf.gz"
    
    ########################################################
    # INTERSECTION (FIXED: COMPRESSED OUTPUT)
    ########################################################
    
    bcftools isec -p "$out_dir" -Oz \
        "$TEMP_DIR/${sample_name}.mutect2.norm.vcf.gz" \
        "$TEMP_DIR/${sample_name}.strelka.norm.vcf.gz"
    
    ########################################################
    # INDEX ISEC OUTPUTS (IMPORTANT)
    ########################################################
    
    for f in "$out_dir"/*.vcf.gz; do
        bcftools index "$f"
    done
    
    ########################################################
    # VALIDATION (FAIL FAST IF BROKEN)
    ########################################################
    
    if [[ ! -f "$out_dir/0002.vcf.gz" ]]; then
        echo "[WARN] No intersection variants for $sample_name"
    fi
    
    ########################################################
    # CLEAN TEMP
    ########################################################
    
    rm -f "$TEMP_DIR/${sample_name}."*
    
    echo "[OK] Done: $sample_name"
done

############################################################
# STEP 3: EXTRACT COMMON + COUNTS
############################################################

echo "=== STEP 3: Extracting common variants + counts ==="

for sample_dir in "${INTERSECTION_DIR}"/*_intersection; do
    
    sample_name=$(basename "$sample_dir" _intersection)
    
    echo "Counting: $sample_name"
    
    cp "$sample_dir/0002.vcf.gz" "$FINAL_DIR/${sample_name}_common.vcf.gz"
    
    mutect2_uniq=$(zgrep -v '^#' "$sample_dir/0000.vcf.gz" | wc -l)
    strelka_uniq=$(zgrep -v '^#' "$sample_dir/0001.vcf.gz" | wc -l)
    common=$(zgrep -v '^#' "$sample_dir/0002.vcf.gz" | wc -l)
    
    echo "$sample_name,$mutect2_uniq,$strelka_uniq,$common" >> "$CSV_FILE"
    
    echo "[OK] Counted: $sample_name"
done

############################################################
# CLEANUP
############################################################

rm -rf "$TEMP_DIR"

echo "=== JOB COMPLETED: $(date) ==="
echo "Results:"
echo " - Common VCFs: $FINAL_DIR"
echo " - Counts CSV: $CSV_FILE"