#!/bin/bash

#PBS -l select=1:ncpus=16:mem=64g
#PBS -l walltime=04:00:00
#PBS -P 11003581
#PBS -N filter2
#PBS -j oe
#PBS -m ae

set -euo pipefail

cd $PBS_O_WORKDIR

module load bcftools/1.15.1

echo "=== FILTER2 STARTED: $(date) ==="

############################################################
# PATHS
############################################################

BASE_DIR="/home/users/nus/ash.ps/scratch/YS-organoid-downstream/analysis1.1"

INPUT_DIR="${BASE_DIR}/03_filtered"

BINOMIAL_DIR="${BASE_DIR}/04_binomial_filtered"

REMOVE_ALL_COMMON_DIR="${BASE_DIR}/05_remove_common_all"

FINAL_DIR="${BASE_DIR}/06_final_hardfiltered"

TEMP_DIR="/tmp/filter2_$$"

mkdir -p \
    "$BINOMIAL_DIR" \
    "$REMOVE_ALL_COMMON_DIR" \
    "$FINAL_DIR" \
    "$TEMP_DIR"

############################################################
# SUMMARY CSV
############################################################

SUMMARY="${BASE_DIR}/variants_filter2_summary.csv"

echo "sample,after_03_filtered,after_binomial_filter,after_remove_all_common,final_after_individual_filter" \
> "$SUMMARY"

############################################################
# STEP 1
# BINOMIAL / HIGH-VAF FILTER
############################################################

echo ""
echo "=== STEP 1: Binomial-like VAF filtering ==="

MIN_VAF=0.20

for vcf in "${INPUT_DIR}"/*.vcf.gz; do

    sample=$(basename "$vcf" .common.filtered.vcf.gz)

    echo ""
    echo "--- $sample ---"

    initial=$(bcftools view -H "$vcf" | wc -l)

    out="${BINOMIAL_DIR}/${sample}.vaf20.vcf.gz"

    bcftools filter \
        --include '(FORMAT/AD[1:1] / FORMAT/DP[1]) >= '"$MIN_VAF" \
        -O z \
        -o "$out" \
        "$vcf"

    bcftools index -f --tbi "$out"

    after_binomial=$(bcftools view -H "$out" | wc -l)

    echo "03_filtered: $initial"
    echo "After VAF20: $after_binomial"

    echo "$sample,$initial,$after_binomial,PENDING,PENDING" \
    >> "$SUMMARY"

done

############################################################
# STEP 2
# FIND COMMON ACROSS ALL SAMPLES
############################################################

echo ""
echo "=== STEP 2: Finding variants common across ALL samples ==="

ALL_VCFS=(${BINOMIAL_DIR}/*.vcf.gz)

bcftools isec \
    -n=4 \
    -c all \
    -w1 \
    -O z \
    -o "${TEMP_DIR}/common_all.vcf.gz" \
    "${ALL_VCFS[@]}"

bcftools index -f --tbi "${TEMP_DIR}/common_all.vcf.gz"

common_all_count=$(bcftools view -H "${TEMP_DIR}/common_all.vcf.gz" | wc -l)

echo "Variants common across ALL samples: $common_all_count"

############################################################
# STEP 3
# REMOVE COMMON ACROSS ALL
############################################################

echo ""
echo "=== STEP 3: Removing variants common across ALL samples ==="

for vcf in "${BINOMIAL_DIR}"/*.vcf.gz; do

    sample=$(basename "$vcf" .vaf20.vcf.gz)

    echo ""
    echo "--- $sample ---"

    out="${REMOVE_ALL_COMMON_DIR}/${sample}.unique.vcf.gz"

    bcftools isec \
        -C \
        -w1 \
        -O z \
        -o "$out" \
        "$vcf" \
        "${TEMP_DIR}/common_all.vcf.gz"

    bcftools index -f --tbi "$out"

    after_all_common=$(bcftools view -H "$out" | wc -l)

    echo "After removing all-common: $after_all_common"

    sed -i "s/^${sample},\([^,]*,[^,]*\),PENDING,PENDING$/${sample},\1,${after_all_common},PENDING/" \
    "$SUMMARY"

done

############################################################
# STEP 4
# FIND COMMON WITHIN INDIVIDUALS
############################################################

echo ""
echo "=== STEP 4: Finding variants common within individuals ==="

###########################################################
# GEJ07
###########################################################

GEJ07_FILES=(
"${REMOVE_ALL_COMMON_DIR}/GEJ07_Wk4_2_vs_GEJ07_Wk4_0.unique.vcf.gz"
"${REMOVE_ALL_COMMON_DIR}/GEJ07_Wk4_5_vs_GEJ07_Wk4_0.unique.vcf.gz"
)

bcftools isec \
    -n=2 \
    -c all \
    -w1 \
    -O z \
    -o "${TEMP_DIR}/GEJ07_common.vcf.gz" \
    "${GEJ07_FILES[@]}"

bcftools index -f --tbi "${TEMP_DIR}/GEJ07_common.vcf.gz"

GEJ07_common_count=$(bcftools view -H "${TEMP_DIR}/GEJ07_common.vcf.gz" | wc -l)

echo "GEJ07 common variants: $GEJ07_common_count"

###########################################################
# GEJ09
###########################################################

GEJ09_FILES=(
"${REMOVE_ALL_COMMON_DIR}/GEJ09_Wk4_2_vs_GEJ09_Wk4_0.unique.vcf.gz"
"${REMOVE_ALL_COMMON_DIR}/GEJ09_Wk4_5_vs_GEJ09_Wk4_0.unique.vcf.gz"
)

bcftools isec \
    -n=2 \
    -c all \
    -w1 \
    -O z \
    -o "${TEMP_DIR}/GEJ09_common.vcf.gz" \
    "${GEJ09_FILES[@]}"

bcftools index -f --tbi "${TEMP_DIR}/GEJ09_common.vcf.gz"

GEJ09_common_count=$(bcftools view -H "${TEMP_DIR}/GEJ09_common.vcf.gz" | wc -l)

echo "GEJ09 common variants: $GEJ09_common_count"

############################################################
# STEP 5
# REMOVE INDIVIDUAL-COMMON
############################################################

echo ""
echo "=== STEP 5: Removing variants common within individuals ==="

for vcf in "${REMOVE_ALL_COMMON_DIR}"/*.vcf.gz; do

    sample=$(basename "$vcf" .unique.vcf.gz)

    echo ""
    echo "--- $sample ---"

    if [[ "$sample" == GEJ07* ]]; then
        common_vcf="${TEMP_DIR}/GEJ07_common.vcf.gz"
    else
        common_vcf="${TEMP_DIR}/GEJ09_common.vcf.gz"
    fi

    out="${FINAL_DIR}/${sample}.hardfiltered.vcf.gz"

    bcftools isec \
        -C \
        -w1 \
        -O z \
        -o "$out" \
        "$vcf" \
        "$common_vcf"

    bcftools index -f --tbi "$out"

    final=$(bcftools view -H "$out" | wc -l)

    echo "Final variants: $final"

    sed -i "s/^${sample},\(.*\),PENDING$/${sample},\1,${final}/" \
    "$SUMMARY"

done

############################################################
# FINAL
############################################################

echo ""
echo "=== FILTER2 COMPLETE: $(date) ==="

echo ""
echo "Directories:"
echo "  04_binomial_filtered  : $BINOMIAL_DIR"
echo "  05_remove_common_all  : $REMOVE_ALL_COMMON_DIR"
echo "  06_final_hardfiltered : $FINAL_DIR"

echo ""
echo "Summary:"
cat "$SUMMARY"

############################################################
# CLEANUP
############################################################

rm -rf "$TEMP_DIR"