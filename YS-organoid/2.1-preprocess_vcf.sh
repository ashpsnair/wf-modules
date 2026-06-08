#!/bin/bash

#PBS -l select=2:ncpus=32:mem=64g
#PBS -l walltime=04:00:00
#PBS -P 11003581
#PBS -N YS-vcf-preprocessing
#PBS -j oe


set -euo pipefail
shopt -s nullglob

cd $PBS_O_WORKDIR

echo "=== VCF PRE-PROCESSING STARTED: $(date) ==="

module load bcftools/1.15.1

############################################################
# PATHS
############################################################

REF="/home/project/11003581/Ref/sarek-refs/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"

BASE_DIR="/home/users/nus/ash.ps/scratch/YS-organoid-downstream/analysis1.1"

RAW_VCF_DIR="${BASE_DIR}/raw-vcfs"

MUTECT2_DIR="${RAW_VCF_DIR}/mutect2"

STRELKA_SNV_DIR="${RAW_VCF_DIR}/strelka/snv"
STRELKA_INDEL_DIR="${RAW_VCF_DIR}/strelka/indel"

############################################################
# OUTPUT DIRECTORIES
############################################################

STRELKA_MERGED_DIR="${BASE_DIR}/01_strelka_merged"
INTERSECTION_DIR="${BASE_DIR}/02_intersections"
FILTERED_DIR="${BASE_DIR}/03_filtered"
LOG_DIR="${BASE_DIR}/logs"

TEMP_DIR="/tmp/vcf_preproc_$$"

############################################################
# CLEAN OLD OUTPUTS
############################################################

rm -rf \
    "$STRELKA_MERGED_DIR" \
    "$INTERSECTION_DIR" \
    "$FILTERED_DIR"

mkdir -p \
    "$STRELKA_MERGED_DIR" \
    "$INTERSECTION_DIR" \
    "$FILTERED_DIR" \
    "$LOG_DIR" \
    "$TEMP_DIR"

############################################################
# SUMMARY CSV
############################################################

SUMMARY_CSV="${BASE_DIR}/preprocessing_summary.csv"

echo "sample,mutect2_pass,strelka_pass,mutect2_only,strelka_only,common_raw,common_filtered" \
> "$SUMMARY_CSV"

############################################################
# SANITY CHECK
############################################################

echo ""
echo "=== SANITY CHECK ==="

echo "--- Mutect2 VCFs ---"
find "$MUTECT2_DIR" -name "*.mutect2.filtered.vcf.gz" | sort

echo "--- Strelka SNV VCFs ---"
find "$STRELKA_SNV_DIR" -name "*.vcf.gz" | sort

echo "--- Strelka INDEL VCFs ---"
find "$STRELKA_INDEL_DIR" -name "*.vcf.gz" | sort

############################################################
# STEP 1
############################################################

echo ""
echo "=== STEP 1: Concatenating Strelka SNVs + INDELs ==="

n_processed=0

while read mutect2_file; do

    n_processed=$((n_processed + 1))

    sample_name=$(basename "$mutect2_file" .mutect2.filtered.vcf.gz)

    snv_file="${STRELKA_SNV_DIR}/${sample_name}.strelka.somatic_snvs.vcf.gz"

    indel_file="${STRELKA_INDEL_DIR}/${sample_name}.strelka.somatic_indels.vcf.gz"

    merged_file="${STRELKA_MERGED_DIR}/${sample_name}.strelka.merged.vcf.gz"

    echo ""
    echo "--- $sample_name ---"

    if [[ ! -f "$snv_file" ]]; then
        echo "[ERROR] Missing SNV file:"
        echo "$snv_file"
        continue
    fi

    if [[ ! -f "$indel_file" ]]; then
        echo "[ERROR] Missing INDEL file:"
        echo "$indel_file"
        continue
    fi

    ########################################################
    # CONCAT
    ########################################################

    bcftools concat \
        -a \
        -O z \
        -o "$TEMP_DIR/${sample_name}.concat.vcf.gz" \
        "$snv_file" \
        "$indel_file"

    ########################################################
    # SORT
    ########################################################

    bcftools sort \
        -O z \
        -o "$merged_file" \
        "$TEMP_DIR/${sample_name}.concat.vcf.gz"

    ########################################################
    # INDEX
    ########################################################

    bcftools index -f --tbi "$merged_file"

    rm -f "$TEMP_DIR/${sample_name}.concat.vcf.gz"

    n=$(bcftools view -H "$merged_file" | wc -l)

    echo "[OK] Strelka merged: $n variants"

done < <(find "$MUTECT2_DIR" -name "*.mutect2.filtered.vcf.gz")

echo ""
echo "STEP 1 processed samples: $n_processed"

############################################################
# STEP 2
############################################################

echo ""
echo "=== STEP 2: PASS Filter + Intersect ==="

while read mutect2_file; do

    sample_name=$(basename "$mutect2_file" .mutect2.filtered.vcf.gz)

    strelka_file="${STRELKA_MERGED_DIR}/${sample_name}.strelka.merged.vcf.gz"

    isec_dir="${INTERSECTION_DIR}/${sample_name}"

    mkdir -p "$isec_dir"

    echo ""
    echo "--- $sample_name ---"

    if [[ ! -f "$strelka_file" ]]; then
        echo "[ERROR] Missing merged Strelka:"
        echo "$strelka_file"
        continue
    fi

    ########################################################
    # PASS FILTER
    ########################################################

    bcftools view \
        -f PASS \
        -O z \
        -o "$TEMP_DIR/${sample_name}.mutect2.pass.vcf.gz" \
        "$mutect2_file"

    bcftools index -f --tbi \
        "$TEMP_DIR/${sample_name}.mutect2.pass.vcf.gz"

    bcftools view \
        -f PASS \
        -O z \
        -o "$TEMP_DIR/${sample_name}.strelka.pass.vcf.gz" \
        "$strelka_file"

    bcftools index -f --tbi \
        "$TEMP_DIR/${sample_name}.strelka.pass.vcf.gz"

    ########################################################
    # NORMALISE
    ########################################################

    bcftools norm \
        -m -any \
        -f "$REF" \
        -O z \
        -o "$TEMP_DIR/${sample_name}.mutect2.norm.vcf.gz" \
        "$TEMP_DIR/${sample_name}.mutect2.pass.vcf.gz"

    bcftools index -f --tbi \
        "$TEMP_DIR/${sample_name}.mutect2.norm.vcf.gz"

    bcftools norm \
        -m -any \
        -f "$REF" \
        -O z \
        -o "$TEMP_DIR/${sample_name}.strelka.norm.vcf.gz" \
        "$TEMP_DIR/${sample_name}.strelka.pass.vcf.gz"

    bcftools index -f --tbi \
        "$TEMP_DIR/${sample_name}.strelka.norm.vcf.gz"

    ########################################################
    # DEDUP
    ########################################################

    bcftools norm \
        -d exact \
        -O z \
        -o "$TEMP_DIR/${sample_name}.mutect2.dedup.vcf.gz" \
        "$TEMP_DIR/${sample_name}.mutect2.norm.vcf.gz"

    bcftools index -f --tbi \
        "$TEMP_DIR/${sample_name}.mutect2.dedup.vcf.gz"

    bcftools norm \
        -d exact \
        -O z \
        -o "$TEMP_DIR/${sample_name}.strelka.dedup.vcf.gz" \
        "$TEMP_DIR/${sample_name}.strelka.norm.vcf.gz"

    bcftools index -f --tbi \
        "$TEMP_DIR/${sample_name}.strelka.dedup.vcf.gz"

    ########################################################
    # INTERSECTION
    ########################################################

    bcftools isec \
        -p "$isec_dir" \
        -O z \
        "$TEMP_DIR/${sample_name}.mutect2.dedup.vcf.gz" \
        "$TEMP_DIR/${sample_name}.strelka.dedup.vcf.gz"

    for f in "$isec_dir"/*.vcf.gz; do
        bcftools index -f --tbi "$f"
    done

    ########################################################
    # COUNTS
    ########################################################

    m2_pass=$(bcftools view -H "$TEMP_DIR/${sample_name}.mutect2.pass.vcf.gz" | wc -l)

    sk_pass=$(bcftools view -H "$TEMP_DIR/${sample_name}.strelka.pass.vcf.gz" | wc -l)

    m2_only=$(bcftools view -H "$isec_dir/0000.vcf.gz" | wc -l)

    sk_only=$(bcftools view -H "$isec_dir/0001.vcf.gz" | wc -l)

    common=$(bcftools view -H "$isec_dir/0002.vcf.gz" | wc -l)

    echo "$sample_name,$m2_pass,$sk_pass,$m2_only,$sk_only,$common,PENDING" \
    >> "$SUMMARY_CSV"

    echo "[OK] Intersected"

done < <(find "$MUTECT2_DIR" -name "*.mutect2.filtered.vcf.gz")

############################################################
# STEP 3
############################################################

echo ""
echo "=== STEP 3: Quality Filtering ==="

MIN_TUMOR_DP=20
MIN_NORMAL_DP=10
MIN_ALT=5
MIN_AF=0.10

for isec_dir in "${INTERSECTION_DIR}"/*; do

    sample_name=$(basename "$isec_dir")

    common_vcf="${isec_dir}/0002.vcf.gz"

    out_vcf="${FILTERED_DIR}/${sample_name}.common.filtered.vcf.gz"

    echo ""
    echo "--- $sample_name ---"

    if [[ ! -f "$common_vcf" ]]; then
        echo "[WARN] Missing intersection VCF"
        continue
    fi

    ########################################################
    # FILTERING
    ########################################################

    bcftools filter \
        --include \
        'FORMAT/DP[1] >= '"$MIN_TUMOR_DP"' &&
         FORMAT/DP[0] >= '"$MIN_NORMAL_DP"' &&
         FORMAT/AD[1:1] >= '"$MIN_ALT"' &&
         (FORMAT/AD[1:1] / FORMAT/DP[1]) >= '"$MIN_AF" \
        -O z \
        -o "$out_vcf" \
        "$common_vcf"

    bcftools index -f --tbi "$out_vcf"

    ########################################################
    # COUNTS
    ########################################################

    before=$(bcftools view -H "$common_vcf" | wc -l)

    after=$(bcftools view -H "$out_vcf" | wc -l)

    sed -i "s/^${sample_name},\(.*\),PENDING$/${sample_name},\1,${after}/" \
    "$SUMMARY_CSV"

    echo "Before=$before After=$after"

done

############################################################
# FINAL
############################################################

echo ""
echo "=== PRE-PROCESSING COMPLETE: $(date) ==="

echo ""
echo "Merged:"
echo "$STRELKA_MERGED_DIR"

echo ""
echo "Intersections:"
echo "$INTERSECTION_DIR"

echo ""
echo "Filtered:"
echo "$FILTERED_DIR"

echo ""
echo "Summary:"
cat "$SUMMARY_CSV"

############################################################
# CLEANUP
############################################################

rm -rf "$TEMP_DIR"