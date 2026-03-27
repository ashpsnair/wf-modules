# By removing reads mapped to non-canonical or extraneous contigs and retaining only standard human chromosomes, 
#!/bin/bash
#PBS -N clean_RNA_BAMs
#PBS -l select=1:ncpus=128:mem=256gb
#PBS -l walltime=04:00:00
#PBS -j oe

set -euo pipefail

module load samtools/1.15.1

INPUT_DIR="/home/users/nus/ash.ps/scratch/data-transfer/TCGA_RNAseq_with_WGS/data_for_Ash_test/RNAseq"
METADATA="/home/users/nus/ash.ps/scratch/RNA-mut/RNA/metadata_rnaseq.txt"
OUTPUT_DIR="/home/users/nus/ash.ps/scratch/RNA-mut/RNA/cleaned_bams"

LOGFILE="${OUTPUT_DIR}/cleaning_$(date +%F_%H-%M-%S).log"

mkdir -p "$OUTPUT_DIR"
cd "$INPUT_DIR"

exec > >(tee -a "$LOGFILE") 2>&1

echo "=== JOB STARTED $(date) ==="

REGIONS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 \
chr20 chr21 chr22 chrX chrY"

########################################
# SAFE RESOURCE SETTINGS
########################################
THREADS_PER_JOB=8
MAX_JOBS=4

########################################
process_bam() {

    BAM_PATH="$1"
    BAM_NAME=$(basename "$BAM_PATH")

    LINE=$(grep -F "$BAM_NAME" "$METADATA" || true)
    if [ -z "$LINE" ]; then
        echo "WARNING: No metadata for $BAM_NAME"
        return
    fi

    CASE_ID=$(echo "$LINE" | cut -f1)
    TYPE=$(echo "$LINE" | cut -f3)

    SAMPLE="${CASE_ID}_${TYPE}"
    FINAL_BAM="${OUTPUT_DIR}/${SAMPLE}.bam"
    FINAL_INDEX="${FINAL_BAM}.csi"

    ########################################
    # Skip if already valid
    ########################################
    if [ -f "$FINAL_BAM" ] && [ -f "$FINAL_INDEX" ]; then
        if samtools quickcheck "$FINAL_BAM" && samtools idxstats "$FINAL_BAM" > /dev/null 2>&1; then
            echo ">>> Skipping $SAMPLE (already valid)"
            return
        else
            echo ">>> Existing $SAMPLE is corrupted → rebuilding"
            rm -f "$FINAL_BAM" "$FINAL_INDEX"
        fi
    fi

    ########################################
    # Unique temp dir
    ########################################
    TMPDIR=$(mktemp -d "${OUTPUT_DIR}/${SAMPLE}_XXXX")
    TMP_BAM="${TMPDIR}/tmp.bam"
    SORTED_BAM="${TMPDIR}/sorted.bam"
    HEADER="${TMPDIR}/header.sam"
    TEST_BAM="${TMPDIR}/test.bam"

    echo ">>> Processing $SAMPLE"

    ########################################
    # Step 1: Extract
    ########################################
    if ! samtools view -@ "$THREADS_PER_JOB" -b "$BAM_PATH" $REGIONS > "$TMP_BAM"; then
        echo "ERROR: extraction failed for $SAMPLE"
        rm -rf "$TMPDIR"
        return
    fi

    ########################################
    # Step 2: Sort
    ########################################
    if ! samtools sort -@ "$THREADS_PER_JOB" -m 2G -o "$SORTED_BAM" "$TMP_BAM"; then
        echo "ERROR: sort failed for $SAMPLE"
        rm -rf "$TMPDIR"
        return
    fi

    ########################################
    # Step 3: Validate BEFORE proceeding
    ########################################
    if ! samtools view "$SORTED_BAM" > /dev/null 2>&1; then
        echo "ERROR: BAM corrupted after sort → $SAMPLE"
        rm -rf "$TMPDIR"
        return
    fi

    ########################################
    # Step 4: Clean header
    ########################################
    samtools view -H "$SORTED_BAM" | \
    grep -E "^@HD|^@PG|^@RG|^@CO|^@SQ\s+SN:(chr([1-9]|1[0-9]|2[0-2]|X|Y))\b" \
    > "$HEADER"

    samtools reheader "$HEADER" "$SORTED_BAM" > "$TEST_BAM"

    ########################################
    # Step 5: Validate AGAIN
    ########################################
    if ! samtools view "$TEST_BAM" > /dev/null 2>&1; then
        echo "ERROR: BAM corrupted after reheader → $SAMPLE"
        rm -rf "$TMPDIR"
        return
    fi

    ########################################
    # Step 6: Move to final
    ########################################
    mv "$TEST_BAM" "$FINAL_BAM"

    ########################################
    # Step 7: CSI index (correct)
    ########################################
    if ! samtools index -c "$FINAL_BAM"; then
        echo "ERROR: indexing failed for $SAMPLE"
        rm -f "$FINAL_BAM"
        rm -rf "$TMPDIR"
        return
    fi

    ########################################
    # Final validation
    ########################################
    samtools quickcheck "$FINAL_BAM"
    samtools idxstats "$FINAL_BAM" > /dev/null

    ########################################
    # Cleanup
    ########################################
    rm -rf "$TMPDIR"

    echo ">>> SUCCESS $SAMPLE"
}

export -f process_bam
export METADATA OUTPUT_DIR THREADS_PER_JOB REGIONS

########################################
# RUN
########################################
find . -name "*.bam" | \
xargs -n 1 -P "$MAX_JOBS" -I {} bash -c 'process_bam "$@"' _ {}

echo "=== ALL DONE $(date) ==="