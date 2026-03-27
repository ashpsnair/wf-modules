#!/bin/bash
#PBS -N BRCA2_slice
#PBS -l select=1:ncpus=64:mem=128gb
#PBS -l walltime=02:00:00
#PBS -j oe

set -euo pipefail

module load samtools/1.15.1

INPUT_BAM="/home/users/nus/ash.ps/scratch/LRK/SA/3031-II/3031-II.sorted.dedup.bam"
OUTPUT_DIR="/home/users/nus/ash.ps/scratch/LRK/bulk-analysis"

mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

SAMPLE=$(basename "$INPUT_BAM" .bam)

echo "Re-indexing BAM..."
samtools index -@ 64 "$INPUT_BAM"

echo "Extracting BRCA2..."
samtools view -@ 64 -b "$INPUT_BAM" \
13:32315474-32400266 \
> "${SAMPLE}_BRCA2.bam"

echo "Indexing output..."
samtools index -@ 64 "${SAMPLE}_BRCA2.bam"

echo "QC:"
samtools flagstat "${SAMPLE}_BRCA2.bam"