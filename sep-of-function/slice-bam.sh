#!/bin/bash
#PBS -l select=1:ncpus=128:mem=256g
#PBS -l walltime=04:00:00
#PBS -P 11003581
#PBS -N slice-bam-rad51
#PBS -j oe

module load samtools

# Reference fasta path (adjust as needed)
REFERENCE="/home/project/11003581/Ref/sarek-refs/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"

# Input base directory containing folders with cram files
BASE_INPUT_DIR="/home/users/nus/ash.ps/scratch/sep-fun-analysis/RAD51_S181P/preprocessing/recalibrated"

# Output directory for chr15 BAMs
OUTPUT_DIR="/home/users/nus/ash.ps/scratch/sep-fun-analysis/RAD51_S181P/preprocessing/chr15-bams"

mkdir -p "$OUTPUT_DIR"

# Find all cram files recursively under BASE_INPUT_DIR
find "$BASE_INPUT_DIR" -type f -name "*.cram" | while read -r cramfile; do
    # Get just filename without path and extension
    filename=$(basename "$cramfile" .cram)
    
    # Output BAM path
    outbam="$OUTPUT_DIR/${filename}_chr15.bam"
    
    echo "Processing $cramfile -> $outbam"
    
    # Extract chr15 reads to BAM
    samtools view -T "$REFERENCE" -r 15 -b -o "$outbam" "$cramfile"
    
    # Optional: Index BAM
    samtools index "$outbam"
done
