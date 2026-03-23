# Cleaned the RNA-seq BAM files because they were aligned to a reference that included additional non-human contigs (e.g., **CMV**), 
# whereas our downstream pipeline (nf-core/sarek) expects strict compatibility with the standard human reference genome (hg38). 
# This mismatch causes errors during analysis since the reference FASTA does not contain these extra contigs present in the BAM headers. 
# By removing reads mapped to non-canonical or extraneous contigs and retaining only standard human chromosomes, 
# we ensure consistency between the BAM files and the reference genome, allowing the pipeline to run correctly while also reducing potential noise from irrelevant or contaminant sequences.




#!/bin/bash
#PBS -N clean_RNA_BAMs
#PBS -l select=1:ncpus=128:mem=256gb
#PBS -l walltime=08:00:00
#PBS -j oe

########################################
# Load modules
########################################
module load samtools/1.15.1

########################################
# Input / Output directories
########################################
INPUT_DIR="/home/users/nus/ash.ps/scratch/RNA-mut/RNA-bams"
OUTPUT_DIR="${INPUT_DIR}/cleaned_bams"

########################################
# Create output directory
########################################
mkdir -p $OUTPUT_DIR

########################################
# Move to input directory
########################################
cd $INPUT_DIR

########################################
# Loop through BAMs
########################################
for BAM in *.bam; do
    
    SAMPLE=$(basename $BAM .bam)
    OUT_BAM="${OUTPUT_DIR}/${SAMPLE}_cleaned.bam"

    echo "Processing $SAMPLE..."

    ########################################
    # Extract contigs excluding CMV
    ########################################
    CONTIGS=$(samtools idxstats $BAM | cut -f1 | grep -v "CMV")

    ########################################
    # Subset BAM
    ########################################
    samtools view -@ 8 -b $BAM $CONTIGS > $OUT_BAM

    ########################################
    # Index
    ########################################
    samtools index $OUT_BAM

    echo "Done $SAMPLE"

done

echo "All BAMs cleaned and ready for Sarek"