#!/bin/bash

#PBS -l select=3:ncpus=64:mem=256g
#PBS -l walltime=12:00:00
#PBS -P 11003581
#PBS -N concat-b1
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

# Set the directory where you want to save the concatenated files
OUTPUT_DIR="/home/users/nus/ash.ps/scratch/NCCS-MASH/analysis/tumor/new-wf/concat-fastqs/"

# Create the output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Declare associative arrays to store fastq file paths for each sample
declare -A READ1_FILES
declare -A READ2_FILES

# Read the sample sheet and store the fastq file paths
while IFS=',' read -r patient status sample lane fastq_1 fastq_2; do
    # Skip the header line
    if [[ "$sample" == "sample" ]]; then
        continue
    fi

    # Store the fastq file paths in the associative arrays
    READ1_FILES[$sample]+="$fastq_1 "
    READ2_FILES[$sample]+="$fastq_2 "
done < samplesheet.csv

# Loop through the unique sample IDs
for sample in "${!READ1_FILES[@]}"; do
    # Create variables for the read 1 and read 2 file names
    SAMPLE_ID=${sample}
    READ1_OUTFILE="${OUTPUT_DIR}/${SAMPLE_ID}_R1.fastq.gz"
    READ2_OUTFILE="${OUTPUT_DIR}/${SAMPLE_ID}_R2.fastq.gz"

    # Check if the output files already exist, if so, skip to the next sample
    if [ -f "$READ1_OUTFILE" ] && [ -f "$READ2_OUTFILE" ]; then
        echo "Output files for sample $SAMPLE_ID already exist. Skipping..."
        continue
    fi

    # Get the number of files being concatenated
    READ1_COUNT=$(echo ${READ1_FILES[$sample]} | wc -w)
    READ2_COUNT=$(echo ${READ2_FILES[$sample]} | wc -w)

    # Calculate the size of input files before concatenation for read 1
    READ1_SIZE_BEFORE=0
    for file in ${READ1_FILES[$sample]}; do
        READ1_SIZE_BEFORE=$((READ1_SIZE_BEFORE + $(stat -c "%s" "$file")))
    done

    # Calculate the size of input files before concatenation for read 2
    READ2_SIZE_BEFORE=0
    for file in ${READ2_FILES[$sample]}; do
        READ2_SIZE_BEFORE=$((READ2_SIZE_BEFORE + $(stat -c "%s" "$file")))
    done

    # Concatenate the fastq files for read 1 and read 2 using cat
    cat ${READ1_FILES[$sample]} > "$READ1_OUTFILE"
    cat ${READ2_FILES[$sample]} > "$READ2_OUTFILE"

    # Get the size of the final concatenated file for read 1
    READ1_SIZE_AFTER=$(stat -c "%s" "$READ1_OUTFILE")

    # Get the size of the final concatenated file for read 2
    READ2_SIZE_AFTER=$(stat -c "%s" "$READ2_OUTFILE")

    echo "Sample: $SAMPLE_ID"
    echo "  Read 1: Concatenated $READ1_COUNT files, Size before: $READ1_SIZE_BEFORE bytes, Size after: $READ1_SIZE_AFTER bytes"
    echo "  Read 2: Concatenated $READ2_COUNT files, Size before: $READ2_SIZE_BEFORE bytes, Size after: $READ2_SIZE_AFTER bytes"
    echo "Concatenated fastq files for sample $SAMPLE_ID"
done

echo "Finished concatenating fastq files for all samples."
