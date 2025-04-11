#!/bin/bash

# Check missing files and filter samplesheet
input="samplesheet.csv"
output="filtered_samplesheet.csv"
missing_log="missing_files.txt"

# Initialize counts and missing files list
missing_fastq1=0
missing_fastq2=0
> "$missing_log"  # Clear previous log

# Process CSV while preserving header
{
    read header  # Save header
    echo "$header" > "$output"
    
    while IFS=, read -r patient status sample lane fastq1 fastq2; do
        missing=0
        
        # Check fastq_1
        if [ ! -f "$fastq1" ]; then
            echo "$fastq1" >> "$missing_log"
            ((missing_fastq1++))
            missing=1
        fi
        
        # Check fastq_2
        if [ ! -f "$fastq2" ]; then
            echo "$fastq2" >> "$missing_log"
            ((missing_fastq2++))
            missing=1
        fi
        
        # Keep row only if both files exist
        [ "$missing" -eq 0 ] && echo "$patient,$status,$sample,$lane,$fastq1,$fastq2" >> "$output"
        
    done
} < "$input"

# Print summary
echo "Missing fastq_1 files: $missing_fastq1"
echo "Missing fastq_2 files: $missing_fastq2"
echo "Full list of missing files saved to: $missing_log"
