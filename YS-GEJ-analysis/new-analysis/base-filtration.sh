#!/bin/bash

#PBS -l select=1:ncpus=64
#PBS -l walltime=2:00:00
#PBS -P 11003581
#PBS -N combined-filter-vcfs
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load bcftools/1.15.1

# Set input and output directories
input_dir="/home/users/nus/ash.ps/scratch/YS-analysis/VCFs/"
output_dir="/home/users/nus/ash.ps/scratch/YS-analysis/base-filtered-vcfs"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Part 1: Filter VCFs based on PASS filter
echo "Starting Part 1: Filtering VCFs based on PASS filter"
for file in "$input_dir"/*.mutect2.filtered.vcf.gz; do
    base_name=$(basename "$file" .mutect2.filtered.vcf.gz)
    output_file="$output_dir/${base_name}_pass_filtered.vcf"
    bcftools view -i 'FILTER="PASS"' "$file" -o "$output_file"
    echo "Processed: $file -> $output_file"
done

