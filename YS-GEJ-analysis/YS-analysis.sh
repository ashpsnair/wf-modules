#!/bin/bash

#PBS -l select=1:ncpus=32
#PBS -l walltime=1:00:00
#PBS -P 11003581
#PBS -N filter-vcfs
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load bcftools/1.15.1

# Set input and output directories
input_dir="/home/users/nus/ash.ps/scratch/YS-analysis/VCFs/"
output_dir="/home/users/nus/ash.ps/scratch/YS-analysis/filter1-vcfs"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop through each .vcf.gz file in the input directory
for file in "$input_dir"/*.vcf; do
    # Get the base name of the file (without path)
    base_name=$(basename "$file" .vcf)
    
    # Define output file name in the output directory
    output_file="$output_dir/${base_name}_filtered.vcf.gz"
    
    # Run bcftools view to filter based on AF column
    bcftools view -i 'FILTER="PASS"'  "$file" -Oz -o "$output_file"
done

echo "Filtering complete. Filtered files are located in: $output_dir"
