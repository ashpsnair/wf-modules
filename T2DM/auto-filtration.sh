####### For PD


#!/bin/bash

#PBS -l select=1:ncpus=64
#PBS -l walltime=1:00:00
#PBS -P 11003581
#PBS -N filter-vcfs
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load bcftools/1.15.1

# Set input and output directories
input_dir="/home/users/nus/ash.ps/scratch/T2DM/PD-analysis/*/variant_calling/mutect2/*/"
output_dir="/home/users/nus/ash.ps/scratch/T2DM/PD-analysis/filtered-vcfs/"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop through each .vcf.gz file in the input directory
find $input_dir -name "*.mutect2.filtered.vcf.gz" | while read -r file; do
    # Get the base name of the file (without path)
    base_name=$(basename "$file" .mutect2.filtered.vcf.gz)

    # Define output file name in the output directory
    output_file="$output_dir/${base_name}_filtered.vcf"

    # Run bcftools view to filter based on AF column
    bcftools view -i 'FILTER="PASS"'  "$file" -o "$output_file"
done

echo "Filtering complete. Filtered files are located in: $output_dir"

