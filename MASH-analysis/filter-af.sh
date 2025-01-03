#!/bin/bash

# Set input and output directories
input_dir="/Users/ash/Downloads/MASH-vcfs/MASH-vcfs"
output_dir="/Users/ash/Downloads/MASH-vcfs/filter-vaf30-MASH"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop through each .vcf.gz file in the input directory
for file in "$input_dir"/*.vcf.gz; do
    # Get the base name of the file (without path)
    base_name=$(basename "$file" .vcf.gz)
    
    # Define output file name in the output directory
    output_file="$output_dir/${base_name}_filtered.vcf.gz"
    
    # Run bcftools view to filter based on AF column
    bcftools view -i 'FORMAT/AF >= 0.3' "$file" -Oz -o "$output_file"
done

echo "Filtering complete. Filtered files are located in: $output_dir"
