#!/bin/bash

# Base directory
base_dir="/Users/ash/Downloads/MASH-vcfs/pooled-patientwise-vaf20"

# Loop through each subdirectory
for folder in "$base_dir"/*; do
    if [ -d "$folder" ]; then
        for subfolder in "$folder"/*; do
            if [ -d "$subfolder/output" ]; then
                output_dir="$subfolder/output"
                # Collect VCF files
                vcf_files=("$output_dir"/*.vcf)
                
                # Check if there are VCF files to process
                if [ ${#vcf_files[@]} -gt 0 ]; then
                    # Create a temporary array for bgzipped files
                    bgzipped_files=()

                    # Bgzip and index each VCF file
                    for vcf in "${vcf_files[@]}"; do
                        bgzip "$vcf" -c > "$vcf.gz"  # Compress the VCF file
                        bcftools index "$vcf.gz"      # Index the compressed file
                        bgzipped_files+=("$vcf.gz")   # Add to the array of bgzipped files
                    done

                    # Extract folder name for output file naming
                    folder_name=$(basename "$subfolder")
                    output_file="$folder_name.merged.vcf.gz"
                    
                    # Merge the indexed bgzipped VCF files
                    bcftools merge "${bgzipped_files[@]}" -Oz -o "$output_file" --force-samples
                    
                    echo "Processed and merged into $output_file"
                else
                    echo "No VCF files found in $output_dir"
                fi
            fi
        done
    fi
done