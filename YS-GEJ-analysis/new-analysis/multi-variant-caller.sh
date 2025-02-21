################### Merge strelka snvs and indels 
#!/bin/bash

#PBS -l select=1:ncpus=64
#PBS -l walltime=1:00:00
#PBS -P 11003581
#PBS -N intersect
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load bcftools/1.15.1

# Set input and output directories
SNV_DIR="/home/users/nus/ash.ps/scratch/YS-GEJ/FINAL/strelka/snvs/"
INDEL_DIR="/home/users/nus/ash.ps/scratch/YS-GEJ/FINAL/strelka/indel/"
OUTPUT_DIR="/home/users/nus/ash.ps/scratch/YS-GEJ/FINAL/strelka/merged/"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through all samples
for snv_file in "$SNV_DIR"*.strelka.somatic_snvs.vcf.gz; do
    # Extract sample name
    sample_name=$(basename "$snv_file" .strelka.somatic_snvs.vcf.gz)
    
    # Find corresponding indel file
    indel_file="$INDEL_DIR${sample_name}.strelka.somatic_indels.vcf.gz"
    
    # Check if indel file exists
    if [ -f "$indel_file" ]; then
        # Merge SNV and indel files
        bcftools merge -o "${OUTPUT_DIR}${sample_name}.merged.vcf.gz" -O z "$snv_file" "$indel_file"
        
        # Index the merged file
        bcftools index "${OUTPUT_DIR}${sample_name}.merged.vcf.gz"
        
        echo "Merged and indexed: ${sample_name}"
    else
        echo "Warning: No matching indel file found for ${sample_name}"
    fi
done

echo "Merging complete. Results stored in ${OUTPUT_DIR}"


################### intersect strelka and mutect2
#!/bin/bash

#PBS -l select=1:ncpus=64
#PBS -l walltime=1:00:00
#PBS -P 11003581
#PBS -N intersect
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load bcftools/1.15.1

# Define directories
mutect2_dir="/home/users/nus/ash.ps/scratch/YS-GEJ/FINAL/mutect2/"
strelka_dir="/home/users/nus/ash.ps/scratch/YS-GEJ/FINAL/strelka/merged/"
output_dir="/home/users/nus/ash.ps/scratch/YS-GEJ/FINAL/common-snvs"
temp_dir="/tmp/vcf_processing"

# Create output and temp directories if they don't exist
mkdir -p "$output_dir" "$temp_dir"

# Process each sample
for mutect2_file in "$mutect2_dir"/*.mutect2.filtered.vcf.gz; do
    # Extract sample name
    sample_name=$(basename "$mutect2_file" .mutect2.filtered.vcf.gz)
    
    # Construct Strelka file name
    strelka_file="$strelka_dir/${sample_name}.merged.vcf.gz"
    
    # Check if Strelka file exists
    if [ ! -f "$strelka_file" ]; then
        echo "Strelka file not found for sample $sample_name. Skipping..."
        continue
    fi
    
    echo "Processing sample: $sample_name"
    
    # 1. Apply PASS filter
    bcftools view -f PASS "$mutect2_file" -Oz -o "$temp_dir/${sample_name}.mutect2_pass.vcf.gz"
    bcftools view -f PASS "$strelka_file" -Oz -o "$temp_dir/${sample_name}.strelka_pass.vcf.gz"
    
    # 2. Index filtered files
    bcftools index "$temp_dir/${sample_name}.mutect2_pass.vcf.gz"
    bcftools index "$temp_dir/${sample_name}.strelka_pass.vcf.gz"
    
    # 3. Perform intersection
    bcftools isec -p "$output_dir/${sample_name}_intersection" -Oz \
        "$temp_dir/${sample_name}.mutect2_pass.vcf.gz" \
        "$temp_dir/${sample_name}.strelka_pass.vcf.gz"
    
    # Clean up temporary files
    rm "$temp_dir/${sample_name}.mutect2_pass.vcf.gz"*
    rm "$temp_dir/${sample_name}.strelka_pass.vcf.gz"*
    
    echo "Completed processing for sample: $sample_name"
done

# Remove temporary directory
rmdir "$temp_dir"

echo "All samples processed. Results are in $output_dir"




############

#!/bin/bash

#PBS -l select=1:ncpus=64
#PBS -l walltime=00:30:00
#PBS -P 11003581
#PBS -N intersect
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load bcftools/1.15.1
# Define directories
intersection_dir="/home/users/nus/ash.ps/scratch/YS-GEJ/FINAL/common-snvs"
csv_file="$output_dir/variant_counts.csv"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Create CSV file with header
echo "samplename,mutect2-uniq,strelka-uniq,common" > "$csv_file"

# Process each sample folder
for sample_dir in "$intersection_dir"/*_intersection; do
    # Extract sample name
    sample_name=$(basename "$sample_dir" _intersection)
    
    echo "Processing sample: $sample_name"
    
    # Copy and rename 0002.vcf.gz
    cp "$sample_dir/0002.vcf.gz" "$output_dir/${sample_name}_common.vcf.gz"
    
    # Count variants in each file
    mutect2_uniq=$(zcat "$sample_dir/0000.vcf.gz" | grep -v '^#' | wc -l)
    strelka_uniq=$(zcat "$sample_dir/0001.vcf.gz" | grep -v '^#' | wc -l)
    common=$(zcat "$sample_dir/0002.vcf.gz" | grep -v '^#' | wc -l)
    
    # Append to CSV file
    echo "$sample_name,$mutect2_uniq,$strelka_uniq,$common" >> "$csv_file"
    
    echo "Completed processing for sample: $sample_name"
done

echo "All samples processed. Results are in $output_dir"
echo "Variant counts are saved in $csv_file"
