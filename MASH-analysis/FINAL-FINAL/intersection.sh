#!/bin/bash
#PBS -l select=1:ncpus=128:mem=128g
#PBS -l walltime=05:00:00
#PBS -P 11003581
#PBS -N intersect-mash
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load bcftools/1.15.1
# Set the base directory for raw files and output directory
base_dir="/home/users/nus/ash.ps/scratch/NCCS-MASH/FINAL/raw_files"
output_dir="/home/users/nus/ash.ps/scratch/NCCS-MASH/FINAL/intersect"

# Path to bgzip and tabix
BGZIP="/home/project/11003581/Tools/bin/bgzip"
TABIX="/home/project/11003581/Tools/bin/tabix"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop through every patient directory inside the base directory
for patient_dir in "$base_dir"/*/; do
    # Extract the patient name from the directory name
    patient_name=$(basename "$patient_dir")
    
    # Define file names for the current patient
    haplotype_vcf="${patient_dir}${patient_name}.haplotypecaller.filtered.vcf.gz"
    strelka_vcf="${patient_dir}${patient_name}.strelka.variants.vcf.gz"
    normal_vcf="${patient_dir}normal-${patient_name}.mutect2.filtered.vcf.gz"
    
    # Check if all required files are present in the directory
    if [[ -f "$haplotype_vcf" && -f "$strelka_vcf" && -f "$normal_vcf" ]]; then
        echo "Processing patient: $patient_name"
        
        # Step 1: Merge germline variant files (haplotypecaller and strelka)
        echo "Merging germline variant files for $patient_name..."
        bcftools merge "$haplotype_vcf" "$strelka_vcf" --force-samples -o "${output_dir}/${patient_name}.germline_merged.vcf.gz" -O z

        # Check if merge was successful
        if [[ ! -f "${output_dir}/${patient_name}.germline_merged.vcf.gz" ]]; then
            echo "Error: Failed to merge germline files for $patient_name. Skipping this patient."
            continue
        fi

        # Step 2: Index the merged VCF file
        echo "Indexing merged germline VCF for $patient_name..."
        "$TABIX" -p vcf "${output_dir}/${patient_name}.germline_merged.vcf.gz"

        # Step 3: Subtract germline variants from normal tissue VCF
        echo "Subtracting germline variants from normal tissue VCF for $patient_name..."
        bcftools isec -n-1 -c all "${output_dir}/${patient_name}.germline_merged.vcf.gz" "$normal_vcf" -p "${output_dir}/${patient_name}_subtracted"

        # Check if the subtracted file exists
        if [[ ! -f "${output_dir}/${patient_name}_subtracted/0000.vcf" ]]; then
            echo "Error: Subtracted VCF file not found for $patient_name. Skipping this patient."
            continue
        fi

        # Step 4: Rename the resulting VCF to "true-normal-{patient_name}.vcf"
        true_normal_vcf="${output_dir}/true-normal-${patient_name}.vcf"
        mv "${output_dir}/${patient_name}_subtracted/0000.vcf" "$true_normal_vcf"
        rm -r "${output_dir}/${patient_name}_subtracted/"

        # Step 5: Compress and index the result using bgzip and tabix from specified path
        echo "Compressing and indexing the result for true-normal-${patient_name}..."
        "$BGZIP" "$true_normal_vcf"
        "$TABIX" -p vcf "${true_normal_vcf}.gz"

        echo "Processing for patient $patient_name is complete."
    else
        echo "One or more required files are missing for patient $patient_name. Skipping this patient."
    fi
done

echo "Germline variants subtracted from all patients. Results stored in $output_dir."
