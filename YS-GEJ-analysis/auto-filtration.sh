######################################################################################################
###################################  Pre-filtration on VCF files (before annotation)  ###############################################
######################################################################################################

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
input_dir="/home/users/nus/ash.ps/scratch/YS-analysis/VCFs/"
output_dir="/home/users/nus/ash.ps/scratch/YS-analysis/filtered-vcfs"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop through each .vcf.gz file in the input directory
for file in "$input_dir"/*.mutect2.filtered.vcf.gz; do
    # Get the base name of the file (without path)
    base_name=$(basename "$file" .mutect2.filtered.vcf.gz)

    # Define output file name in the output directory
    output_file="$output_dir/${base_name}_filtered.vcf"

    # Run bcftools view to filter based on AF column
    bcftools view -i 'FILTER="PASS"'  "$file" -o "$output_file"
done

echo "Filtering complete. Filtered files are located in: $output_dir"

######################################################################################################
###################################  pop Filtration vcfs  ###############################################
######################################################################################################

#!/bin/bash
#PBS -l select=1:ncpus=32
#PBS -l walltime=1:00:00
#PBS -P 11003581
#PBS -N pop-filter-vcfs
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load bcftools/1.15.1

INPUT_DIR="/home/users/nus/ash.ps/scratch/YS-analysis/annotated/"
OUTPUT_DIR="/home/users/nus/ash.ps/scratch/YS-analysis/pop-filter-vcfs"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Find all .hg38_multianno.vcf files in subfolders
find "$INPUT_DIR" -name "*.hg38_multianno.vcf" | while read -r vcf_file; do
    # Generate output filename
    output_file="$OUTPUT_DIR/$(basename "${vcf_file%.vcf}_pop_filt.vcf")"

    # Process each file
    bcftools view "$vcf_file" | awk '
    BEGIN {FS="\t"; OFS="\t"}
    /^#/ {print; next}
    {
        split($8, info, ";");
        exac_all="."; exac_eas="."; exac_sas="."; gnomad_genome_all="."; gnomad_genome_eas="."; all_sites_2015_08="."; esp6500siv2_all=".";
        for (i in info) {
            split(info[i], pair, "=");
            if (pair[1] == "ExAC_ALL") exac_all = pair[2];
            if (pair[1] == "ExAC_EAS") exac_eas = pair[2];
            if (pair[1] == "ExAC_SAS") exac_sas = pair[2];
            if (pair[1] == "gnomAD_genome_ALL") gnomad_genome_all = pair[2];
            if (pair[1] == "gnomAD_genome_EAS") gnomad_genome_eas = pair[2];
            if (pair[1] == "ALL.sites.2015_08") all_sites_2015_08 = pair[2];
            if (pair[1] == "esp6500siv2_all") esp6500siv2_all = pair[2];
        }
        if ((exac_all == "." || exac_all <= 0.01) && 
            (exac_eas == "." || exac_eas <= 0.01) && 
            (exac_sas == "." || exac_sas <= 0.01) && 
            (gnomad_genome_all == "." || gnomad_genome_all <= 0.01) && 
            (gnomad_genome_eas == "." || gnomad_genome_eas <= 0.01) && 
            (all_sites_2015_08 == "." || all_sites_2015_08 <= 0.01) && 
            (esp6500siv2_all == "." || esp6500siv2_all <= 0.01)) 
            print
    }' > "$output_file"

    echo "Processed: $vcf_file -> $output_file"
done


######################################################################################################
###################################  pop Filtration multianno files  ###############################################
######################################################################################################

#!/bin/bash
#PBS -l select=1:ncpus=32
#PBS -l walltime=1:00:00
#PBS -P 11003581
#PBS -N pop-filter-txt
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

INPUT_DIR="/home/users/nus/ash.ps/scratch/YS-analysis/annotated/"
OUTPUT_DIR="/home/users/nus/ash.ps/scratch/YS-analysis/pop-filter-multianno"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Find all .txt files in subfolders
find "$INPUT_DIR" -name "*.hg38_multianno.txt" | while read -r txt_file; do
    # Generate output filename
    output_file="$OUTPUT_DIR/$(basename "${txt_file%.txt}_pop_filt.txt")"

    # Process each file
    awk 'BEGIN {FS=OFS="\t"}
    NR==1 {print; next}  # Print header line
    {
        exac_all=$12; exac_eas=$15; exac_sas=$19; gnomad_genome_all=$20; gnomad_genome_eas=$24; all_sites_2015_08=$28; esp6500siv2_all=$29;
        if ((exac_all == "." || exac_all <= 0.01) && 
            (exac_eas == "." || exac_eas <= 0.01) && 
            (exac_sas == "." || exac_sas <= 0.01) && 
            (gnomad_genome_all == "." || gnomad_genome_all <= 0.01) && 
            (gnomad_genome_eas == "." || gnomad_genome_eas <= 0.01) && 
            (all_sites_2015_08 == "." || all_sites_2015_08 <= 0.01) && 
            (esp6500siv2_all == "." || esp6500siv2_all <= 0.01)) 
            print
    }' "$txt_file" > "$output_file"

    echo "Processed: $txt_file -> $output_file"
done

