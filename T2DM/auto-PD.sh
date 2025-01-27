#!/bin/bash

#PBS -l select=1:ncpus=128:mem=128g
#PBS -l walltime=05:00:00
#PBS -P 11003581
#PBS -N auto-process
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load bcftools/1.15.1

base_dir="/home/users/nus/ash.ps/scratch/T2DM/PD-analysis/"


###############################################
### pop filtration vcf files
###############################################

INPUT_DIR="${base_dir}/annotated/"
OUTPUT_DIR="${base_dir}/pop-filter-vcfs"

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


###############################################
### pop filtration on multianno files
###############################################

INPUT_DIR="${base_dir}/annotated/"
OUTPUT_DIR="${base_dir}/pop-filter-multianno"

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


###############################################
### pop filtration on multianno files
###############################################

module load r/4.2.0



# Run the R command directly within the bash script using Rscript
Rscript -e "
library(maftools)

# Set input and output directories
input_dir <- '/home/users/nus/ash.ps/scratch/YS-tumor-only/normal-analysis/pop-filter-multianno/'
output_dir <- '/home/users/nus/ash.ps/scratch/YS-tumor-only/normal-analysis/mafs/'

# Create output directory if it doesn't exist
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Function to process files and create MAF
create_maf <- function(files, output_name) {
  if (length(files) > 0) {
    multi_maf <- annovarToMaf(
      files,
      Center = NULL,
      refBuild = 'hg38',
      tsbCol = NULL,
      ens2hugo = TRUE,
      basename = NULL,
      sep = '\t',
      MAFobj = FALSE,
      sampleAnno = NULL
    )
    output_file <- file.path(output_dir, paste0(output_name, '.maf'))
    write.table(multi_maf, file = output_file, sep = '\t', quote = FALSE, row.names = FALSE)
    cat('Created', output_name, 'MAF with', length(files), 'files\n')
  } else {
    cat('No files found for', output_name, '\n')
  }
}

# Get list of all .txt files
all_files <- list.files(path = input_dir, pattern = '\\.hg38_multianno_pop_filt\\.txt$', full.names = TRUE)

# Create combined MAF
create_maf(all_files, 'combined')

# Create combined_w014 MAF (all except T14)
w014_files <- all_files[!grepl('N14', all_files)]
create_maf(w014_files, 'combined_w014')

# Create high_maf
high_maf_files <- all_files[grepl('N11|N04|N10|N07', all_files)]
create_maf(high_maf_files, 'high_maf')

# Create intermediate_maf
intermediate_maf_files <- all_files[grepl('N13|N14|N08|N02|N09|N05', all_files)]
create_maf(intermediate_maf_files, 'intermediate_maf')

# Create low_maf
low_maf_files <- all_files[grepl('N06|N03|N01|N12', all_files)]
create_maf(low_maf_files, 'low_maf')

cat('Processing complete.\n')
"
