#!/bin/bash

#PBS -l select=3:ncpus=64:mem=256g
#PBS -l walltime=10:00:00
#PBS -P 11003581
#PBS -N auto-annovar-multiomics
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR
module load bcftools/1.15.1

base_dir="/home/users/nus/ash.ps/scratch/mulitomics/analysis/"

###############################################
### annotation
###############################################

# Define the input and output directories
input_dir="${base_dir}/Step4_VariantCalling/"
output_dir="${base_dir}/annotated/"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

for vcf_file in "$input_dir"/*_filtered.calling.step1.tsv; do
  # Extract the filename without the extension
  samplename=$(basename "$vcf_file" _filtered.calling.step1.tsv)

  # Print the sample name
  echo "Processing sample: $samplename"

  # Create a corresponding output directory
  mkdir -p "$output_dir/$samplename"
  
  # Temporary VCF path
  vcf_compatible="$output_dir/$samplename/$samplename.vcf"

  # Convert SComatic TSV to minimal, clean VCF
  {
    echo "##fileformat=VCFv4.2"
    echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    awk -F'\t' 'BEGIN {OFS="\t"}
      !/^#/ && $4 ~ /^[ACGT]+$/ && $5 ~ /^[ACGT]+$/ {
        print $1, $2, ".", $4, $5, ".", "PASS", "."
      }' "$vcf_file"
  } > "$vcf_compatible"

  # Run the table_annovar.pl script with the converted VCF
  perl /home/project/11003581/Tools/annovar/table_annovar.pl "$vcf_compatible" \
      /home/project/11003581/Tools/annovar/humandb/ \
      -buildver hg38 \
      -out "$output_dir/$samplename/$samplename" \
      -protocol refGene,cosmic70,exac03,gnomad_genome,ALL.sites.2015_08,SAS.sites.2015_08,EAS.sites.2015_08,esp6500siv2_all,avsnp151,icgc28,dbnsfp30a,intervar_20180118,clinvar_20240611 \
      -operation g,f,f,f,f,f,f,f,f,f,f,f,f \
      -remove -vcfinput -polish -nastring .

done

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
### create maf files
###############################################

module load r/4.2.0

Rscript -e "
library(maftools)

# Set input and output directories
input_dir <- '/home/users/nus/ash.ps/scratch/mulitomics/analysis/pop-filter-multianno/'
output_dir <- '/home/users/nus/ash.ps/scratch/mulitomics/analysis/mafs/'

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Helper function
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
    cat('Created', output_name, 'MAF with', length(files), 'files\\n')
  } else {
    cat('No files found for', output_name, '\\n')
  }
}

# All files
all_files <- list.files(path = input_dir, pattern = '\\\\.hg38_multianno_pop_filt\\\\.txt$', full.names = TRUE)

# Group-specific filters
high_maf_files <- all_files[grepl('C10|C11|C12', all_files)]
intermediate_maf_files <- all_files[grepl('C9|C8|C7', all_files)]
low_maf_files <- all_files[grepl('C4|C5|C6', all_files)]
null_maf_files <- all_files[grepl('C1|C2|C3', all_files)]

# Generate MAFs
create_maf(high_maf_files, 'high_maf')
create_maf(intermediate_maf_files, 'intermediate_maf')
create_maf(low_maf_files, 'low_maf')
create_maf(null_maf_files, 'null_maf')

cat('All MAF files generated.\\n')
"
