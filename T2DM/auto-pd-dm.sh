#!/bin/bash

#PBS -l select=1:ncpus=128:mem=128g
#PBS -l walltime=05:00:00
#PBS -P 11003581
#PBS -N auto-pd-dm
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load bcftools/1.15.1

base_dir="/home/users/nus/ash.ps/scratch/T2DM/PD-DM-analysis/"

# Set input and output directories
input_dir=${base_dir}/VCFs/
output_dir=${base_dir}/filtered-vcfs

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

###############################################
### annotation
###############################################

# Define the input and output directories
input_dir="${base_dir}/filtered-vcfs"
output_dir="${base_dir}/annotated/"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

for vcf_file in "$input_dir"/*_filtered.vcf; do
  # Extract the filename without the '.mutect2.vcf' extension
  samplename=$(basename "$vcf_file" _filtered.vcf)

  # Print the sample name
  echo "Processing sample: $samplename"

  # Create a corresponding output directory
  mkdir -p "$output_dir/$samplename"

  # Run the table_annovar.pl script with the appropriate parameters
  perl /home/project/11003581/Tools/annovar/table_annovar.pl "$vcf_file" \
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
### creating maf files
###############################################

module load r/4.2.0

Rscript -e "
library(maftools)

# Set input and output directories
input_dir <- '/home/users/nus/ash.ps/scratch/T2DM/PD-DM-analysis/pop-filter-multianno/'
output_dir <- '/home/users/nus/ash.ps/scratch/T2DM/PD-DM-analysis/mafs/'

# Create output directory if it doesn't exist
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Get list of .txt files in the input directory
annovar_outputs <- list.files(path = input_dir, pattern = '\\\\.hg38_multianno_pop_filt\\\\.txt$', full.names = TRUE)

if (length(annovar_outputs) > 0) {
  # Run the annovarToMaf function
  multi_maf <- annovarToMaf(
    annovar_outputs,
    Center = NULL,
    refBuild = 'hg38',
    tsbCol = NULL,
    ens2hugo = TRUE,
    basename = NULL,
    sep = '\t',
    MAFobj = FALSE,
    sampleAnno = NULL
  )
  
  # Create output filename
  output_file <- file.path(output_dir, 'combined_tumor.maf')
  
  # Write the MAF to file
  write.table(multi_maf, file = output_file, sep = '\t', quote = FALSE, row.names = FALSE)
  
  cat('Processed', length(annovar_outputs), 'files and saved combined MAF to', output_file, '\n')
} else {
  cat('No .hg38_multianno_pop_filt.txt files found in', input_dir, '\n')
}

cat('Processing complete.\n')
"
