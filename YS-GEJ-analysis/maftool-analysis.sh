

######### FIltering the multianno files ##########

#!/bin/bash
#PBS -l select=1:ncpus=32:mem=64g
#PBS -l walltime=01:00:00
#PBS -P 11003581
#PBS -N run-filter-multianno
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load python/3.12.1-gcc11

cat << EOF > filter-multianno.py
import os
import glob
import pandas as pd

INPUT_DIR = "/home/users/nus/ash.ps/scratch/YS-analysis/annovar-out"
OUTPUT_DIR = "/home/users/nus/ash.ps/scratch/YS-analysis/multianno-filtered"

# Create output directory if it doesn't exist
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Function to apply filtering criteria
def filter_row(row):
    conditions = [
        row['gnomad41_genome_AF'] == '.' or float(row['gnomad41_genome_AF']) <= 0.01,
        row['gnomad41_genome_AF_eas'] == '.' or float(row['gnomad41_genome_AF_eas']) <= 0.01,
        row['gnomad41_genome_AF_sas'] == '.' or float(row['gnomad41_genome_AF_sas']) <= 0.01,
        row['ExAC_ALL'] == '.' or float(row['ExAC_ALL']) <= 0.01,
        row['ExAC_SAS'] == '.' or float(row['ExAC_SAS']) <= 0.01,
        row['esp6500siv2_all'] == '.' or float(row['esp6500siv2_all']) <= 0.01
    ]
       
    return all(conditions)

# Process each .hg38_multianno.txt file
for txt_file in glob.glob(os.path.join(INPUT_DIR, "**", "*.hg38_multianno.txt"), recursive=True):
    print(f"Processing: {txt_file}")
    
    # Generate output filename
    output_file = os.path.join(OUTPUT_DIR, os.path.basename(txt_file).replace('.txt', '_filtered.txt'))
    
    # Read the file
    df = pd.read_csv(txt_file, sep='\t', low_memory=False)
    
    # Apply the filter
    filtered_df = df[df.apply(filter_row, axis=1)]
    
    # Save the filtered data
    filtered_df.to_csv(output_file, sep='\t', index=False)
    
    print(f"Saved filtered file: {output_file}")

print("All files processed.")

EOF

# Run the Python script
python filter-multianno.py

# Clean up the temporary Python script
rm filter-multianno.py


#####################

############# categorizing the files based on tumor ############

#!/bin/bash
#PBS -l select=1:ncpus=32
#PBS -l walltime=01:00:00
#PBS -P 11003581
#PBS -N categorize
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR


INPUT_OUTPUT_DIR="/home/users/nus/ash.ps/scratch/YS-analysis/multianno-filtered"

# Create category subdirectories
mkdir -p "$INPUT_OUTPUT_DIR/Normal" "$INPUT_OUTPUT_DIR/Tumor" "$INPUT_OUTPUT_DIR/Tumor-only"

# Loop through all .vcf files in the directory
for vcf_file in "$INPUT_OUTPUT_DIR"/*.txt; do
    filename=$(basename "$vcf_file")

    if [[ $filename == *"_vs_"* ]]; then
        mv "$vcf_file" "$INPUT_OUTPUT_DIR/Tumor/"
    elif [[ $filename == T* ]]; then
        mv "$vcf_file" "$INPUT_OUTPUT_DIR/Tumor-only/"
    elif [[ $filename == N* ]]; then
        mv "$vcf_file" "$INPUT_OUTPUT_DIR/Normal/"
    else
        echo "Unknown category for file: $filename"
    fi
done

echo "Categorization complete."



##################### Running MAFTOOLS ##########################

#!/bin/bash

#PBS -l select=1:ncpus=32:mem=128g
#PBS -l walltime=03:30:00
#PBS -P 11003581
#PBS -N run-create-maf
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load r/4.2.0

# Run the R command directly within the bash script using Rscript
Rscript -e "
library(maftools)

# Set input and output directories
input_dir <- '/home/users/nus/ash.ps/scratch/YS-analysis/multianno-filtered/Tumor'
output_dir <- '/home/users/nus/ash.ps/scratch/YS-analysis/mafs/Tumor'

# Create output directory if it doesn't exist
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Get list of subdirectories
subdirs <- list.dirs(input_dir, full.names = TRUE, recursive = FALSE)

# Process each subdirectory
for (subdir in subdirs) {
  # Get the subfolder name
  subfolder_name <- basename(subdir)
  
  # Get the list of .hg38_multianno.txt files in the current subdirectory
  annovar_outputs <- list.files(path = subdir, pattern = '\\\\.hg38_multianno_filtered\\\\.txt$', full.names = TRUE)
  
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
    output_file <- file.path(output_dir, paste0(subfolder_name, '_tumor.maf'))
    
    # Write the MAF to file
    write.table(multi_maf, file = output_file, sep = '\t', quote = FALSE, row.names = FALSE)
    
    cat('Processed', subfolder_name, 'and saved MAF to', output_file, '\n')
  } else {
    cat('No .hg38_multianno.txt files found in', subfolder_name, '\n')
  }
}

cat('Processing complete.\n')



"