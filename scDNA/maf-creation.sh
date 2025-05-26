#!/bin/bash

#PBS -l select=3:ncpus=64:mem=256g
#PBS -l walltime=04:00:00
#PBS -P 11003581
#PBS -N scdna-mafs
#PBS -j oe
#PBS -M ash.ps@nus.edu.sg

cd $PBS_O_WORKDIR

module load r/4.2.0

# Run the R command directly within the bash script using Rscript
Rscript -e "
library(maftools)

# Set input and output directories
input_dir <- '/home/users/nus/ash.ps/scratch/scDNA/analysis/annotation/pop-filter-multianno/3401'
output_dir <- '/home/users/nus/ash.ps/scratch/scDNA/analysis/annotation/mafs'

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
all_files <- list.files(path = input_dir, pattern = '_filtered_pop_filt\\\\.txt$', full.names = TRUE)

# Create combined MAF
create_maf(all_files, '3401')

cat('Processing complete.\n')
"
