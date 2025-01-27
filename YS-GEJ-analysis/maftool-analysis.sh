##################### Creating MAF files from annovar outputs ##########################

#!/bin/bash

#PBS -l select=1:ncpus=64:mem=128g
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
input_dir <- '/home/users/nus/ash.ps/scratch/YS-analysis/pop-filter-multianno/'
output_dir <- '/home/users/nus/ash.ps/scratch/YS-analysis/mafs/'

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





#!/bin/bash

#PBS -l select=1:ncpus=64:mem=128g
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
input_dir <- '/home/users/nus/ash.ps/scratch/YS-analysis/pop-filter-multianno/'
output_dir <- '/home/users/nus/ash.ps/scratch/YS-analysis/mafs/'

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
w014_files <- all_files[!grepl('T14', all_files)]
create_maf(w014_files, 'combined_w014')

# Create high_maf
high_maf_files <- all_files[grepl('T11|T04|T10|T07', all_files)]
create_maf(high_maf_files, 'high_maf')

# Create intermediate_maf
intermediate_maf_files <- all_files[grepl('T13|T14|T08|T02|T09|T05', all_files)]
create_maf(intermediate_maf_files, 'intermediate_maf')

# Create low_maf
low_maf_files <- all_files[grepl('T06|T03|T01|T12', all_files)]
create_maf(low_maf_files, 'low_maf')

cat('Processing complete.\n')
"
