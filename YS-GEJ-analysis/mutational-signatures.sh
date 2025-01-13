############# categorizing the files based on tumor ############

#!/bin/bash
#PBS -l select=1:ncpus=32
#PBS -l walltime=01:00:00
#PBS -P 11003581
#PBS -N categorize
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR


INPUT_OUTPUT_DIR="/home/users/nus/ash.ps/scratch/YS-analysis/filter-annovar"

# Create category subdirectories
mkdir -p "$INPUT_OUTPUT_DIR/Normal" "$INPUT_OUTPUT_DIR/Tumor" "$INPUT_OUTPUT_DIR/Tumor-only"

# Loop through all .vcf files in the directory
for vcf_file in "$INPUT_OUTPUT_DIR"/*.vcf; do
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


####################### Making folders of categories ############

#!/bin/bash
#PBS -l select=1:ncpus=32
#PBS -l walltime=01:00:00
#PBS -P 11003581
#PBS -N categorize
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load miniforge3
conda activate /home/project/11003581/conda-envs/sigprofile

import os
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen

input_dir = "/home/users/nus/ash.ps/scratch/YS-analysis/filter-annovar/Tumor"

# Iterate through each folder in the input directory
for category_folder in os.listdir(input_dir):
    category_path = os.path.join(input_dir, category_folder)
    
    # Check if it's a directory
    if os.path.isdir(category_path):
        print(f"Processing category: {category_folder}")
        
        # Run SigProfilerMatrixGeneratorFunc for this category
        matrices = matGen.SigProfilerMatrixGeneratorFunc(
            category_folder,  # Use folder name as Category_name
            "GRCh38",
            category_path,    # Use the full path to the category folder
            plot=True,
            exome=False,
            bed_file=None,
            chrom_based=False,
            tsb_stat=False,
            seqInfo=False,
            cushion=100
        )
        
        print(f"Completed processing for category: {category_folder}")

print("All categories processed.")












