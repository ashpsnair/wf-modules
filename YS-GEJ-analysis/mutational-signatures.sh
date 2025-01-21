'''

module load miniforge3
conda create --prefix /home/project/11003581/conda-envs/sigprofile python=3.9.12

conda activate /home/project/11003581/conda-envs/sigprofile
pip install SigProfilerExtractor

conda install r-base r-devtools r-reticulate -c conda-forge -y

python
from SigProfilerMatrixGenerator import install as genInstall
genInstall.install('GRCh38', rsync=True, bash=True)

However, it wasnt working fr=or me so i manually copied the grch38 form local system
copied required files to /home/project/11003581/conda-envs/sigprofile/lib/python3.9/site-packages/SigProfilerMatrixGenerator/references/chromosomes/tsb/

'''


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


####################### Generating matrix files for each of the categories ############

#!/bin/bash
#PBS -l select=1:ncpus=128:mem=128gb
#PBS -l walltime=01:00:00
#PBS -P 11003581
#PBS -N run-matrixgen
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load miniforge3
conda activate /home/project/11003581/conda-envs/sigprofile

# Create a temporary Python script
cat << EOF > run_sigprofiler.py
import os
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen

input_dir = "/home/users/nus/ash.ps/scratch/YS-analysis/multianno-filtered/"

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
EOF

# Run the Python script
python run_sigprofiler.py

# Clean up the temporary Python script
rm run_sigprofiler.py


################## Extracting signatures using sigprofile extractor #############

#!/bin/bash
#PBS -l select=1:ncpus=128:mem=128gb
#PBS -l walltime=10:00:00
#PBS -P 11003581
#PBS -N run-sigprofile-extractor
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load miniforge3
conda activate /home/project/11003581/conda-envs/sigprofile

# Change to the directory where the script is submitted from
cd $PBS_O_WORKDIR

# Create a temporary Python script
cat << EOF > run_sigprofile-extractor.py

import os
import glob
from SigProfilerExtractor import sigpro as sig

# Set the path to the directory containing your folders
input_dir = "/home/users/nus/ash.ps/scratch/YS-analysis/filter-annovar/Tumor"

# Get all immediate subdirectories
subfolders = [f.path for f in os.scandir(input_dir) if f.is_dir()]

for folder in subfolders:
    print(f"Processing folder: {folder}")
    
    # Change directory to the current folder
    os.chdir(folder)

    # Find the SBS96.all file
    sbs_file = glob.glob(os.path.join("output", "SBS", "*SBS96.all"))
    input_matrix= sbs_file[0]

    # Run SigProfilerExtractor
    sig.sigProfilerExtractor("matrix", "sigpro_extract", input_matrix, reference_genome="GRCh38", opportunity_genome="GRCh38", minimum_signatures=1, maximum_signatures=10, nmf_replicates=100)

    print(f"SigProfilerExtractor analysis complete for {folder}. Results are in: {output_dir}")

print("All folders processed.")

EOF

# Run the Python script
python run_sigprofile-extractor.py

# Clean up the temporary Python script
rm run_sigprofile-extractor.py

#######################################################

############# categorizing the files based on tumor ############

#!/bin/bash
#PBS -l select=1
#PBS -l walltime=00:05:00
#PBS -P 11003581
#PBS -N categorize
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR


INPUT_OUTPUT_DIR="/home/users/nus/ash.ps/scratch/YS-analysis/VCFs/"

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





############################################ TRIAL for automating sigprofile extractor ###############################
#!/bin/bash
#PBS -l select=1:ncpus=128:mem=128gb
#PBS -l walltime=10:00:00
#PBS -P 11003581
#PBS -N run-sigprofile-extractor
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load miniforge3
conda activate /home/project/11003581/conda-envs/sigprofile

# Change to the directory where the script is submitted from
cd $PBS_O_WORKDIR

# Create a temporary Python script
