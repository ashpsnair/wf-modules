#!/bin/bash

#PBS -N auto-mut-sig
#PBS -l select=1:ncpus=128:mem=128gb
#PBS -l walltime=04:00:00
#PBS -j oe

# Change to the directory where the job was submitted
cd $PBS_O_WORKDIR

# Load necessary modules
module load miniforge3
conda activate /home/project/11003581/conda-envs/sigprofile

# Set variables
base_dir="/path/to/your/base/directory"
sample_name="your_sample_name"

# Run the combined analysis script
python combined_analysis.py $base_dir $sample_name

# Create result directory
result_dir="${sample_name}_mut_sig_results"
mkdir -p "$result_dir"

# Create and populate sigprofile subfolder
mkdir -p "$result_dir/sigprofile"
cp "$base_dir/sigpro_extract/SBS96/All_solutions_stat.csv" "$result_dir/sigprofile/"
cp "$base_dir/sigpro_extract/SBS96/SBS96_selection_plot.pdf" "$result_dir/sigprofile/"
cp -r "$base_dir/sigpro_extract/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Activities" "$result_dir/sigprofile/"
cp -r "$base_dir/sigpro_extract/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Signatures" "$result_dir/sigprofile/"
cp -r "$base_dir/sigpro_extract/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Solution_Stats" "$result_dir/sigprofile/"

# Create and populate sigminer subfolder
mkdir -p "$result_dir/sigminer"
cp "$base_dir/sigminer/similarity_heatmap.png" "$result_dir/sigminer/"
cp -r "$base_dir/sigminer/Decompose_Solution/Activities" "$result_dir/sigminer/"
cp -r "$base_dir/sigminer/Decompose_Solution/Signatures" "$result_dir/sigminer/"
cp -r "$base_dir/sigminer/Decompose_Solution/Solution_Stats" "$result_dir/sigminer/"

echo "Analysis complete. Results are in $result_dir"

# End of script
