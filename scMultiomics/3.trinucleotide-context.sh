#!/bin/bash

#PBS -l select=2:ncpus=64:mem=128g
#PBS -l walltime=10:00:00
#PBS -P 11003581
#PBS -j oe
#PBS -N run-scomatic

# ------------------------------------------------------------------------------
# SAFETY & LOGGING
# ------------------------------------------------------------------------------
set -euo pipefail
cd "$PBS_O_WORKDIR"

LOGFILE="$PBS_O_WORKDIR/scomatic-run.log"
exec > >(tee -a "$LOGFILE") 2>&1

echo "[INFO] SComatic pipeline started at $(date)"

# ------------------------------------------------------------------------------
# ENVIRONMENT SETUP
# ------------------------------------------------------------------------------
echo "[INFO] Activating environment"
source /home/project/11003581/Tools/miniforge3/bin/activate
conda activate /home/project/11003581/conda-envs/scMultiomics

module load gcc
module load intel/2024.0
module load libfabric/1.11.0.4.125
module load hdf5/1.12.1-parallel-icc23-cmpi
module load gsl/2.7.1-gcc11
module load r/4.2.0
module load samtools/1.15.1
module load bedtools/2.30.0
module load python/3.12.1-gcc11

# ------------------------------------------------------------------------------
# PATHS AND PARAMETERS
# ------------------------------------------------------------------------------
project="bc"
SCOMATIC="/home/project/11003581/Tools/SComatic"

base_dir="/home/users/nus/ash.ps/scratch/mulitomics/10x_data/data/breast-cancer"
output_dir="/home/users/nus/ash.ps/scratch/mulitomics/10x_data/analysis/"
output_dir2="$output_dir/Step2_BaseCellCounts"
output_dir3=$output_dir/Step3_BaseCellCountsMerged


mkdir -p $output_dir3

# ------------------------------------------------------------------------------
# ABsecellcount-merge
# ------------------------------------------------------------------------------

python $SCOMATIC/scripts/MergeCounts/MergeBaseCellCounts.py --tsv_folder ${output_dir2} \
  --outfile ${output_dir3}/${project}.BaseCellCounts.AllCellTypes.tsv

# ------------------------------------------------------------------------------
# TRINUCLEOTIDE CONTEXT
# ------------------------------------------------------------------------------

output_dir8=$output_dir/TrinucleotideContext
output_dir4=$output_dir/Step4_VariantCalling
mkdir -p $output_dir8

# List all available step1 files into the input list
ls ${output_dir4}/bc.*_filtered.calling.step1.tsv > ${output_dir8}/step1_files.txt

# Confirm the files have been listed
echo "[INFO] Using the following input files for trinucleotide context:"
cat ${output_dir8}/step1_files.txt

# Run the trinucleotide context background script with all files
python $SCOMATIC/scripts/TrinucleotideBackground/TrinucleotideContextBackground.py \
    --in_tsv ${output_dir8}/step1_files.txt \
    --out_file ${output_dir8}/TrinucleotideBackground.txt

