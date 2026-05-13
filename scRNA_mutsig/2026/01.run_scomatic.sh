#!/bin/bash

#PBS -l select=2:ncpus=64:mem=128g
#PBS -l walltime=10:00:00
#PBS -P 11003581
#PBS -j oe
#PBS -N scomatic-clean

# ------------------------------------------------------------------------------
# SAFETY & LOGGING
# ------------------------------------------------------------------------------
set -euo pipefail
cd "$PBS_O_WORKDIR"

LOGFILE="$PBS_O_WORKDIR/scomatic-clean.log"
exec > >(tee -a "$LOGFILE") 2>&1

echo "[INFO] Starting clean SComatic run at $(date)"

# ------------------------------------------------------------------------------
# ENVIRONMENT (KEEP EXACTLY AS YOUR WORKING SETUP)
# ------------------------------------------------------------------------------
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
# PATHS
# ------------------------------------------------------------------------------
SCOMATIC="/home/project/11003581/Tools/SComatic"

base_dir="/home/users/nus/ash.ps/scratch/scRNA-mutsig/data"
output_dir="/home/users/nus/ash.ps/scratch/scRNA-mutsig/analysis"

output_dir2="$output_dir/Step2_BaseCellCounts"
output_dir4="$output_dir/Step4_VariantCalling"

REF="/home/project/11003581/Ref/Homo_sapiens/GATK/hg38/Homo_sapiens_assembly38.fasta"
editing="$SCOMATIC/RNAediting/AllEditingSites.hg38.txt"

threads=128

mkdir -p "$output_dir2" "$output_dir4"

# ------------------------------------------------------------------------------
# INPUT BAM
# ------------------------------------------------------------------------------
bam_file="$base_dir/Breast_Cancer_3p_possorted_genome_bam.bam"

echo "[INFO] Using BAM: $bam_file"

# ------------------------------------------------------------------------------
# STEP 1 — BaseCellCounter
# ------------------------------------------------------------------------------
echo "[STEP 1] BaseCellCounter"

python "$SCOMATIC/scripts/BaseCellCounter/BaseCellCounter.py" \
    --bam "$bam_file" \
    --ref "$REF" \
    --chrom all \
    --out_folder "$output_dir2" \
    --min_bq 20 \
    --min_ac 1 \
    --min_af 0.01 \
    --min_dp 2 \
    --min_cc 1 \
    --nprocs "$threads"

echo "[STEP 1 DONE]"

# ------------------------------------------------------------------------------
# STEP 2 — Variant calling step1
# ------------------------------------------------------------------------------
echo "[STEP 2] Step1 variant calling"

python "$SCOMATIC/scripts/BaseCellCalling/BaseCellCalling.step1.py" \
    --infile "$output_dir2/"*.tsv \
    --outfile "$output_dir4/merged" \
    --ref "$REF"

echo "[STEP 2 DONE]"

# ------------------------------------------------------------------------------
# STEP 3 — Variant calling step2 (OPTIONAL — may be empty)
# ------------------------------------------------------------------------------
echo "[STEP 3] Step2 variant calling (light filter)"

python "$SCOMATIC/scripts/BaseCellCalling/BaseCellCalling.step2.py" \
    --infile "$output_dir4/merged.calling.step1.tsv" \
    --outfile "$output_dir4/merged" \
    --editing "$editing"

echo "[DONE] SComatic run completed at $(date)"