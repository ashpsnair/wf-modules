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
project="Example"
SCOMATIC="/home/project/11003581/Tools/SComatic"

base_dir="/home/users/nus/ash.ps/scratch/scRNA-mutsig/data"
output_dir="/home/users/nus/ash.ps/scratch/scRNA-mutsig/analysis"
meta_celltype="/home/users/nus/ash.ps/scratch/scRNA-mutsig/data/cell_annotation_AMPK.txt"

output_dir1="$output_dir/Step1_BamCellTypes"
filtered_out="$output_dir1/filtered"
output_dir2="$output_dir/Step2_BaseCellCounts"
output_dir4="$output_dir/Step4_VariantCalling"
REF="/home/project/11003581/Ref/Homo_sapiens/GATK/hg38/Homo_sapiens_assembly38.fasta"
editing="$SCOMATIC/RNAediting/AllEditingSites.hg38.txt"
PON="$SCOMATIC/PoNs/PoN.scRNAseq.hg38.tsv"
threads=128

mkdir -p "$output_dir1" "$filtered_out" "$output_dir2" "$output_dir4"

bam_file=$(find "$base_dir" -name "*.bam" -type f)

# ------------------------------------------------------------------------------
# STEP 1: SPLIT BAM BY CELL TYPE
# ------------------------------------------------------------------------------
echo "[STEP 1.1] Splitting BAM by cell type at $(date)"
python "$SCOMATIC/scripts/SplitBam/SplitBamCellTypes.py" \
    --bam "$bam_file" \
    --meta "$meta_celltype" \
    --id "$project" \
    --n_trim 5 \
    --max_nM 5 \
    --max_NH 1 \
    --outdir "$output_dir1"

echo "[STEP 1.1] Done"

# ------------------------------------------------------------------------------
# STEP 1.5: FILTER & INDEX SPLIT BAM FILES
# ------------------------------------------------------------------------------
chromosomes="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"

echo "[STEP 1.2] Filtering and indexing BAMs at $(date)"
for input_bam in "$output_dir1"/*.bam; do
    base_name=$(basename "$input_bam" .bam)
    output_bam="${filtered_out}/${base_name}_filtered.bam"

    echo "[INFO] Filtering $base_name"
    samtools view -h "$input_bam" $chromosomes | samtools sort -@ "$threads" -o "$output_bam"
    samtools index "$output_bam"
    echo "[INFO] Indexed $output_bam"
done

echo "[STEP 1.2] Done"

echo "[STEP 1] Done"
