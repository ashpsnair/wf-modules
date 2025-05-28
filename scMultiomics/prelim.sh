#!/bin/bash

#PBS -l select=2:ncpus=64:mem=128g
#PBS -l walltime=10:00:00
#PBS -P 11003581
#PBS -j oe
#PBS -N run-scMulti

# ------------------------------------------------------------------------------
# SAFETY AND LOGGING SETUP
# ------------------------------------------------------------------------------
set -euo pipefail

cd "$PBS_O_WORKDIR"

# Define a custom absolute log file
LOGFILE="$PBS_O_WORKDIR/scomatic-run.log"
exec > >(tee -a "$LOGFILE") 2>&1

echo "[INFO] Logging started at $(date)"
echo "[INFO] Job started in: $PBS_O_WORKDIR"

# ------------------------------------------------------------------------------
# ENVIRONMENT SETUP
# ------------------------------------------------------------------------------
echo "[INFO] Activating conda environment: scMultiomics"
source /home/project/11003581/Tools/miniforge3/bin/activate
conda activate /home/project/11003581/conda-envs/scMultiomics

echo "[INFO] Loading required modules"
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
# PATH CONFIGURATION
# ------------------------------------------------------------------------------
project="bc"
base_dir="/home/users/nus/ash.ps/scratch/mulitomics/data/breast-cancer"
output_dir="/home/users/nus/ash.ps/scratch/mulitomics/analysis"
SCOMATIC="/home/project/11003581/Tools/SComatic"

echo "[INFO] Locating input files"
h5_file=$(find "$base_dir" -name "*.h5" -type f)
bam_file=$(find "$base_dir" -name "*.bam" -type f)
bai_file=$(find "$base_dir" -name "*.bam.bai" -type f)

# ------------------------------------------------------------------------------
# STEP 1: SPLIT BAM BY CELL TYPE
# ------------------------------------------------------------------------------
echo "[STEP 1] Splitting BAM by cell type"
output_dir1="$output_dir/Step1_BamCellTypes"
mkdir -p "$output_dir1"

python $SCOMATIC/scripts/SplitBam/SplitBamCellTypes.py \
    --bam "$bam_file" \
    --meta /home/users/nus/ash.ps/scratch/mulitomics/data/cell_type.tsv \
    --id "$project" \
    --n_trim 5 \
    --max_nM 5 \
    --max_NH 1 \
    --outdir "$output_dir1"

echo "[STEP 1] BAM splitting completed at $(date)"

# ------------------------------------------------------------------------------
# STEP 1.5: FILTER AND INDEX SPLIT BAM FILES (MULTITHREADED)
# ------------------------------------------------------------------------------
echo "[STEP 1.5] Filtering and indexing BAM files at $(date)"

filtered_out="$output_dir1/filtered"
mkdir -p "$filtered_out"

chromosomes="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"
threads=8

for input_bam in "$output_dir1"/*.bam; do
    base_name=$(basename "$input_bam")
    output_bam="${filtered_out}/${base_name%.bam}_filtered.bam"

    echo "[INFO] Filtering $base_name"
    samtools view -h "$input_bam" $chromosomes | samtools sort -@ "$threads" -o "$output_bam"

    echo "[INFO] Indexing $output_bam"
    samtools index "$output_bam"

    echo "[INFO] Completed: $base_name"
done

echo "[STEP 1.5] All BAM files filtered and indexed at $(date)"

# ------------------------------------------------------------------------------
# STEP 2–3: BASE COUNTING AND VARIANT CALLING
# ------------------------------------------------------------------------------
REF="/home/project/11003581/Ref/Homo_sapiens/GATK/hg38/Homo_sapiens_assembly38.fasta"
editing="$SCOMATIC/RNAediting/AllEditingSites.hg38.txt"
PON="$SCOMATIC/PoNs/PoN.scRNAseq.hg38.tsv"

output_dir2="$output_dir/Step2_BaseCellCounts"
output_dir4="$output_dir/Step4_VariantCalling"

mkdir -p "$output_dir2" "$output_dir4"

for bam in "$filtered_out"/*.bam; do
    cell_type=$(basename "$bam" | awk -F'.' '{print $(NF-1)}')
    temp="$output_dir2/temp_${cell_type}"
    mkdir -p "$temp"

    echo "[STEP 2] Running BaseCellCounter for $cell_type at $(date)"
    python $SCOMATIC/scripts/BaseCellCounter/BaseCellCounter.py \
        --bam "$bam" \
        --ref "$REF" \
        --chrom all \
        --out_folder "$output_dir2" \
        --min_bq 30 \
        --tmp_dir "$temp" \
        --nprocs 128

    echo "[STEP 3.1] Variant calling step 1 for $cell_type"
    python $SCOMATIC/scripts/BaseCellCalling/BaseCellCalling.step1.py \
        --infile "${output_dir2}/${cell_type}.BaseCellCounts.tsv" \
        --outfile "${output_dir4}/${cell_type}" \
        --ref "$REF"

    echo "[STEP 3.2] Variant calling step 2 for $cell_type"
    python $SCOMATIC/scripts/BaseCellCalling/BaseCellCalling.step2.py \
        --infile "${output_dir4}/${cell_type}.calling.step1.tsv" \
        --outfile "${output_dir4}/${cell_type}" \
        --editing "$editing" \
        --pon "$PON"

    echo "[INFO] Cleaning temp folder for $cell_type"
    rm -rf "$temp"
done

# ------------------------------------------------------------------------------
# COMPLETE
# ------------------------------------------------------------------------------
echo "[COMPLETE] SComatic workflow finished at $(date)"
exit 0
