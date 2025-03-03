#!/bin/bash

#PBS -l select=2:ncpus=64:mem=128g
#PBS -l walltime=10:00:00
#PBS -P 11003581
#PBS -j oe
#PBS -N run-auto-scMulti


# Change to the directory where the job was submitted
cd $PBS_O_WORKDIR

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

# Set variables
project="breast-cancer"
base_dir="/home/users/nus/ash.ps/scratch/mulitomics/data/breast-cancer"
output_dir="/home/users/nus/ash.ps/scratch/mulitomics/analysis/breast-cancer"
SCOMATIC=/home/project/11003581/Tools/SComatic/
REF=/home/project/11003581/Ref/Homo_sapiens/GATK/hg38/Homo_sapiens_assembly38.fasta

output_dir1=$output_dir/Step1_BamCellTypes/filtered
output_dir2=$output_dir/Step2_BaseCellCounts

# Create directories
mkdir -p "$output_dir2/AMPK_negative"
mkdir -p "$output_dir2/AMPK_positive"

# Cell types
cell_types=("AMPK_negative" "AMPK_positive")

# Loop through cell types
for cell_type in "${cell_types[@]}"; do
    bam="$output_dir1/${project}.${cell_type}_filtered.bam"  # Adjusted BAM path
    echo "Processing $bam"

    # Temp folder
    temp="$output_dir2/${cell_type}/temp"
    mkdir -p "$temp"

    # Run BaseCellCounter
    python "$SCOMATIC/scripts/BaseCellCounter/BaseCellCounter.py" \
        --bam "$bam" \
        --ref "$REF" \
        --chrom all \
        --out_folder "$output_dir2/${cell_type}" \
        --min_bq 30 \
        --tmp_dir "$temp" \
        --nprocs 16  # Important: Adjust based on your resources

    rm -rf "$temp"
done

echo "Step 2 complete"

################################################################
############# Step 3: Merging base count matrices ###############
################################################################

output_dir3=$output_dir/Step3_BaseCellCountsMerged

# Create output directory
mkdir -p "$output_dir3"

# Loop through cell types
for cell_type in "${cell_types[@]}"; do
    echo "Merging counts for $cell_type"
    python "$SCOMATIC/scripts/MergeCounts/MergeBaseCellCounts.py" \
        --tsv_folder "$output_dir2/${cell_type}" \
        --outfile "$output_dir3/${project}.BaseCellCounts.${cell_type}.tsv"
done

echo "Step 3 complete"


################################################################
############# Step 4: Detection of somatic mutations ###############
################################################################

output_dir4=$output_dir/Step4_VariantCalling

# Create output directories
mkdir -p "$output_dir4/AMPK_negative"
mkdir -p "$output_dir4/AMPK_positive"

editing=$SCOMATIC/RNAediting/AllEditingSites.hg38.txt
PON=$SCOMATIC/PoNs/PoN.scRNAseq.hg38.tsv

# Loop through cell types
for cell_type in "${cell_types[@]}"; do
    echo "Calling variants for $cell_type"

    infile="$output_dir3/${project}.BaseCellCounts.${cell_type}.tsv"
    outfile_prefix="$output_dir4/${cell_type}/${project}"

    # Step 4.1
    python "$SCOMATIC/scripts/BaseCellCalling/BaseCellCalling.step1.py" \
        --infile "$infile" \
        --outfile "$outfile_prefix" \
        --ref "$REF"

    # Step 4.2
    python "$SCOMATIC/scripts/BaseCellCalling/BaseCellCalling.step2.py" \
        --infile "${outfile_prefix}.calling.step1.tsv" \
        --outfile "$outfile_prefix" \
        --editing "$editing" \
        --pon "$PON"

    # (Optional) Intersect with high-quality regions
    bedtools intersect -header -a "${outfile_prefix}.calling.step2.tsv" -b "$SCOMATIC/bed_files_of_interest/UCSC.k100_umap.without.repeatmasker.bed" | awk '$1 ~ /^#/ || $6 == "PASS"' > "${outfile_prefix}.calling.step2.pass.tsv"
done

echo "Step 4 complete"






