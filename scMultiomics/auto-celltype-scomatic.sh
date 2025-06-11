#!/bin/bash

#PBS -l select=2:ncpus=64:mem=128g
#PBS -l walltime=10:00:00
#PBS -P 11003581
#PBS -j oe
#PBS -N run-celltype-scomatic

set -euo pipefail
cd "$PBS_O_WORKDIR"

LOGFILE="$PBS_O_WORKDIR/scomatic-run.log"
exec > >(tee -a "$LOGFILE") 2>&1

echo "[INFO] SComatic pipeline started at $(date)"

# -- Environment Setup --------------------------------------------------------
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

# -- Paths & Params -----------------------------------------------------------
project="bc"
SCOMATIC="/home/project/11003581/Tools/SComatic"

base_dir="/home/users/nus/ash.ps/scratch/mulitomics/10x_data/data"
output_dir="/home/users/nus/ash.ps/scratch/mulitomics/10x_data/analysis/"
meta_celltype="/home/users/nus/ash.ps/scratch/mulitomics/10x_data/data/cell_barcode_annotations.tsv"
REF="/home/project/11003581/Ref/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa"
editing="$SCOMATIC/RNAediting/AllEditingSites.hg38.txt"
PON="$SCOMATIC/PoNs/PoN.scRNAseq.hg38.tsv"
threads=128

bam_file=$(find "$base_dir" -name "*.bam" -type f)

# ------------------------------------------------------------------------------
# STEP 1: SPLIT BAM BY CELL TYPE
# ------------------------------------------------------------------------------
echo "[STEP 1.1] Splitting BAM by cell type at $(date)"

output_dir1="$output_dir/Step1_BamCellTypes"
mkdir -p "$output_dir1" 

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
# STEP 2: BASE COUNTING (per cell type)
# ------------------------------------------------------------------------------
echo "[STEP 2] Base counting for all filtered BAMs at $(date)"

output_dir2="$output_dir/Step2_BaseCellCounts"
mkdir -p "$output_dir2"

for bam in $output_dir1/*.bam; do
    cell_type=$(basename "$bam" .bam)
    temp="$output_dir2/temp_${cell_type}"
    mkdir -p "$temp"

    echo "[STEP 2] BaseCellCounter for $cell_type"
    python "$SCOMATIC/scripts/BaseCellCounter/BaseCellCounter.py" \
        --bam "$bam" \
        --ref "$REF" \
        --chrom all \
        --out_folder "$output_dir2" \
        --min_bq 30 \
        --tmp_dir "$temp" \
        --nprocs 128

    echo "[INFO] Cleanup: $temp"
    rm -rf "$temp"
done
echo "[STEP 2] All BaseCellCounting done at $(date)"

# ------------------------------------------------------------------------------
# STEP 3: VARIANT CALLING (per cell type)
# ------------------------------------------------------------------------------
echo "[STEP 3] Variant calling per cell type at $(date)"

output_dir3="$output_dir/Step3_VariantCalling"
mkdir -p "$output_dir3"

for tsv in $output_dir2/*.tsv; do
    cell_type=$(basename "$tsv" .tsv)
    # Step 3.1: Variant calling step1 for this cell type
    python $SCOMATIC/scripts/BaseCellCalling/BaseCellCalling.step1.py \
        --infile "$tsv" \
        --outfile "$output_dir3/${cell_type}" \
        --ref "$REF"
    # Step 3.2: Variant calling step2 for this cell type
    python $SCOMATIC/scripts/BaseCellCalling/BaseCellCalling.step2.py \
        --infile "$output_dir3/${cell_type}.calling.step1.tsv" \
        --outfile "$output_dir3/${cell_type}" \
        --editing $editing \
        --pon $PON
done
echo "[STEP 3] Variant calling per cell type done at $(date)"

# ------------------------------------------------------------------------------
# STEP 4: CELL TYPE CALLABLE SITES (optional, recommended by SComatic)
# ------------------------------------------------------------------------------
output_dir4=$output_dir/CellTypeCallableSites
mkdir -p $output_dir4

for call1 in $output_dir3/*.calling.step1.tsv; do
    cell_type=$(basename "$call1" .calling.step1.tsv)
    python $SCOMATIC/scripts/GetCallableSites/GetAllCallableSites.py \
        --infile "$call1" \
        --outfile "$output_dir4/${cell_type}" \
        --max_cov 150 --min_cell_types 2
done

# ------------------------------------------------------------------------------
# STEP 5: Computing the number of callable sites per CELL (optional)
# ------------------------------------------------------------------------------
output_dir5=$output_dir/UniqueCellCallableSites
mkdir -p $output_dir5

for bam in $output_dir1/*.bam; do  
    cell_type=$(basename "$bam" .bam)
    temp=$output_dir5/temp_${cell_type}
    mkdir -p $temp
    # You may want to use each cell's variant file here instead of merged
    python $SCOMATIC/scripts/SitesPerCell/SitesPerCell.py \
        --bam $bam \
        --infile "$output_dir3/${cell_type}.calling.step1.tsv" \
        --ref $REF \
        --out_folder $output_dir5 --tmp_dir $temp --nprocs 1
    rm -rf $temp
done

# ------------------------------------------------------------------------------
# STEP 6: Computing the genotype for each cell at the variant sites (optional)
# ------------------------------------------------------------------------------
output_dir6=$output_dir/SingleCellAlleles
mkdir -p $output_dir6

for bam in $output_dir1/*.bam; do  
    cell_type=$(basename "$bam" .bam)
    temp=$output_dir6/temp_${cell_type}
    mkdir -p $temp

    python $SCOMATIC/scripts/SingleCellGenotype/SingleCellGenotype.py --bam $bam  \
        --infile "$output_dir3/${cell_type}.calling.step2.tsv" \
        --nprocs 1   \
        --meta $meta_celltype   \
        --outfile ${output_dir6}/${cell_type}.single_cell_genotype.tsv  \
        --tmp_dir $temp  \
        --ref $REF

    rm -rf $temp
done

# ------------------------------------------------------------------------------
# STEP 7: TRINUCLEOTIDE CONTEXT (optional, if needed)
# ------------------------------------------------------------------------------
output_dir7=$output_dir/TrinucleotideContext
mkdir -p $output_dir7

# List all step1 variant tsvs (per cell type) in a text file
ls $output_dir3/*.calling.step1.tsv > $output_dir7/step1_files.txt

python $SCOMATIC/scripts/TrinucleotideBackground/TrinucleotideContextBackground.py \
        --in_tsv $output_dir7/step1_files.txt \
        --out_file $output_dir7/TrinucleotideBackground.txt

echo "[DONE] SComatic full pipeline completed at $(date)"
