#!/bin/bash

#PBS -l select=3:ncpus=64:mem=256g
#PBS -l walltime=24:00:00
#PBS -P 11003581
#PBS -N run-somatic-scDNA
#PBS -j oe

source /home/project/11003581/Tools/miniforge3/bin/activate
conda activate /home/project/11003581/conda-envs/sccaller

cd $PBS_O_WORKDIR

# Set paths
BAM_BASE_DIR="/home/users/nus/ash.ps/scratch/scDNA/DNA/Secondary-Analysis-DNA"
BULK_BAM="/home/users/nus/ash.ps/scratch/scDNA/Bulk-data/Secondary-Analysis/2204/2204.sorted.bqsr.dedup.bam"
FASTA="/home/project/11003581/Ref/Homo_sapiens/Ensembl/Homo_sapiens.GRCh38.dna_rm.primary_assembly.fa"
SNP_IN="/home/users/nus/ash.ps/scratch/scDNA/analysis/somatic/2204.g.vcf"
OUT_BASE="/home/users/nus/ash.ps/scratch/scDNA/analysis/somatic"
SCCALLER_SCRIPT="/home/users/nus/ash.ps/scratch/scDNA/analysis/logs/sccaller_v2.0.0.py"

# Find all subfolders starting with 2204
for CELL_DIR in "$BAM_BASE_DIR"/2204*; do
    # Get cell name (basename)
    CELL_NAME=$(basename "$CELL_DIR")

    # Define the BAM file for this cell
    CELL_BAM="$CELL_DIR/$CELL_NAME.sorted.bqsr.dedup.bam"

    # Check if BAM file exists
    if [ ! -f "$CELL_BAM" ]; then
        echo "Warning: BAM file not found for $CELL_NAME, skipping."
        continue
    fi

    # Create output and log directories
    OUT_DIR="$OUT_BASE/$CELL_NAME"
    mkdir -p "$OUT_DIR"

    # Run sccaller and redirect logs
    python "$SCCALLER_SCRIPT" \
        --bam "$CELL_BAM" \
        --bulk "$BULK_BAM" \
        --fasta "$FASTA" \
        --output "$OUT_DIR/somatic-$CELL_NAME.vcf" \
        --snp_type hsnp \
        --snp_in "$SNP_IN" \
        --cpu_num 128 \
        --engine samtools \
        > "$OUT_DIR/sccaller_${CELL_NAME}.log" 2>&1

done
