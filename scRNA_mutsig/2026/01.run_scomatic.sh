#!/bin/bash

#PBS -l select=2:ncpus=64:mem=128g
#PBS -l walltime=10:00:00
#PBS -P 11003581
#PBS -j oe
#PBS -N scomatic-step2-3

set -euo pipefail
cd "$PBS_O_WORKDIR"

LOGFILE="$PBS_O_WORKDIR/scomatic-step2-3.log"
exec > >(tee -a "$LOGFILE") 2>&1

echo "[INFO] Resuming Step 2 + 3 at $(date)"

# ------------------------------------------------------------------------------
# ENVIRONMENT
# ------------------------------------------------------------------------------
source /home/project/11003581/Tools/miniforge3/bin/activate
conda activate /home/project/11003581/conda-envs/scMultiomics

module load samtools/1.15.1
module load python/3.12.1-gcc11
module load bedtools/2.30.0

echo "[INFO] Starting GLOBAL SComatic pipeline"

############################################################
# INPUTS
############################################################

base_dir="/home/users/nus/ash.ps/scratch/scRNA-mutsig"
bam="$base_dir/data/Breast_Cancer_3p_possorted_genome_bam.bam"
meta="$base_dir/cell_annotation_AMPK.txt"

SCOMATIC="/home/project/11003581/Tools/SComatic"

REF="/home/project/11003581/Ref/Homo_sapiens/GATK/hg38/Homo_sapiens_assembly38.fasta"
editing="$SCOMATIC/RNAediting/AllEditingSites.hg38.txt"

outdir="$base_dir/scomatic_global"
mkdir -p $outdir

############################################################
# STEP 1 — Split BAM by decile (ONLY FOR LABELING)
############################################################

echo "[STEP 1] Splitting BAM by AMPK deciles"

python $SCOMATIC/scripts/SplitBam/SplitBamCellTypes.py \
  --bam "$bam" \
  --meta "$meta" \
  --id "AMPK" \
  --outdir "$outdir/Step1_split"

echo "[STEP 1 DONE]"

############################################################
# STEP 2 — Merge BAMs back (IMPORTANT TRICK)
############################################################

echo "[STEP 2] Merging BAMs back for GLOBAL calling"

samtools merge -@ 32 $outdir/merged.bam $outdir/Step1_split/*.bam
samtools index $outdir/merged.bam

echo "[STEP 2 DONE]"

############################################################
# STEP 3 — BaseCellCounter (GLOBAL)
############################################################

echo "[STEP 3] Running BaseCellCounter"

python $SCOMATIC/scripts/BaseCellCounter/BaseCellCounter.py \
  --bam "$outdir/merged.bam" \
  --ref "$REF" \
  --chrom all \
  --out_folder "$outdir/Step2_counts" \
  --min_bq 20 \
  --min_ac 1 \
  --min_af 0.01 \
  --min_dp 2 \
  --min_cc 1 \
  --nprocs 64

echo "[STEP 3 DONE]"

############################################################
# STEP 4 — Variant calling (GLOBAL)
############################################################

echo "[STEP 4.1] Step1"

python $SCOMATIC/scripts/BaseCellCalling/BaseCellCalling.step1.py \
  --infile "$outdir/Step2_counts/"*.tsv \
  --outfile "$outdir/Step4_variants/merged" \
  --ref "$REF"

echo "[STEP 4.2] Step2 (LIGHT FILTER)"

python $SCOMATIC/scripts/BaseCellCalling/BaseCellCalling.step2.py \
  --infile "$outdir/Step4_variants/merged.calling.step1.tsv" \
  --outfile "$outdir/Step4_variants/merged" \
  --editing "$editing"

echo "[DONE] GLOBAL SComatic complete"