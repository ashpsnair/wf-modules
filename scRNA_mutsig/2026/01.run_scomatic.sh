#!/bin/bash

#PBS -l select=2:ncpus=64:mem=128g
#PBS -l walltime=12:00:00
#PBS -P 11003581
#PBS -j oe
#PBS -N scomatic_pipeline

set -euo pipefail
cd "$PBS_O_WORKDIR"

echo "Pipeline started at $(date)"

###############################################
# ENVIRONMENT
###############################################

source /home/project/11003581/Tools/miniforge3/bin/activate
conda activate /home/project/11003581/conda-envs/scMultiomics

module load samtools/1.15.1
module load bedtools/2.30.0
module load python/3.12.1-gcc11
module load r/4.2.0

###############################################
# PATHS
###############################################

SCOMATIC=/home/project/11003581/Tools/SComatic

DATA=/home/users/nus/ash.ps/scratch/mulitomics/example/data
ANALYSIS=/home/users/nus/ash.ps/scratch/mulitomics/example/analysis

META=${DATA}/Example.cell_barcode_annotations.tsv

REF=/home/project/11003581/Ref/Homo_sapiens/GATK/hg38/Homo_sapiens_assembly38.fasta
EDITING=${SCOMATIC}/RNAediting/AllEditingSites.hg38.txt
PON=${SCOMATIC}/PoNs/PoN.scRNAseq.hg38.tsv

THREADS=128

###############################################
# OUTPUT STRUCTURE
###############################################

STEP1=${ANALYSIS}/Step1_CellTypeBams
STEP2=${ANALYSIS}/Step2_BaseCellCounts
STEP3=${ANALYSIS}/Step3_MergedCounts
STEP4=${ANALYSIS}/Step4_VariantCalling
TRI=${ANALYSIS}/TrinucleotideContext

mkdir -p $STEP1 $STEP2 $STEP3 $STEP4 $TRI

###############################################
# INPUT BAM
###############################################

BAM=$(find $DATA -name "*.bam")

###############################################
# STEP 1 — SPLIT BAM BY CELL TYPE
###############################################

echo "Splitting BAM by cell type..."

python ${SCOMATIC}/scripts/SplitBam/SplitBamCellTypes.py \
    --bam $BAM \
    --meta $META \
    --id bc \
    --n_trim 5 \
    --max_nM 5 \
    --max_NH 1 \
    --outdir $STEP1


###############################################
# STEP 1.5 — SORT + INDEX SPLIT BAMS
###############################################

echo "Sorting and indexing split BAMs..."

for bam in ${STEP1}/*.bam
do

name=$(basename $bam .bam)

samtools sort -@ $THREADS -o ${STEP1}/${name}_sorted.bam $bam

samtools index ${STEP1}/${name}_sorted.bam

done


###############################################
# STEP 2 — BASE CELL COUNTS
###############################################

echo "Running BaseCellCounter..."

for bam in ${STEP1}/*_sorted.bam
do

cell=$(basename $bam _sorted.bam)

TMP=${STEP2}/tmp_${cell}

mkdir -p $TMP

python ${SCOMATIC}/scripts/BaseCellCounter/BaseCellCounter.py \
    --bam $bam \
    --ref $REF \
    --chrom all \
    --out_folder $STEP2 \
    --min_bq 20 \
    --min_dp 2 \
    --min_cc 1 \
    --min_ac 1 \
    --min_af 0.01 \
    --tmp_dir $TMP \
    --nprocs $THREADS

rm -rf $TMP

done


###############################################
# STEP 3 — MERGE COUNTS
###############################################

echo "Merging base cell counts..."

python ${SCOMATIC}/scripts/MergeCounts/MergeBaseCellCounts.py \
    --tsv_folder $STEP2 \
    --outfile ${STEP3}/bc.BaseCellCounts.AllCellTypes.tsv


###############################################
# STEP 4 — VARIANT CALLING
###############################################

echo "Running SComatic variant calling..."

python ${SCOMATIC}/scripts/BaseCellCalling/BaseCellCalling.step1.py \
    --infile ${STEP3}/bc.BaseCellCounts.AllCellTypes.tsv \
    --outfile ${STEP4}/bc \
    --ref $REF


python ${SCOMATIC}/scripts/BaseCellCalling/BaseCellCalling.step2.py \
    --infile ${STEP4}/bc.calling.step1.tsv \
    --outfile ${STEP4}/bc \
    --editing $EDITING \
    --pon $PON


###############################################
# TRINUCLEOTIDE BACKGROUND
###############################################

echo "Calculating trinucleotide context..."

ls ${STEP4}/*.calling.step1.tsv > ${TRI}/step1_files.txt

python ${SCOMATIC}/scripts/TrinucleotideBackground/TrinucleotideContextBackground.py \
    --in_tsv ${TRI}/step1_files.txt \
    --out_file ${TRI}/TrinucleotideBackground.txt

echo "Pipeline finished at $(date)"