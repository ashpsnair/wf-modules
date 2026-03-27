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

# ------------------------------------------------------------------------------
# PATHS
# ------------------------------------------------------------------------------
SCOMATIC="/home/project/11003581/Tools/SComatic"

filtered_out="/home/users/nus/ash.ps/scratch/scRNA-mutsig/analysis/Step1_BamCellTypes/filtered"
output_dir2="/home/users/nus/ash.ps/scratch/scRNA-mutsig/analysis/Step2_BaseCellCounts"
output_dir4="/home/users/nus/ash.ps/scratch/scRNA-mutsig/analysis/Step4_VariantCalling"

REF="/home/project/11003581/Ref/Homo_sapiens/GATK/hg38/Homo_sapiens_assembly38.fasta"
editing="$SCOMATIC/RNAediting/AllEditingSites.hg38.txt"
PON="$SCOMATIC/PoNs/PoN.scRNAseq.hg38.tsv"

mkdir -p "$output_dir2" "$output_dir4"

# ------------------------------------------------------------------------------
# LOOP THROUGH EXISTING BAMs (FROM STEP 1)
# ------------------------------------------------------------------------------
for bam in "$filtered_out"/*.bam; do
    
    cell_name=$(basename "$bam" .bam)   # KEEP ORIGINAL NAME (important)

    echo "--------------------------------------"
    echo "[INFO] Processing: $cell_name"

    basecount_out="${output_dir2}/${cell_name}.tsv"
    step1_out="${output_dir4}/${cell_name}.calling.step1.tsv"
    step2_out="${output_dir4}/${cell_name}.calling.step2.tsv"

    temp="${output_dir2}/temp_${cell_name}"
    mkdir -p "$temp"

    ############################################################
    # STEP 2: BaseCellCounter (skip if done)
    ############################################################
    if [ ! -f "$basecount_out" ]; then
        echo "[STEP 2] Running BaseCellCounter"

        python "$SCOMATIC/scripts/BaseCellCounter/BaseCellCounter.py" \
            --bam "$bam" \
            --ref "$REF" \
            --chrom all \
            --out_folder "$output_dir2" \
            --min_bq 20 \
            --min_ac 1 \
            --min_af 0.01 \
            --min_dp 2 \
            --min_cc 1 \
            --tmp_dir "$temp" \
            --nprocs 128
    else
        echo "[SKIP] BaseCellCounter already exists"
    fi

    ############################################################
    # STEP 3.1
    ############################################################
    if [ ! -f "$step1_out" ]; then
        echo "[STEP 3.1] Variant calling step 1"

        python "$SCOMATIC/scripts/BaseCellCalling/BaseCellCalling.step1.py" \
            --infile "$basecount_out" \
            --outfile "${output_dir4}/${cell_name}" \
            --ref "$REF"
    else
        echo "[SKIP] Step1 already exists"
    fi

    ############################################################
    # STEP 3.2
    ############################################################
    if [ ! -f "$step2_out" ]; then
        echo "[STEP 3.2] Variant calling step 2"

        python "$SCOMATIC/scripts/BaseCellCalling/BaseCellCalling.step2.py" \
            --infile "$step1_out" \
            --outfile "${output_dir4}/${cell_name}" \
            --editing "$editing" \
            --pon "$PON"
    else
        echo "[SKIP] Step2 already exists"
    fi

    ############################################################
    # CLEAN TEMP
    ############################################################
    rm -rf "$temp"

    echo "[DONE] $cell_name"

done

echo "[INFO] Step 2 + 3 completed at $(date)"











#########TEST

############################################################
# DEBUG: Test without aggressive filtering
############################################################
source /home/project/11003581/Tools/miniforge3/bin/activate
conda activate /home/project/11003581/conda-envs/scMultiomics

module load samtools/1.15.1
module load python/3.12.1-gcc11
module load bedtools/2.30.0

# ------------------------------------------------------------------------------
# PATHS
# ------------------------------------------------------------------------------
SCOMATIC="/home/project/11003581/Tools/SComatic"

filtered_out="/home/users/nus/ash.ps/scratch/scRNA-mutsig/analysis/Step1_BamCellTypes/filtered"
output_dir2="/home/users/nus/ash.ps/scratch/scRNA-mutsig/analysis/Step2_BaseCellCounts"
output_dir4="/home/users/nus/ash.ps/scratch/scRNA-mutsig/analysis/Step4_VariantCalling"

REF="/home/project/11003581/Ref/Homo_sapiens/GATK/hg38/Homo_sapiens_assembly38.fasta"
editing="$SCOMATIC/RNAediting/AllEditingSites.hg38.txt"
PON="$SCOMATIC/PoNs/PoN.scRNAseq.hg38.tsv"

cd $output_dir4
test_sample="ampk.D1_filtered"
infile="${test_sample}.calling.step1.tsv"

python "$SCOMATIC/scripts/BaseCellCalling/BaseCellCalling.step2.py" \
  --infile "$infile" \
  --outfile "${test_sample}_NOFILTER"