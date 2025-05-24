#!/bin/bash

#PBS -N bulk-somatic-VC
#PBS -l select=3:ncpus=64:mem=256g          
#PBS -l walltime=24:00:00
#PBS -P 11003581
#PBS -j oe
#PBS -M ash.ps@nus.edu.sg

# ===== ENVIRONMENT =====
module purge
module load nextflow/24.10.5
module load miniforge3                 

export NXF_CONDA_CACHEDIR=/home/project/11003581/Tools/nf-conda-cache
export NXF_CONDA_USE_MAMBA=true

# ===== RUN SAREK =====
nextflow run nf-core/sarek -r 3.5.1 \
    -profile conda \
    --input /home/users/nus/ash.ps/scratch/scDNA/analysis/bulk-somatic/samplesheet.csv \
    --step variant_calling \
    --tools mutect2,manta,ascat \
    --genome GATK.GRCh38 \
    --outdir /home/users/nus/enambis/scratch/NCCS-MASH/bulk/ \
    --igenomes_base /home/project/11003581/Ref/sarek-refs/



### samplesheet.csv
patient,sex,status,sample,bam,bai
2204,XX,1,2204_bulk,/home/users/nus/ash.ps/scratch/scDNA/Bulk-data/Secondary-Analysis/2204/2204.sorted.bqsr.dedup.bam,/home/users/nus/ash.ps/scratch/scDNA/Bulk-data/Secondary-Analysis/2204/2204.sorted.bqsr.dedup.bam.bai
3401,XX,1,3401_bulk,/home/users/nus/ash.ps/scratch/scDNA/Bulk-data/Secondary-Analysis/3401/3401.sorted.bqsr.dedup.bam,/home/users/nus/ash.ps/scratch/scDNA/Bulk-data/Secondary-Analysis/3401/3401.sorted.bqsr.dedup.bam.bai

