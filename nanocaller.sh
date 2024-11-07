#!/bin/bash

#PBS -l select=1:ngpus=1
#PBS -l walltime=10:00:00
#PBS -P personal-ash.ps
#PBS -N trial-clairs
#PBS -j oe

module load miniforge3
conda activate myenv

#conda install -c bioconda nanocaller

#sample details
ref_fasta="/home/project/11003581/Ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
p3l_bam="/home/users/nus/ash.ps/scratch/P3L-lab-analysis/sorted_bams/p3292l_hac_sorted.bam"
wt_bam="/home/users/nus/ash.ps/scratch/P3L-lab-analysis/sorted_bams/WT_MCF10A_hac_sorted.bam"

NanoCaller --bam $p3l_bam --ref $ref_fasta --cpu 10 --mode snps
NanoCaller --bam $wt_bam --ref $ref_fasta --cpu 10 --mode snps