#!/bin/bash

#PBS -l select=1
#PBS -l walltime=1:00:00
#PBS -P 11003581
#PBS -N trial-mosdepth
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load miniforge3
module load gcc
module load python/3.12.1-gcc11
conda activate myenv

#output_dir=/home/project/11003581/Data/Ash/P3L-lab-analysis/mosdepth/p3l
#bamfile=/home/project/11003581/Data/Ash/P3L-lab-analysis/sorted_bams/p3292l_hac_sorted.bam

#mosdepth -t 8 -x -b 1000 -Q 20 $output_dir $bamfile

spectre CNVCaller \
  --coverage /home/project/11003581/Data/Ash/P3L-lab-analysis/mosdepth/p3l.regions.bed.gz \
  --sample-id p3l \
  --output-dir /home/project/11003581/Data/Ash/P3L-lab-analysis/spectre_out/ \
  --reference r/home/project/11003581/Ref/hg38.fa.gz


