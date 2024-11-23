#!/bin/bash

#PBS -l select=1:ncpus=32:mem=128g
#PBS -l walltime=1:00:00
#PBS -P 11003581
#PBS -N p3l-coverage
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load samtools

#samtools view -b /home/project/11003581/Data/Ash/P3L-lab-analysis/sorted_bams/p3292l_hac_sorted.bam chr13 chr12 > /home/project/11003581/Data/Ash/P3L-lab-analysis/coverage/chr12_13_reads.bam

#samtools depth /home/project/11003581/Data/Ash/P3L-lab-analysis/sorted_bams/p3292l_hac_sorted.bam > /home/project/11003581/Data/Ash/P3L-lab-analysis/coverage/p3l_coverage.txt

samtools coverage /home/project/11003581/Data/Ash/P3L-lab-analysis/coverage/chr12_13_reads.bam -o /home/project/11003581/Data/Ash/P3L-lab-analysis/coverage/cov-out