#!/bin/bash

#PBS -l select=1:ncpus=64:mem=256g
#PBS -l walltime=6:00:00
#PBS -P 11003581
#PBS -N trial-gatk
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load samtools

/home/project/11003581/Tools/bwa-mem2-2.2.1_x64-linux/bwa-mem2 mem -t 128 \
    /home/project/11003581/Ref/Homo_sapiens/GATK/hg38/Homo_sapiens_assembly38.fasta \
    /home/project/11003581/Data/YS-analysis/raw_data/N01_1.fq.gz \
    /home/project/11003581/Data/YS-analysis/raw_data/N01_2.fq.gz >/home/users/nus/ash.ps/scratch/YS-analysis/bwa2-bams/p08-normal.sam

samtools view -bS /home/users/nus/ash.ps/scratch/YS-analysis/bwa2-bams/p08-normal.sam > /home/users/nus/ash.ps/scratch/YS-analysis/bwa2-bams/p08-normal.bam