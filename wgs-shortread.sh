#!/bin/bash

#PBS -l select=1:ncpus=64:mem=128g:ngpus=1
#PBS -l walltime=10:00:00
#PBS -P 11003581
#PBS -N wgs-auto
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load gcc
module load python/3.12.1-gcc11

/home/project/11003581/Tools/bwa-mem2-2.2.1_x64-linux/bwa-mem2 mem /home/project/11003581/Ref/Homo_sapiens/GATK/hg38/Homo_sapiens_assembly38.fasta \
/home/project/11003581/Data/YS-analysis/raw_data/N08_1.fq.gz \
/home/project/11003581/Data/YS-analysis/raw_data/N08_2.fq.gz > /home/users/nus/ash.ps/scratch/wgs-auto/N08.sam

'''
ERROR! Unable to open the file: /home/project/11003581/Ref/Homo_sapiens/GATK/hg38/Homo_sapiens_assembly38.fasta.bwt.2bit.64
'''
###gatk- mutect

gatk Mutect2 \
     -R /home/project/11003581/Ref/Homo_sapiens/GATK/hg38/Homo_sapiens_assembly38.fasta \
     -I /home/users/nus/ash.ps/scratch/YS-batch2/preprocessing/recalibrated/T08/T08.recal.cram \
     -I /home/users/nus/ash.ps/scratch/YS-batch2/preprocessing/recalibrated/N08/N08.recal.cram \
     -normal T08_vs_N08 \
     -O /home/users/nus/ash.ps/scratch/wgs-auto/P08_somatic.vcf.gz
