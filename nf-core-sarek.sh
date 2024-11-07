#!/bin/bash

#PBS -l select=1:ncpus=64:mem=256g
#PBS -l walltime=20:00:00
#PBS -P 11003581
#PBS -N trial-nf-sarek
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR


module load java/17.0.6-jdk
module load singularity

/home/project/11003581/Tools/nextflow  run nf-core/sarek \
   -profile singularity \
   --input /home/project/11003581/Data/Ash/nf-core/samplesheet.csv \
   --outdir /home/project/11003581/Data/Ash/nf-core/




'''

patient,sample,lane,fastq_1,fastq_2
N01,DKDN240040899-1A_HVNMJDSXC,L1,/home/project/11003581/Data/YS/01.RawData/N01/N01_DKDN240040899-1A_HVNMJDSXC_L1_1.fq.gz,/home/project/11003581/Data/YS/01.RawData/N01/N01_DKDN240040899-1A_HVNMJDSXC_L1_2.fq.gz


'''

#resuming run

