#!/bin/bash

#PBS -l select=3:ncpus=64:mem=256g
#PBS -l walltime=24:00:00
#PBS -P 11003581
#PBS -N sarek-PD-b1
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR


module load java/17.0.6-jdk
module load singularity

/home/project/11003581/Tools/nextflow run nf-core/sarek -r 3.4.4 \
   -profile singularity \
   --input /home/users/nus/ash.ps/scratch/T2DM/analysis-PD/samplesheet.csv \
   --outdir /home/users/nus/ash.ps/scratch/T2DM/analysis-PD \
   --tools mutect2,manta \
   --pon /home/project/11003581/Ref/pons/1000g_pon.hg38.vcf.gz \
   --pon_tbi /home/project/11003581/Ref/pons/1000g_pon.hg38.vcf.gz.tbi

#  PD trial
patient,status,sample,lane,fastq_1,fastq_2
SRR29413853,1,SRR29413853,lane_1,/home/users/nus/ash.ps/scratch/T2DM/PD/fastq/SRR29413853_pass_1.fastq.gz,/home/users/nus/ash.ps/scratch/T2DM/PD/fastq/SRR29413853_pass_2.fastq.gz
SRR29413854,1,SRR29413854,lane_1,/home/users/nus/ash.ps/scratch/T2DM/PD/fastq/SRR29413854_pass_1.fastq.gz,/home/users/nus/ash.ps/scratch/T2DM/PD/fastq/SRR29413854_pass_2.fastq.gz

