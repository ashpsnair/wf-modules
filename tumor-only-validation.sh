patient,sex,status,sample,lane,fastq_1,fastq_2
patient01,XY,1,T01,lane_1,/home/project/11003581/Data/YS-analysis/raw_data/T01_1.fq.gz,/home/project/11003581/Data/YS-analysis/raw_data/T01_2.fq.gz
patient02,XY,1,T02,lane_1,/home/project/11003581/Data/YS-analysis/raw_data/T02_1.fq.gz,/home/project/11003581/Data/YS-analysis/raw_data/T02_2.fq.gz
patient03,XX,1,T03,lane_1,/home/project/11003581/Data/YS-analysis/raw_data/T03_1.fq.gz,/home/project/11003581/Data/YS-analysis/raw_data/T03_2.fq.gz
patient04,XY,1,T04,lane_1,/home/project/11003581/Data/YS-analysis/raw_data/T04_1.fq.gz,/home/project/11003581/Data/YS-analysis/raw_data/T04_2.fq.gz
patient05,XY,1,T05,lane_1,/home/project/11003581/Data/YS-analysis/raw_data/T05_1.fq.gz,/home/project/11003581/Data/YS-analysis/raw_data/T05_2.fq.gz
patient06,XX,1,T06,lane_1,/home/project/11003581/Data/YS-analysis/raw_data/T06_1.fq.gz,/home/project/11003581/Data/YS-analysis/raw_data/T06_2.fq.gz
patient07,XY,1,T07,lane_1,/home/project/11003581/Data/YS-analysis/raw_data/T07_1.fq.gz,/home/project/11003581/Data/YS-analysis/raw_data/T07_2.fq.gz
patient08,XY,1,T08,lane_1,/home/project/11003581/Data/YS-analysis/raw_data/T08_1.fq.gz,/home/project/11003581/Data/YS-analysis/raw_data/T08_2.fq.gz
patient09,XY,1,T09,lane_1,/home/project/11003581/Data/YS-analysis/raw_data/T09_1.fq.gz,/home/project/11003581/Data/YS-analysis/raw_data/T09_2.fq.gz
patient10,XY,1,T10,lane_1,/home/project/11003581/Data/YS-analysis/raw_data/T10_1.fq.gz,/home/project/11003581/Data/YS-analysis/raw_data/T10_2.fq.gz
patient11,XY,1,T11,lane_1,/home/project/11003581/Data/YS-analysis/raw_data/T11_1.fq.gz,/home/project/11003581/Data/YS-analysis/raw_data/T11_2.fq.gz
patient12,XY,1,T12,lane_1,/home/project/11003581/Data/YS-analysis/raw_data/T12_1.fq.gz,/home/project/11003581/Data/YS-analysis/raw_data/T12_2.fq.gz
patient13,XY,1,T13,lane_1,/home/project/11003581/Data/YS-analysis/raw_data/T13_1.fq.gz,/home/project/11003581/Data/YS-analysis/raw_data/T13_2.fq.gz
patient14,XY,1,T14,lane_1,/home/project/11003581/Data/YS-analysis/raw_data/T14_1.fq.gz,/home/project/11003581/Data/YS-analysis/raw_data/T14_2.fq.gz


#!/bin/bash

#PBS -l select=3:ncpus=64:mem=128g
#PBS -l walltime=24:00:00
#PBS -P 11003581
#PBS -N tumor-only-validation
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR


module load java/17.0.6-jdk
module load singularity

/home/project/11003581/Tools/nextflow run nf-core/sarek -r 3.4.4 \
   -profile singularity \
   --input /home/users/nus/ash.ps/scratch/tumor-only-validation/samplesheet.csv \
   --outdir /home/users/nus/ash.ps/scratch/tumor-only-validation/ \
   --tools mutect2,snpeff,vep,msisensorpro \
   --pon /home/project/11003581/Ref/pons/1000g_pon.hg38.vcf.gz \
   --pon_tbi /home/project/11003581/Ref/pons/1000g_pon.hg38.vcf.gz.tbi