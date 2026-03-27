

#samplesheet.csv

patient,sex,status,sample,bam,bai
TCGA-L7-A6VZ,XX,1,TCGA-L7-A6VZ,/home/users/nus/ash.ps/scratch/RNA-mut/RNA-bams/cleaned_bams/TCGA-L7-A6VZ_cleaned.bam,/home/users/nus/ash.ps/scratch/RNA-mut/RNA-bams/cleaned_bams/TCGA-L7-A6VZ_cleaned.bam.bai
TCGA-L5-A4OP,XX,1,TCGA-L5-A4OP,/home/users/nus/ash.ps/scratch/RNA-mut/RNA-bams/cleaned_bams/TCGA-L5-A4OP_cleaned.bam,/home/users/nus/ash.ps/scratch/RNA-mut/RNA-bams/cleaned_bams/TCGA-L5-A4OP_cleaned.bam.bai
TCGA-IG-A7DP,XX,1,TCGA-IG-A7DP,/home/users/nus/ash.ps/scratch/RNA-mut/RNA-bams/cleaned_bams/TCGA-IG-A7DP_cleaned.bam,/home/users/nus/ash.ps/scratch/RNA-mut/RNA-bams/cleaned_bams/TCGA-IG-A7DP_cleaned.bam.bai
TCGA-L5-A43E,XX,1,TCGA-L5-A43E,/home/users/nus/ash.ps/scratch/RNA-mut/RNA-bams/cleaned_bams/TCGA-L5-A43E_cleaned.bam,/home/users/nus/ash.ps/scratch/RNA-mut/RNA-bams/cleaned_bams/TCGA-L5-A43E_cleaned.bam.bai

TCGA-L5-A4OR,XX,1,TCGA-L5-A4OR,/home/users/nus/ash.ps/scratch/RNA-mut/RNA-bams/cleaned_bams/TCGA-L5-A4OR_cleaned.bam,/home/users/nus/ash.ps/scratch/RNA-mut/RNA-bams/cleaned_bams/TCGA-L5-A4OR_cleaned.bam.bai
TCGA-L5-A4OR,XX,0,TCGA-L5-A4OR,/home/users/nus/ash.ps/scratch/RNA-mut/RNA-bams/cleaned_bams/TCGA-L5-A4OR-Normal_cleaned.bam,/home/users/nus/ash.ps/scratch/RNA-mut/RNA-bams/cleaned_bams/TCGA-L5-A4OR-Normal_cleaned.bam.bai

patient,sex,status,sample,bam,bai
TCGA-IG-A7DP,XX,1,TCGA-IG-A7DP,/home/users/nus/ash.ps/scratch/RNA-mut/RNA/bams/cleaned_bams/TCGA-IG-A7DP_cleaned.bam,/home/users/nus/ash.ps/scratch/RNA-mut/RNA/bams/cleaned_bams/TCGA-IG-A7DP_cleaned.bam.bai



#!/bin/bash

#PBS -l select=3:ncpus=64:mem=128g
#PBS -l walltime=24:00:00
#PBS -P 11003581
#PBS -N RNA-mut-sarek
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load nextflow/24.10.5
module load miniforge3

export NXF_CONDA_CACHEDIR=/home/project/11003581/Tools/nf-conda-cache
export NXF_CONDA_USE_MAMBA=true

WORKDIR=/home/users/nus/ash.ps/scratch/RNA-mut/RNA/sarek/

nextflow run nf-core/sarek -r 3.5.1 \
   -profile mamba \
   --input $WORKDIR/samplesheet.csv \
   --outdir $WORKDIR \
   --step variant_calling \
   --joint_calling \
   --tools mutect2,manta \
   --fasta /home/project/11003581/Ref/sarek-refs/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta \
   --fai /home/project/11003581/Ref/sarek-refs/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta.fai \
   --dict /home/project/11003581/Ref/sarek-refs/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.dict \
   --pon /home/project/11003581/Ref/sarek-refs/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000g_pon.hg38.vcf.gz \
   --pon_tbi /home/project/11003581/Ref/sarek-refs/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000g_pon.hg38.vcf.gz.tbi



