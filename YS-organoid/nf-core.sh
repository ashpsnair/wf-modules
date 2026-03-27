#samplesheet.csv
patient,sex,status,sample,lane,fastq_1,fastq_2
GEJ07,XX,0,GEJ07_Wk0_0,1,/home/project/11003581/kaijing/rawdata/YS_project/merged_rawdata/GEJ07_Wk0_0_R1.fq.gz,/home/project/11003581/kaijing/rawdata/YS_project/merged_rawdata/GEJ07_Wk0_0_R2.fq.gz
GEJ07,XX,1,GEJ07_Wk4_0,1,/home/project/11003581/kaijing/rawdata/YS_project/merged_rawdata/GEJ07_Wk4_0_R1.fq.gz,/home/project/11003581/kaijing/rawdata/YS_project/merged_rawdata/GEJ07_Wk4_0_R2.fq.gz
GEJ07,XX,1,GEJ07_Wk4_2,1,/home/project/11003581/kaijing/rawdata/YS_project/merged_rawdata/GEJ07_Wk4_2_R1.fq.gz,/home/project/11003581/kaijing/rawdata/YS_project/merged_rawdata/GEJ07_Wk4_2_R2.fq.gz
GEJ07,XX,1,GEJ07_Wk4_5,1,/home/project/11003581/kaijing/rawdata/YS_project/merged_rawdata/GEJ07_Wk4_5_R1.fq.gz,/home/project/11003581/kaijing/rawdata/YS_project/merged_rawdata/GEJ07_Wk4_5_R2.fq.gz
GEJ09,XX,0,GEJ09_Wk0_0,1,/home/project/11003581/kaijing/rawdata/YS_project/merged_rawdata/GEJ09_Wk0_0_R1.fq.gz,/home/project/11003581/kaijing/rawdata/YS_project/merged_rawdata/GEJ09_Wk0_0_R2.fq.gz
GEJ09,XX,1,GEJ09_Wk4_0,1,/home/project/11003581/kaijing/rawdata/YS_project/merged_rawdata/GEJ09_Wk4_0_R1.fq.gz,/home/project/11003581/kaijing/rawdata/YS_project/merged_rawdata/GEJ09_Wk4_0_R2.fq.gz
GEJ09,XX,1,GEJ09_Wk4_2,1,/home/project/11003581/kaijing/rawdata/YS_project/merged_rawdata/GEJ09_Wk4_2_R1.fq.gz,/home/project/11003581/kaijing/rawdata/YS_project/merged_rawdata/GEJ09_Wk4_2_R2.fq.gz
GEJ09,XX,1,GEJ09_Wk4_5,1,/home/project/11003581/kaijing/rawdata/YS_project/merged_rawdata/GEJ09_Wk4_5_R1.fq.gz,/home/project/11003581/kaijing/rawdata/YS_project/merged_rawdata/GEJ09_Wk4_5_R2.fq.gz

#!/bin/bash

#PBS -l select=3:ncpus=64:mem=128g
#PBS -l walltime=24:00:00
#PBS -P 11003581
#PBS -N YS-organoid-sing
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load nextflow/24.10.5
module load java/17.0.6-jdk
module load singularity

export NXF_SINGULARITY_CACHEDIR=//home/users/nus/ash.ps/scratch/sing/cache/
export SINGULARITY_CACHEDIR=/home/users/nus/ash.ps/scratch/sing/cache/

WORKDIR=/home/users/nus/ash.ps/scratch/sing/

nextflow run nf-core/sarek -r 3.5.1 \
   -profile singularity \
   --input $WORKDIR/samplesheet.csv \
   --outdir $WORKDIR \
   --genome GATK.GRCh38 \
   --joint_calling \
   --tools mutect2,strelka,mpileup,manta,ascat \
   --pon /home/project/11003581/Ref/sarek-refs/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000g_pon.hg38.vcf.gz \
   --pon_tbi /home/project/11003581/Ref/sarek-refs/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000g_pon.hg38.vcf.gz.tbi

