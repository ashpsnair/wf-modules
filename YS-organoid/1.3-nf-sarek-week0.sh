patient,sex,status,sample,cram,crai
GEJ07_Wk0_0,XX,0,GEJ07_Wk0_0,/home/users/nus/ash.ps/scratch/sing//preprocessing/recalibrated/GEJ07_Wk0_0/GEJ07_Wk0_0.recal.cram,/home/users/nus/ash.ps/scratch/sing//preprocessing/recalibrated/GEJ07_Wk0_0/GEJ07_Wk0_0.recal.cram.crai
GEJ09_Wk0_0,XX,0,GEJ09_Wk0_0,/home/users/nus/ash.ps/scratch/sing//preprocessing/recalibrated/GEJ09_Wk0_0/GEJ09_Wk0_0.recal.cram,/home/users/nus/ash.ps/scratch/sing//preprocessing/recalibrated/GEJ09_Wk0_0/GEJ09_Wk0_0.recal.cram.crai

#!/bin/bash

#PBS -l select=1:ncpus=128:mem=256g
#PBS -l walltime=12:00:00
#PBS -P 11003581
#PBS -N YS-organoid-vc
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load nextflow/24.10.5
module load java/17.0.6-jdk
module load singularity

export NXF_SINGULARITY_CACHEDIR=/home/users/nus/ash.ps/scratch/YS-organoid-week0/cache/
export SINGULARITY_CACHEDIR=/home/users/nus/ash.ps/scratch/YS-organoid-week0/cache/

WORKDIR=/home/users/nus/ash.ps/scratch/YS-organoid/

nextflow run nf-core/sarek -r 3.5.1 \
   -profile singularity \
   --input $WORKDIR/samplesheet.csv \
   --outdir $WORKDIR \
   --genome GATK.GRCh38 \
   --tools mutect2 \
   --step variant_calling \
   --pon /home/project/11003581/Ref/sarek-refs/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000g_pon.hg38.vcf.gz \
   --pon_tbi /home/project/11003581/Ref/sarek-refs/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000g_pon.hg38.vcf.gz.tbi
