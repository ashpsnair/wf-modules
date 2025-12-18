#!/bin/bash

#PBS -l select=1:ncpus=128:mem=256g
#PBS -l walltime=24:00:00
#PBS -P 11003581
#PBS -N run-GEJ09
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load nextflow/24.10.5
module load miniforge3

export NXF_CONDA_CACHEDIR=/home/project/11003581/Tools/nf-conda-cache
export NXF_CONDA_USE_MAMBA=true

WORKDIR=/home/users/nus/ash.ps/scratch/YS-organoid/GEJ09

nextflow run nf-core/sarek -r 3.5.1 \
   -profile conda \
   --input $WORKDIR/samplesheet.csv \
   --outdir $WORKDIR \
   --genome GATK.GRCh38 \
   --tools mutect2,strelka,manta,cnvkit \
   --joint_mutect2 \
   --pon /home/project/11003581/Ref/pons/somatic-hg38_1000g_pon.hg38.vcf.gz \
   --pon_tbi /home/project/11003581/Ref/pons/somatic-hg38_1000g_pon.hg38.vcf.gz.tbi
