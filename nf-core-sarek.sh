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

/home/project/11003581/Tools/nextflow run nf-core/sarek \
   -profile singularity \
   --input /home/project/11003581/Data/YS-analysis/samplesheet.csv \
   --outdir /home/project/11003581/Data/YS-analysis/ \
   -resume
