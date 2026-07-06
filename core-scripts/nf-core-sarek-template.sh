#!/bin/bash

#PBS -l select=1:ncpus=64:mem=128g
#PBS -l walltime=24:00:00
#PBS -P 11003581
#PBS -N <PROJECT_NAME>
#PBS -M ash.ps@nus.edu.sg
#PBS -j oe

cd $PBS_O_WORKDIR
module purge
module load java/17.0.6-jdk
module load nextflow/25.10.4
module load singularity

export NXF_HOME=/home/project/11003581/Tools/nf-home
export NXF_SINGULARITY_CACHEDIR=/home/users/nus/ash.ps/scratch/sing/cache/
export SINGULARITY_CACHEDIR=/home/users/nus/ash.ps/scratch/sing/cache/

WORKDIR=/home/users/nus/ash.ps/scratch/sep-fun/nf-sarek/<PROJECT_NAME>

nextflow run nf-core/sarek -r 3.9.0 \
   -profile singularity \
   --input $WORKDIR/samplesheet.csv \
   --outdir $WORKDIR \
   --tools mutect2,strelka,ascat,manta \
   -params-file /home/project/11003581/Tools/nf-sarek-params.yamls
