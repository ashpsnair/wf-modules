#!/bin/bash

#PBS -l select=1:ncpus=128:mem=256g
#PBS -l walltime=24:00:00
#PBS -P 11003581
#PBS -N vc-RAD51_S181P
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load nextflow/24.10.5
module load miniforge3

export NXF_CONDA_CACHEDIR=/home/project/11003581/Tools/nf-conda-cache
export NXF_CONDA_USE_MAMBA=true

WORKDIR=/home/users/nus/ash.ps/scratch/sep-fun-analysis/RAD51_S181P

nextflow run nf-core/sarek -r 3.5.1 \
   -profile conda \
   --input $WORKDIR/csv/recalibrated.csv \
   --outdir $WORKDIR \
   --genome GATK.GRCh38 \
   --step variant_calling \
   --tools mutect2,cnvkit,manta,tiddit \
   --joint_mutect2 \
   -resume
