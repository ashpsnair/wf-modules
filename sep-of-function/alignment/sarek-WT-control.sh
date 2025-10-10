
patient,sex,status,sample,lane,fastq_1,fastq_2
WT,XX,1,WT_P30_B,1,/home/project/11003581/Data/sep-function/b1-fastqs/WT-control/WT_P30_B/WT_P30_B_EDHG210027024-1a_HVCKJDSX2_L1_1.fq.gz,/home/project/11003581/Data/sep-function/b1-fastqs/WT-control/WT_P30_B/WT_P30_B_EDHG210027024-1a_HVCKJDSX2_L1_2.fq.gz


#!/bin/bash

#PBS -l select=1:ncpus=128:mem=256g
#PBS -l walltime=10:00:00
#PBS -P 11003581
#PBS -N run-WT-control
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load nextflow/24.10.5
module load miniforge3

export NXF_CONDA_CACHEDIR=/home/project/11003581/Tools/nf-conda-cache
export NXF_CONDA_USE_MAMBA=true

WORKDIR=/home/users/nus/ash.ps/scratch/sep-fun-analysis/CONTROL

nextflow run nf-core/sarek -r 3.5.1 \
   -profile mamba \
   --input $WORKDIR/samplesheet.csv \
   --outdir $WORKDIR