#!/bin/bash

#PBS -l select=1:ngpus=1
#PBS -l walltime=10:00:00
#PBS -P personal-ash.ps
#PBS -N cellranger-trial
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR


/home/project/11003581/Tools/cellranger-arc-2.0.2/bin/cellranger-arc count --id=10k_pbmc \
                       --reference=/home/project/11003581/Ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                       --libraries=/home/users/nus/ash.ps/scratch/mulitomics/analysis/10k_PBMC_Multiome_nextgem_Chromium_X_library.csv \
                       --localcores=128 \
                       --localmem=128
