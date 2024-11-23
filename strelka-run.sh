#!/bin/bash

#PBS -l select=1:ncpus=64
#PBS -l walltime=6:00:00
#PBS -P 11003581
#PBS -N YS8-strelka
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

STRELKA_INSTALL_PATH=/home/project/11003581/Tools/strelka-2.9.2.centos6_x86_64/

# configuration
${STRELKA_INSTALL_PATH}/bin/configureStrelkaSomaticWorkflow.py \
    --normalBam /home/users/nus/ash.ps/scratch/YS-batch2/preprocessing/recalibrated/N08/N08.recal.cram \
    --tumorBam /home/users/nus/ash.ps/scratch/YS-batch2/preprocessing/recalibrated/T08/T08.recal.cram \
    --referenceFasta /home/project/11003581/Ref/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/version0.6.0/genome.fa \
    --runDir YS8_somatic

# execution on a single local machine with 20 parallel jobs
demo_somatic/runWorkflow.py -m local -j 20
