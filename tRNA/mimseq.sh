#!/bin/bash

#PBS -l select=1:ngpus=1:ncpus=64
#PBS -l walltime=15:00:00
#PBS -P 11003581
#PBS -N run-mimseq
#PBS -j oe

module load miniforge3
conda activate /home/project/11003581/conda-envs/mimseq
module load r/4.2.0

cd /home/project/11003581/Data/JQQ-analysis

mimseq --species Hsap --cluster-id 0.97 --threads 32 --min-cov 0.0005 \
    --max-mismatches 0.075 --control-condition WT \
    -n WTvsDox --out-dir WTvsDox \
    --max-multi 4 --remap --remap-mismatches 0.05 sample.txt


