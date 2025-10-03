#samplesheet
/home/users/nus/ash.ps/scratch/JQQ/tRNA-Sep/Data/X401SC25086193-Z01-F002/01.RawData/tRF_1/tRF_1_FKDL250259013-1A_HWWTCDRX5_L1.fq.gz tRNA1
/home/users/nus/ash.ps/scratch/JQQ/tRNA-Sep/Data/X401SC25086193-Z01-F002/01.RawData/tRF_3/tRF_3_FKDL250259015-1A_HWWTCDRX5_L1.fq.gz tRNA2
/home/users/nus/ash.ps/scratch/JQQ/tRNA-Sep/Data/X401SC25086193-Z01-F002/01.RawData/tRF_2/tRF_2_FKDL250259014-1A_HWWTCDRX5_L1.fq.gz tRNA3
#note that it has to be tab separated


###############################################################################################
################### MIM SEQ script  #######################
###############################################################################################

#!/bin/bash

#PBS -l select=1:ncpus=64:mem=128g
#PBS -l walltime=05:00:00
#PBS -P 11003581
#PBS -N run-mimseq
#PBS -j oe

module load miniforge3
conda activate /home/project/11003581/conda-envs/mimseq
module load r/4.2.0

export PATH="/home/project/11003581/Tools:$PATH"

cd /home/users/nus/ash.ps/scratch/JQQ/tRNA-Sep/analysis

mimseq --species Hsap --cluster-id 0.97 --threads 64 --min-cov 0.0005 \
    --max-mismatches 0.075 --control-condition tRNA1 \
    -n trna-compare --out-dir trna-compare \
    --max-multi 4 --remap --remap-mismatches 0.05 samplesheet.txt
