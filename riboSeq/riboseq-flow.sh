#!/bin/bash

#PBS -l select=3:ncpus=64:mem=256g
#PBS -l walltime=24:00:00
#PBS -P 11003581
#PBS -N jqq-riboseq-ira
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load java/17.0.6-jdk
module load singularity

export NXF_SINGULARITY_CACHEDIR=/home/project/11003581/Tools/singularity-cache/
export SINGULARITY_CACHEDIR=/home/project/11003581/Tools/singularity-cache/


/home/project/11003581/Tools/nextflow run iraiosub/riboseq-flow

nextflow run iraiosub/riboseq-flow -r v1.1.1 \
    -profile singularity \
    --skip_umi_extract \
    --input /home/users/nus/ash.ps/scratch/JQQ/analysis/riboseq-ira/samplesheet.csv \
    --outdir /home/users/nus/ash.ps/scratch/JQQ/analysis/riboseq/ \
    --adapter_threeprime AGATCGGAAGAGCACACGTCTGAACTCCATCAC \
    --org GRCh38 \
    --strandedness forward \
    --feature CDS



### Samplesheet
sample,fastq
Ribo_0,/home/users/nus/ash.ps/scratch/JQQ/JQQ-rRNA+blkRNA/Ribo_0/Ribo_0_FKDL250009332-1A_HNNCCDRX5_L2.fq.gz
Ribo_21,/home/users/nus/ash.ps/scratch/JQQ/JQQ-rRNA+blkRNA/Ribo_21/Ribo_21_FKDL250009335-1A_HNNCCDRX5_L2.fq.gz
Ribo_3,/home/users/nus/ash.ps/scratch/JQQ/JQQ-rRNA+blkRNA/Ribo_3/Ribo_3_FKDL250009333-1A_HNNCCDRX5_L2.fq.gz,auto
Ribo_5,/home/users/nus/ash.ps/scratch/JQQ/JQQ-rRNA+blkRNA/Ribo_5/Ribo_5_FKDL250009334-1A_HNNCCDRX5_L2.fq.gz,auto
Ribo_Rem5,/home/users/nus/ash.ps/scratch/JQQ/JQQ-rRNA+blkRNA/Ribo_Rem5/Ribo_Rem5_FKDL250009336-1A_HNNCCDRX5_L2.fq.gz
