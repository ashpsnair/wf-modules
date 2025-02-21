#!/bin/bash

#PBS -l select=3:ncpus=64:mem=256g
#PBS -l walltime=24:00:00
#PBS -P 11003581
#PBS -N jqq-rna-with-riboseq-ira
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
    --input /home/users/nus/ash.ps/scratch/JQQ/analysis/rna-with-ribo-ira/samplesheet.csv \
    --outdir /home/users/nus/ash.ps/scratch/JQQ/analysis/rna-with-ribo-ira/ \
    --org GRCh38 \
    --strandedness forward \
    --feature CDS



### Samplesheet

sample,fastq
RNA_0,/home/users/nus/ash.ps/scratch/JQQ/JQQ-rRNA+blkRNA/RNA_0/RNA_0_FKDL250009337-1A_HNNCCDRX5_L2.fq.gz
RNA_21,/home/users/nus/ash.ps/scratch/JQQ/JQQ-rRNA+blkRNA/RNA_21/RNA_21_FKDL250009340-1A_HNNCCDRX5_L2.fq.gz
RNA_3,/home/users/nus/ash.ps/scratch/JQQ/JQQ-rRNA+blkRNA/RNA_3/RNA_3_FKDL250009338-1A_HNNCCDRX5_L2.fq.gz
RNA_5,/home/users/nus/ash.ps/scratch/JQQ/JQQ-rRNA+blkRNA/RNA_5/RNA_5_FKDL250009339-1A_HNNCCDRX5_L2.fq.gz
RNA_Rem5,/home/users/nus/ash.ps/scratch/JQQ/JQQ-rRNA+blkRNA/RNA_Rem5/RNA_Rem5_FKDL250009341-1A_HLVKJDRX5_L1.fq.gz
RNA_Rem5,/home/users/nus/ash.ps/scratch/JQQ/JQQ-rRNA+blkRNA/RNA_Rem5/RNA_Rem5_FKDL250009341-1A_HLVNHDRX5_L1.fq.gz
RNA_Rem5,/home/users/nus/ash.ps/scratch/JQQ/JQQ-rRNA+blkRNA/RNA_Rem5/RNA_Rem5_FKDL250009341-1A_HNNCCDRX5_L2.fq.gz
