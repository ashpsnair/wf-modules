#!/bin/bash

#PBS -l select=3:ncpus=64:mem=256g
#PBS -l walltime=24:00:00
#PBS -P 11003581
#PBS -N jqq-ira-ribo
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module unload java
module load java/17.0.6-jdk
module load singularity/3.10.0
module load nextflow/25.04.6
module load miniforge3

export CACHE=/home/project/11003581/Tools/nf-conda-cache
export PATH=$CACHE/bin:$PATH                   # makes micromamba discoverable
export MAMBA_ROOT_PREFIX=$CACHE/root           # where micromamba keeps “base”
export NXF_MAMBA_EXE=$CACHE/bin/micromamba     # tell Nextflow exactly which mamba
export NXF_CONDA_USE_MAMBA=true
export NXF_CONDA_CACHEDIR=$CACHE/envs          # per-process envs live here
export CONDA_PKGS_DIRS=$CACHE/pkgs             # shared package tarballs


WORKDIR=/home/users/nus/ash.ps/scratch/JQQ/analysis-july/ira-ribo

nextflow run iraiosub/riboseq-flow -r v1.1.1 \
    -profile mamba \
    --skip_umi_extract \
    --contaminants_fasta \
    --input $WORKDIR/samplesheet.csv \
    --outdir $WORKDIR \
    --adapter_threeprime AGATCGGAAGAGCACACGTCTGAACTCCATCAC \
    --org GRCh38 \
    --strandedness forward \
    --feature CDS

    


### Samplesheet
sample,fastq
Ribo_rep1_0d,/home/project/11003581/Data/JQQ-july/Ribo_rep1_0d/Ribo_rep1_0d_FKDL250224978-1A_HYF7GDRX5_L2.fq.gz
Ribo_rep2_0d,/home/project/11003581/Data/JQQ-july/Ribo_rep2_0d/Ribo_rep2_0d_FKDL250220578-1A_HYF7GDRX5_L2.fq.gz
Ribo_rep1_5d,/home/project/11003581/Data/JQQ-july/Ribo_rep1_5d/Ribo_rep1_5d_FKDL250224979-1A_HYF7GDRX5_L2.fq.gz
Ribo_rep2_5d,/home/project/11003581/Data/JQQ-july/Ribo_rep2_5d/Ribo_rep2_5d_FKDL250220579-1A_HYF7GDRX5_L2.fq.gz
Ribo_rep1_20d,/home/project/11003581/Data/JQQ-july/Ribo_rep1_20d/Ribo_rep1_20d_FKDL250224980-1A_HYF7GDRX5_L2.fq.gz
Ribo_rep2_20d,/home/project/11003581/Data/JQQ-july/Ribo_rep2_20d/Ribo_rep2_20d_FKDL250220580-1A_HYF7GDRX5_L2.fq.gz

