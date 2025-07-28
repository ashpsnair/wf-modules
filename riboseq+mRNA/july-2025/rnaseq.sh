#!/bin/bash

#PBS -l select=1:ncpus=128:mem=256g
#PBS -l walltime=24:00:00
#PBS -P 11003581
#PBS -N run-rnaseq
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

WORKDIR=/home/users/nus/ash.ps/scratch/JQQ/analysis-july/rnaseq

/home/project/11003581/Tools/nextflow run nf-core/rnaseq \
    -profile singularity \
    --input $WORKDIR/samplesheet.csv \
    --outdir $WORKDIR \
    --genome GRCh38 \
    --trimmer fastp \
    --remove_ribo_rna \
    --aligner star_salmon


## Samplesheet

sample,fastq_1,strandedness
RNA_rep1_0d,/home/project/11003581/Data/JQQ-july/RNA_rep1_0d/RNA_rep1_0d_FKDL250224972-1A_HYF7GDRX5_L1.fq.gz,forward
RNA_rep2_0d,/home/project/11003581/Data/JQQ-july/RNA_rep2_0d/RNA_rep2_0d_FKDL250224975-1A_HYF7GDRX5_L1.fq.gz,forward
RNA_rep1_5d,/home/project/11003581/Data/JQQ-july/RNA_rep1_5d/RNA_rep1_5d_FKDL250224973-1A_HYF7GDRX5_L1.fq.gz,forward
RNA_rep2_5d,/home/project/11003581/Data/JQQ-july/RNA_rep2_5d/RNA_rep2_5d_FKDL250224976-1A_HYF7GDRX5_L1.fq.gz,forward
RNA_rep1_20d,/home/project/11003581/Data/JQQ-july/RNA_rep1_20d/RNA_rep1_20d_FKDL250224974-1A_HYF7GDRX5_L1.fq.gz,forward
RNA_rep2_20d,/home/project/11003581/Data/JQQ-july/RNA_rep2_20d/RNA_rep2_20d_FKDL250224977-1A_HYF7GDRX5_L1.fq.gz,forward

    