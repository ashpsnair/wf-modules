#!/bin/bash

#PBS -l select=1:ncpus=128:mem=256g
#PBS -l walltime=24:00:00
#PBS -P 11003581
#PBS -N run-riboseq-mamba
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module unload java
module load java/17.0.6-jdk
module load singularity/3.10.0
module load nextflow/25.04.6
module load miniforge3

export CACHE=/home/project/11003581/Tools/nf-conda-cache
export MAMBA_NO_CONFIRM=1
export NXF_MAMBA_CLI_OPTS='-y'
export NXF_CONDA_CLI_OPTS='-y'
export PATH=$CACHE/bin:$PATH                   # makes micromamba discoverable
export MAMBA_ROOT_PREFIX=$CACHE/root           # where micromamba keeps “base”
export NXF_MAMBA_EXE=$CACHE/bin/micromamba     # tell Nextflow exactly which mamba
export NXF_CONDA_USE_MAMBA=true
export NXF_CONDA_CACHEDIR=$CACHE/envs          # per-process envs live here
export CONDA_PKGS_DIRS=$CACHE/pkgs             # shared package tarballs

WORKDIR=/home/users/nus/ash.ps/scratch/JQQ/analysis-july/re-sent/analysis/
GTF=/home/users/nus/ash.ps/scratch/refs/Homo_sapiens.GRCh38.114.gtf
FASTA=/home/users/nus/ash.ps/scratch/refs/Homo_sapiens.GRCh38.dna.primary_assembly.fa

nextflow run nf-core/riboseq -r dev \
   -profile mamba \
   --input $WORKDIR/samplesheet.csv \
   --outdir $WORKDIR \
   --gtf $GTF \
   --fasta $FASTA \
   --skip_ribotricer true \
   --remove_ribo_rna true \
   --email ash.ps@nus.edu.sg



### samplesheet.csv
sample,fastq_1,fastq_2,strandedness,type,treatment
ribo_1,/home/users/nus/ash.ps/scratch/JQQ/analysis-july/re-sent/data/X401SC25086193-Z01-F001/01.RawData/ribo_1/ribo_1_FKDL250251250-1A_HW772DRX5_L1.fq.gz,,forward,riboseq
ribo_2,/home/users/nus/ash.ps/scratch/JQQ/analysis-july/re-sent/data/X401SC25086193-Z01-F001/01.RawData/ribo_2/ribo_2_FKDL250251251-1A_HW772DRX5_L1.fq.gz,,forward,riboseq
ribo_2,/home/users/nus/ash.ps/scratch/JQQ/analysis-july/re-sent/data/X401SC25086193-Z01-F001/01.RawData/ribo_2/ribo_2_FKDL250251251-1A_HW7CVDRX5_L2.fq.gz,,forward,riboseq
ribo_3,/home/users/nus/ash.ps/scratch/JQQ/analysis-july/re-sent/data/X401SC25086193-Z01-F001/01.RawData/ribo_3/ribo_3_FKDL250251252-1A_HW772DRX5_L1.fq.gz,,forward,riboseq
ribo_3,/home/users/nus/ash.ps/scratch/JQQ/analysis-july/re-sent/data/X401SC25086193-Z01-F001/01.RawData/ribo_3/ribo_3_FKDL250251252-1A_HW7CVDRX5_L2.fq.gz,,forward,riboseq
ribo_4,/home/users/nus/ash.ps/scratch/JQQ/analysis-july/re-sent/data/X401SC25086193-Z01-F001/01.RawData/ribo_4/ribo_4_FKDL250251253-1A_HW772DRX5_L1.fq.gz,,forward,riboseq
ribo_5,/home/users/nus/ash.ps/scratch/JQQ/analysis-july/re-sent/data/X401SC25086193-Z01-F001/01.RawData/ribo_5/ribo_5_FKDL250251254-1A_HW772DRX5_L1.fq.gz,,forward,riboseq
ribo_6,/home/users/nus/ash.ps/scratch/JQQ/analysis-july/re-sent/data/X401SC25086193-Z01-F001/01.RawData/ribo_6/ribo_6_FKDL250251255-1A_HW772DRX5_L1.fq.gz,,forward,riboseq


