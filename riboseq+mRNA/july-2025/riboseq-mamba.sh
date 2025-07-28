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
export PATH=$CACHE/bin:$PATH                   # makes micromamba discoverable
export MAMBA_ROOT_PREFIX=$CACHE/root           # where micromamba keeps “base”
export NXF_MAMBA_EXE=$CACHE/bin/micromamba     # tell Nextflow exactly which mamba
export NXF_CONDA_USE_MAMBA=true
export NXF_CONDA_CACHEDIR=$CACHE/envs          # per-process envs live here
export CONDA_PKGS_DIRS=$CACHE/pkgs             # shared package tarballs

WORKDIR=/home/users/nus/ash.ps/scratch/JQQ/analysis-july/mamba
GTF=/home/users/nus/ash.ps/scratch/refs/Homo_sapiens.GRCh38.114.gtf
FASTA=/home/users/nus/ash.ps/scratch/refs/Homo_sapiens.GRCh38.dna.primary_assembly.fa

nextflow run nf-core/riboseq -r dev \
   -profile mamba \
   --input $WORKDIR/samplesheet.csv \
   --contrasts $WORKDIR/contrasts.csv \
   --outdir $WORKDIR \
   --gtf $GTF \
   --fasta $FASTA \
   --skip_ribotricer true \
   --remove_ribo_rna true \
   -resume

### samplesheet.csv
sample,fastq_1,fastq_2,strandedness,type,treatment
Ribo_rep1_0d,/home/project/11003581/Data/JQQ-july/Ribo_rep1_0d/Ribo_rep1_0d_FKDL250224978-1A_HYF7GDRX5_L2.fq.gz,,forward,riboseq,control
Ribo_rep2_0d,/home/project/11003581/Data/JQQ-july/Ribo_rep2_0d/Ribo_rep2_0d_FKDL250220578-1A_HYF7GDRX5_L2.fq.gz,,forward,riboseq,control
Ribo_rep1_5d,/home/project/11003581/Data/JQQ-july/Ribo_rep1_5d/Ribo_rep1_5d_FKDL250224979-1A_HYF7GDRX5_L2.fq.gz,,forward,riboseq,day5
Ribo_rep2_5d,/home/project/11003581/Data/JQQ-july/Ribo_rep2_5d/Ribo_rep2_5d_FKDL250220579-1A_HYF7GDRX5_L2.fq.gz,,forward,riboseq,day5
Ribo_rep1_20d,/home/project/11003581/Data/JQQ-july/Ribo_rep1_20d/Ribo_rep1_20d_FKDL250224980-1A_HYF7GDRX5_L2.fq.gz,,forward,riboseq,day20
Ribo_rep2_20d,/home/project/11003581/Data/JQQ-july/Ribo_rep2_20d/Ribo_rep2_20d_FKDL250220580-1A_HYF7GDRX5_L2.fq.gz,,forward,riboseq,day20
RNA_rep1_0d,/home/project/11003581/Data/JQQ-july/RNA_rep1_0d/RNA_rep1_0d_FKDL250224972-1A_HYF7GDRX5_L1.fq.gz,,forward,rnaseq,control
RNA_rep2_0d,/home/project/11003581/Data/JQQ-july/RNA_rep2_0d/RNA_rep2_0d_FKDL250224975-1A_HYF7GDRX5_L1.fq.gz,,forward,rnaseq,control
RNA_rep1_5d,/home/project/11003581/Data/JQQ-july/RNA_rep1_5d/RNA_rep1_5d_FKDL250224973-1A_HYF7GDRX5_L1.fq.gz,,forward,rnaseq,day5
RNA_rep2_5d,/home/project/11003581/Data/JQQ-july/RNA_rep2_5d/RNA_rep2_5d_FKDL250224976-1A_HYF7GDRX5_L1.fq.gz,,forward,rnaseq,day5
RNA_rep1_20d,/home/project/11003581/Data/JQQ-july/RNA_rep1_20d/RNA_rep1_20d_FKDL250224974-1A_HYF7GDRX5_L1.fq.gz,,forward,rnaseq,day20
RNA_rep2_20d,/home/project/11003581/Data/JQQ-july/RNA_rep2_20d/RNA_rep2_20d_FKDL250224977-1A_HYF7GDRX5_L1.fq.gz,,forward,rnaseq,day20

##contrasts.csv
id,variable,reference,target
day5_vs_day0,treatment,control,day5
day20_vs_day0,treatment,control,day20
