module load singularity/4.3.1
module load nextflow/25.04.6
module unload java
module load java/17.0.6-jdk

#Set memory caps in your shell (recommended by nf-core):
export NXF_OPTS='-Xms1g -Xmx4g'

export SINGULARITY_BINDPATH=/data/projects/11003581,/home/users/nus/ash.ps,/scratch,/home/project,/home/project/11003581/Tools/singularity-cache

export NXF_SINGULARITY_CACHEDIR=/home/project/11003581/Tools/singularity-cache
export SINGULARITY_CACHEDIR=/home/project/11003581/Tools/singularity-cache



#!/bin/bash

#PBS -l select=3:ncpus=64:mem=256g
#PBS -l walltime=24:00:00
#PBS -P 11003581
#PBS -N run-riboseq
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

nextflow run nf-core/riboseq \
   -profile mamba \
   --input $WORKDIR/samplesheet.csv \
   --contrasts $WORKDIR/contrasts.csv \
   --outdir $WORKDIR \
   --gtf $GTF \
   --fasta $FASTA \
   --remove_ribo_rna true




sample,fastq,condition,replicate
Ribo_rep1_0d,/home/project/11003581/Data/JQQ-july/Ribo_rep1_0d/Ribo_rep1_0d_FKDL250224978-1A_HYF7GDRX5_L2.fq.gz,day0,1


replace ,/home/users/nus/ash.ps/scratch/JQQ/JQQ-ribsome-profile-july2025 with /home/project/11003581/Data/JQQ-july/ in the samplesheet.csv