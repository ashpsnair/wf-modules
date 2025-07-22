#!/bin/bash

#PBS -l select=3:ncpus=64:mem=256g
#PBS -l walltime=24:00:00
#PBS -P 11003581
#PBS -N run-riboseq
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

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

# ───────── NXF FILESYSTEM TUNING ─────────
export NXF_WORK=/home/users/nus/ash.ps/scratch/JQQ/analysis-july/riboseq-nf-core/work
export NXF_TEMP=/home/users/nus/ash.ps/scratch/JQQ/analysis-july/riboseq-nf-core/.nxf_temp

# Optional: Faster hashing for NFS/GPFS
export NXF_FILE_HASH_ALGO=SHA-1

# Optional: Limit thread pool queue to avoid FS lockups
export NXF_OPTS='-Dnxf.file.transfer.buffer.size=1M -Dnxf.file.transfer.threads=2'

WORKDIR=/home/users/nus/ash.ps/scratch/JQQ/analysis-july/riboseq-nf-core
GTF=/home/users/nus/ash.ps/scratch/refs/Homo_sapiens.GRCh38.114.gtf
FASTA=/home/users/nus/ash.ps/scratch/refs/Homo_sapiens.GRCh38.dna.primary_assembly.fa

nextflow run nf-core/riboseq \
   -profile conda \
   --input $WORKDIR/samplesheet.csv \
   --contrasts $WORKDIR/contrasts.csv \
   --outdir $WORKDIR \
   --gtf $GTF \
   --fasta $FASTA \
   --remove_ribo_rna true \
   --email ash.ps@nus.edu.sg




sample,fastq,condition,replicate
Ribo_rep1_0d,/home/project/11003581/Data/JQQ-july//Ribo_rep1_0d/Ribo_rep1_0d_FKDL250224978-1A_HYF7GDRX5_L2.fq.gz,day0,1
