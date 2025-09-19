#!/bin/bash

#PBS -l select=1:ncpus=128:mem=256g
#PBS -l walltime=24:00:00
#PBS -P 11003581
#PBS -N an-conda
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module unload java
module load java/17.0.6-jdk
module load singularity/3.10.0
module load nextflow/25.04.6
module load miniforge3

WORKDIR=/home/users/nus/ash.ps/scratch/JQQ/analysis-july/re-sent/analysis-1/
GTF=/home/users/nus/ash.ps/scratch/refs/Homo_sapiens.GRCh38.114.gtf
FASTA=/home/users/nus/ash.ps/scratch/refs/Homo_sapiens.GRCh38.dna.primary_assembly.fa

nextflow run nf-core/riboseq -r dev \
   -profile conda \
   --input $WORKDIR/samplesheet.csv \
   --outdir $WORKDIR \
   --gtf $GTF \
   --fasta $FASTA \
   --skip_ribotricer true \
   --remove_ribo_rna true
