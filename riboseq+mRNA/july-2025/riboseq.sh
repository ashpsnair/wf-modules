#samplesheet
sample,fastq_1,strandedness,type,treatment,pair
Ribo_rep1_0d,/home/users/nus/ash.ps/scratch/JQQ/JQQ-ribsome-profile-july2025/Ribo_rep1_0d/Ribo_rep1_0d_FKDL250224978-1A_HYF7GDRX5_L2.fq.gz,auto,riboseq,control,1
Ribo_rep2_0d,/home/users/nus/ash.ps/scratch/JQQ/JQQ-ribsome-profile-july2025/Ribo_rep2_0d/Ribo_rep2_0d_FKDL250220578-1A_HYF7GDRX5_L2.fq.gz,auto,riboseq,control,2
Ribo_rep1_5d,/home/users/nus/ash.ps/scratch/JQQ/JQQ-ribsome-profile-july2025/Ribo_rep1_5d/Ribo_rep1_5d_FKDL250224979-1A_HYF7GDRX5_L2.fq.gz,auto,riboseq,day5,3
Ribo_rep2_5d,/home/users/nus/ash.ps/scratch/JQQ/JQQ-ribsome-profile-july2025/Ribo_rep2_5d/Ribo_rep2_5d_FKDL250220579-1A_HYF7GDRX5_L2.fq.gz,auto,riboseq,day5,4
Ribo_rep1_20d,/home/users/nus/ash.ps/scratch/JQQ/JQQ-ribsome-profile-july2025/Ribo_rep1_20d/Ribo_rep1_20d_FKDL250224980-1A_HYF7GDRX5_L2.fq.gz,auto,riboseq,day20,5
Ribo_rep2_20d,/home/users/nus/ash.ps/scratch/JQQ/JQQ-ribsome-profile-july2025/Ribo_rep2_20d/Ribo_rep2_20d_FKDL250220580-1A_HYF7GDRX5_L2.fq.gz,auto,riboseq,day20,6
RNA_rep1_0d,/home/users/nus/ash.ps/scratch/JQQ/JQQ-ribsome-profile-july2025/RNA_rep1_0d/RNA_rep1_0d_FKDL250224972-1A_HYF7GDRX5_L1.fq.gz,auto,rnaseq,control,1
RNA_rep2_0d,/home/users/nus/ash.ps/scratch/JQQ/JQQ-ribsome-profile-july2025/RNA_rep2_0d/RNA_rep2_0d_FKDL250224975-1A_HYF7GDRX5_L1.fq.gz,auto,rnaseq,control,2
RNA_rep1_5d,/home/users/nus/ash.ps/scratch/JQQ/JQQ-ribsome-profile-july2025/RNA_rep1_5d/RNA_rep1_5d_FKDL250224973-1A_HYF7GDRX5_L1.fq.gz,auto,rnaseq,day5,3
RNA_rep2_5d,/home/users/nus/ash.ps/scratch/JQQ/JQQ-ribsome-profile-july2025/RNA_rep2_5d/RNA_rep2_5d_FKDL250224976-1A_HYF7GDRX5_L1.fq.gz,auto,rnaseq,day5,4
RNA_rep1_20d,/home/users/nus/ash.ps/scratch/JQQ/JQQ-ribsome-profile-july2025/RNA_rep1_20d/RNA_rep1_20d_FKDL250224974-1A_HYF7GDRX5_L1.fq.gz,auto,rnaseq,day20,5
RNA_rep2_20d,/home/users/nus/ash.ps/scratch/JQQ/JQQ-ribsome-profile-july2025/RNA_rep2_20d/RNA_rep2_20d_FKDL250224977-1A_HYF7GDRX5_L1.fq.gz,auto,rnaseq,day20,6


### contrasts
id,variable,reference,target
day5_vs_day0,treatment,control,day5
day20_vs_day0,treatment,control,day20


#script

#!/bin/bash

#PBS -l select=3:ncpus=64:mem=256g
#PBS -l walltime=10:00:00
#PBS -P 11003581
#PBS -N run-riboseq
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load java/17.0.6-jdk
module load miniforge3

export CACHE=/home/project/11003581/Tools/nf-conda-cache
export PATH=$CACHE/bin:$PATH                   # makes micromamba discoverable
export MAMBA_ROOT_PREFIX=$CACHE/root           # where micromamba keeps “base”
export NXF_MAMBA_EXE=$CACHE/bin/micromamba     # tell Nextflow exactly which mamba
export NXF_CONDA_USE_MAMBA=true
export NXF_CONDA_CACHEDIR=$CACHE/envs          # per-process envs live here
export CONDA_PKGS_DIRS=$CACHE/pkgs             # shared package tarballs

WORKDIR=/home/users/nus/ash.ps/scratch/JQQ/analysis-july/riboseq-nf-core
GTF=/home/project/11003581/Ref/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf
FASTA=/home/project/11003581/Ref/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa

/home/project/11003581/Tools/nextflow.24.10.6 run nf-core/riboseq \
   -profile conda \
   --input $WORKDIR/samplesheet.csv \
   --contrasts $WORKDIR/contrasts.csv \
   --outdir $WORKDIR \
   --gtf $GTF \
   --fasta $FASTA \
   --remove_ribo_rna true \
   --extra_featurecounts_args "-t CDS" \
   --email ash.ps@nus.edu.sg


