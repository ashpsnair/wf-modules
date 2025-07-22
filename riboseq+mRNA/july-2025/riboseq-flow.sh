#!/bin/bash

#PBS -l select=3:ncpus=64:mem=256g
#PBS -l walltime=24:00:00
#PBS -P 11003581
#PBS -N jqq-riboseq-ira
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load java/17.0.6-jdk
module load singularity/3.10.0
module load nextflow/22.04.0 

export NXF_SINGULARITY_CACHEDIR=/home/project/11003581/Tools/singularity-cache/
export SINGULARITY_CACHEDIR=/home/project/11003581/Tools/singularity-cache/

WORKDIR=/home/users/nus/ash.ps/scratch/JQQ/analysis-july/ira-ribo
ENSEMBL44=/home/users/nus/ash.ps/scratch/refs

nextflow run iraiosub/riboseq-flow -r v1.1.1 \
    -profile singularity \
    --skip_umi_extract \
    --contaminants_fasta \
    --input $WORKDIR/samplesheet.csv \
    --outdir $WORKDIR \
    --adapter_threeprime AGATCGGAAGAGCACACGTCTGAACTCCATCAC \
    --fasta $ENSEMBL44/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
    --gtf $ENSEMBL44/Homo_sapiens.GRCh38.114.gtf.gz \
    --strandedness forward \
    --feature CDS

    


### Samplesheet

sample,fastq,condition,replicate
Ribo_rep1_0d,/home/users/nus/ash.ps/scratch/JQQ/JQQ-ribsome-profile-july2025/Ribo_rep1_0d/Ribo_rep1_0d_FKDL250224978-1A_HYF7GDRX5_L2.fq.gz,day0,1
Ribo_rep2_0d,/home/users/nus/ash.ps/scratch/JQQ/JQQ-ribsome-profile-july2025/Ribo_rep2_0d/Ribo_rep2_0d_FKDL250220578-1A_HYF7GDRX5_L2.fq.gz,day0,2
Ribo_rep1_5d,/home/users/nus/ash.ps/scratch/JQQ/JQQ-ribsome-profile-july2025/Ribo_rep1_5d/Ribo_rep1_5d_FKDL250224979-1A_HYF7GDRX5_L2.fq.gz,day5,1
Ribo_rep2_5d,/home/users/nus/ash.ps/scratch/JQQ/JQQ-ribsome-profile-july2025/Ribo_rep2_5d/Ribo_rep2_5d_FKDL250220579-1A_HYF7GDRX5_L2.fq.gz,day5,2
Ribo_rep1_20d,/home/users/nus/ash.ps/scratch/JQQ/JQQ-ribsome-profile-july2025/Ribo_rep1_20d/Ribo_rep1_20d_FKDL250224980-1A_HYF7GDRX5_L2.fq.gz,day20,1
Ribo_rep2_20d,/home/users/nus/ash.ps/scratch/JQQ/JQQ-ribsome-profile-july2025/Ribo_rep2_20d/Ribo_rep2_20d_FKDL250220580-1A_HYF7GDRX5_L2.fq.gz,day20,2