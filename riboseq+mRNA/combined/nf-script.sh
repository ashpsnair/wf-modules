##samplesheet

sample,fastq_1,fastq_2,strandedness,type,sample_description,pair,treatment
RNA_0,/home/users/nus/ash.ps/scratch/JQQ/JQQ-rRNA+blkRNA/RNA_0/RNA_0_FKDL250009337-1A_HNNCCDRX5_L2.fq.gz,,forward,rnaseq,RNA_0,1,control
RNA_Rem5,/home/users/nus/ash.ps/scratch/JQQ/JQQ-rRNA+blkRNA/RNA_Rem5/RNA_Rem5.fq.gz,,forward,rnaseq,RNA_Rem5,2,control
RNA_3,/home/users/nus/ash.ps/scratch/JQQ/JQQ-rRNA+blkRNA/RNA_3/RNA_3_FKDL250009338-1A_HNNCCDRX5_L2.fq.gz,,forward,rnaseq,RNA_3,3,treated
RNA_5,/home/users/nus/ash.ps/scratch/JQQ/JQQ-rRNA+blkRNA/RNA_5/RNA_5_FKDL250009339-1A_HNNCCDRX5_L2.fq.gz,,forward,rnaseq,RNA_5,4,treated
RNA_21,/home/users/nus/ash.ps/scratch/JQQ/JQQ-rRNA+blkRNA/RNA_21/RNA_21_FKDL250009340-1A_HNNCCDRX5_L2.fq.gz,,forward,rnaseq,RNA_21,5,treated
Ribo_0,/home/users/nus/ash.ps/scratch/JQQ/JQQ-rRNA+blkRNA/Ribo_0/Ribo_0_FKDL250009332-1A_HNNCCDRX5_L2.fq.gz,,forward,riboseq,Ribo_0,1,control
Ribo_Rem5,/home/users/nus/ash.ps/scratch/JQQ/JQQ-rRNA+blkRNA/Ribo_Rem5/Ribo_Rem5_FKDL250009336-1A_HNNCCDRX5_L2.fq.gz,,forward,riboseq,Ribo_Rem5,2,control
Ribo_3,/home/users/nus/ash.ps/scratch/JQQ/JQQ-rRNA+blkRNA/Ribo_3/Ribo_3_FKDL250009333-1A_HNNCCDRX5_L2.fq.gz,,forward,riboseq,Ribo_3,3,treated
Ribo_5,/home/users/nus/ash.ps/scratch/JQQ/JQQ-rRNA+blkRNA/Ribo_5/Ribo_5_FKDL250009334-1A_HNNCCDRX5_L2.fq.gz,,forward,riboseq,Ribo_5,4,treated
Ribo_21,/home/users/nus/ash.ps/scratch/JQQ/JQQ-rRNA+blkRNA/Ribo_21/Ribo_21_FKDL250009335-1A_HNNCCDRX5_L2.fq.gz,,forward,riboseq,Ribo_21,5,treated


#!/bin/bash

#PBS -l select=3:ncpus=64:mem=256g
#PBS -l walltime=24:00:00
#PBS -P 11003581
#PBS -N combined-riboseq
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load java/17.0.6-jdk
module load singularity

export NXF_SINGULARITY_CACHEDIR=/home/project/11003581/Tools/singularity-cache/
export SINGULARITY_CACHEDIR=/home/project/11003581/Tools/singularity-cache/

/home/project/11003581/Tools/nextflow run nf-core/riboseq \
   -profile singularity \
   --input /home/users/nus/ash.ps/scratch/JQQ/analysis/combined-analysis/samplesheet.csv \
   --outdir /home/users/nus/ash.ps/scratch/JQQ/analysis/combined-analysis \
   --gtf /home/project/11003581/Ref/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf \
   --fasta /home/project/11003581/Ref/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa \
   --remove_ribo_rna \
   --email ash.ps@nus.edu.sg

