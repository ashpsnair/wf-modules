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
#PBS -l walltime=24:00:00
#PBS -P 11003581
#PBS -N run-riboseq-singu
#PBS -j oe

cd /scratch/users/nus/ash.ps/JQQ/analysis-july/sing0riboseq

module load java/17.0.6-jdk
module load singularity/3.10.0
module load nextflow/25.04.6

export NXF_SINGULARITY_CACHEDIR=/home/project/11003581/Tools/singularity-cache/
export SINGULARITY_CACHEDIR=/home/project/11003581/Tools/singularity-cache/

WORKDIR=/scratch/users/nus/ash.ps/JQQ/analysis-july/sing0riboseq
GTF=/home/users/nus/ash.ps/scratch/refs/Homo_sapiens.GRCh38.114.gtf.gz
FASTA=/home/users/nus/ash.ps/scratch/refs/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

nextflow run nf-core/riboseq -r 1.1.0 \
   -profile singularity \
   -c $WORKDIR/hpc.config \
   --input $WORKDIR/samplesheet.csv \
   --contrasts $WORKDIR/contrasts.csv \
   --outdir $WORKDIR \
   --gtf $GTF \
   --fasta $FASTA \
   --remove_ribo_rna true


