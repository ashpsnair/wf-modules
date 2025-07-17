## Samplesheet

sample,fastq_1,strandedness
RNA_0,/home/users/nus/ash.ps/scratch/JQQ/JQQ-rRNA+blkRNA/RNA_0/RNA_0_FKDL250009337-1A_HNNCCDRX5_L2.fq.gz,auto
RNA_21,/home/users/nus/ash.ps/scratch/JQQ/JQQ-rRNA+blkRNA/RNA_21/RNA_21_FKDL250009340-1A_HNNCCDRX5_L2.fq.gz,auto
RNA_3,/home/users/nus/ash.ps/scratch/JQQ/JQQ-rRNA+blkRNA/RNA_3/RNA_3_FKDL250009338-1A_HNNCCDRX5_L2.fq.gz,auto
RNA_5,/home/users/nus/ash.ps/scratch/JQQ/JQQ-rRNA+blkRNA/RNA_5/RNA_5_FKDL250009339-1A_HNNCCDRX5_L2.fq.gz,auto
RNA_Rem5,/home/users/nus/ash.ps/scratch/JQQ/JQQ-rRNA+blkRNA/RNA_Rem5/RNA_Rem5_FKDL250009341-1A_HLVKJDRX5_L1.fq.gz,auto
RNA_Rem5,/home/users/nus/ash.ps/scratch/JQQ/JQQ-rRNA+blkRNA/RNA_Rem5/RNA_Rem5_FKDL250009341-1A_HLVNHDRX5_L1.fq.gz,auto
RNA_Rem5,/home/users/nus/ash.ps/scratch/JQQ/JQQ-rRNA+blkRNA/RNA_Rem5/RNA_Rem5_FKDL250009341-1A_HNNCCDRX5_L2.fq.gz,auto

#!/bin/bash

#PBS -l select=3:ncpus=64:mem=256g
#PBS -l walltime=24:00:00
#PBS -P 11003581
#PBS -N jqq-rnaseq
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load java/17.0.6-jdk
module load singularity

export NXF_SINGULARITY_CACHEDIR=/home/project/11003581/Tools/singularity-cache/
export SINGULARITY_CACHEDIR=/home/project/11003581/Tools/singularity-cache/

/home/project/11003581/Tools/nextflow run nf-core/rnaseq \
    -profile singularity \
    --input /home/users/nus/ash.ps/scratch/JQQ/analysis/rnaseq/samplesheet.csv \
    --outdir /home/users/nus/ash.ps/scratch/JQQ/analysis/rnaseq/ \
    --genome GRCh38 \
    --trimmer fastp \
    --remove_ribo_rna \
    --aligner star_salmon



    