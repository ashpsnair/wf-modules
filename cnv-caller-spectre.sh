#!/bin/bash

#PBS -l select=1:mem=128gb
#PBS -l walltime=1:00:00
#PBS -P 11003581
#PBS -N trial-mosdepth
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load miniforge3
conda activate /home/project/11003581/conda-envs/spectre

#output_dir=/home/project/11003581/Data/Ash/P3L-lab-analysis/mosdepth/p3l
#bamfile=/home/project/11003581/Data/Ash/P3L-lab-analysis/sorted_bams/p3292l_hac_sorted.bam

#mosdepth -t 8 -x -b 1000 -Q 20 $output_dir $bamfile

/home/project/11003581/conda-envs/spectre/bin/spectre CNVCaller \
  --coverage /home/project/11003581/Data/Ash/P3L-lab-analysis/mosdepth/p3l.regions.bed.gz \
  --sample-id p3l \
  --only-chr chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY \
  --output-dir /home/project/11003581/Data/Ash/P3L-lab-analysis/spectre_out/ \
  --reference /home/project/11003581/Ref/hg38-bgzip.fa.gz \
  --cancer


'''
#Installation of spectre

conda create --prefix /home/project/11003581/conda-envs/spectre pip -y
module load gcc
module load python/3.11.7-gcc11
conda install bioconda::mosdepth
pip install spectre-cnv

'''





#running for wile type

#!/bin/bash

#PBS -l select=1:mem=128gb
#PBS -l walltime=1:00:00
#PBS -P 11003581
#PBS -N wt-spectre
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load miniforge3
conda activate /home/project/11003581/conda-envs/spectre


output_dir=/home/project/11003581/Data/Ash/P3L-lab-analysis/mosdepth/wt
bamfile=/home/project/11003581/Data/Ash/P3L-lab-analysis/sorted_bams/WT_MCF10A_hac_sorted.bam

mosdepth -t 8 -x -b 1000 -Q 20 $output_dir $bamfile

/home/project/11003581/conda-envs/spectre/bin/spectre CNVCaller \
  --coverage /home/project/11003581/Data/Ash/P3L-lab-analysis/mosdepth/wt.regions.bed.gz \
  --sample-id wt \
  --only-chr chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY \
  --output-dir /home/project/11003581/Data/Ash/P3L-lab-analysis/spectre_out/ \
  --reference /home/project/11003581/Ref/hg38-bgzip.fa.gz