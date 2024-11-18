#!/bin/bash

#PBS -l select=1:ngpus=2
#PBS -l walltime=10:00:00
#PBS -P 11003581
#PBS -N process-p3l-lab
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load python
module load miniforge3
conda activate /home/project/11003581/conda-envs/spectre

##Tools Directory
dorado_bin_path="/home/project/11003581/Tools/dorado-0.7.3-linux-x64/bin"
minimap2_path="/home/project/11003581/Tools/minimap2-2.28_x64-linux/"

#sample details
p3l_pod5_dir="/home/users/nus/ash.ps/scratch/MCF10A-lab/pod5/P3L/"

#ref genome
ref_fasta="/home/project/11003581/Ref/hg38.fa.gz"

#analysis directory
output_dir="/home/users/nus/ash.ps/scratch/re-P3L/"

####### Running dorado
mkdir -p $output_dir/dorado_output

#running dorado
$dorado_bin_path/dorado basecaller hac $p3l_pod5_dir \
        --trim all \
        --reference $ref_fasta > $output_dir/dorado_output/p3292l_hac.bam

#Sorting & Indexing
mkdir $output_dir/sorted_bams
samtools sort $output_dir/dorado_output/p3292l_hac.bam -o $output_dir/sorted_bams/p3292l_hac_sorted.bam
samtools index $output_dir/sorted_bams/p3292l_hac_sorted.bam

#Running mosdepth
mkdir $output_dir/mosdepth

mosdepth -t 8 -x -b 1000 -Q 20 $output_dir/mosdepth $output_dir/sorted_bams/p3292l_hac_sorted.bam

'''
#########bgzip the fasta
module load bcftools

'''

#Running spectre
mkdir -p $output_dir/spectre 

/home/project/11003581/conda-envs/spectre/bin/spectre CNVCaller \
  --coverage $output_dir/mosdepth/p3l.regions.bed.gz \
  --sample-id p3l \
  --only-chr chr13 \
  --output-dir $output_dir/spectre \
  --reference /home/project/11003581/Ref/hg38.fa.gz

