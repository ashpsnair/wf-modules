#!/bin/bash

#PBS -l select=1:ncpus=10
#PBS -l walltime=10:00:00
#PBS -P 11003581
#PBS -N process-p3l-lab
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load miniforge3
conda activate myenv

#sample details
ref_fasta="/home/project/11003581/Ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
p3l_bam="/home/users/nus/ash.ps/scratch/P3L-lab-analysis/sorted_bams/p3292l_hac_sorted.bam"
wt_bam="/home/users/nus/ash.ps/scratch/P3L-lab-analysis/sorted_bams/WT_MCF10A_hac_sorted.bam"

cnvpytor -root normal_sample.pytor -rd $wt_bam -T $ref_fasta
cnvpytor -root tumor_sample.pytor -rd $p3l_bam -T $ref_fasta

######### Running for tumor normal ##########

#generating read depth files
cnvpytor -root normal_sample.pytor -rd $wt_bam -T $ref_fasta
cnvpytor -root tumor_sample.pytor -rd $p3l_bam -T $ref_fasta

# Call CNVs for both samples
cnvpytor -root normal_sample.pytor -call
cnvpytor -root tumor_sample.pytor -call

#export to vcfs
cnvpytor -root tumor_sample.pytor -view > tumor_cnvs.vcf
cnvpytor -root normal_sample.pytor -view > normal_cnvs.vcf
