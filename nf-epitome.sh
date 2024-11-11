#!/bin/bash

#PBS -l select=1:ncpus=64:mem=256g
#PBS -l walltime=20:00:00
#PBS -P 11003581
#PBS -N trial-epi2me
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR


module load java/17.0.6-jdk
module load singularity


/home/project/11003581/Tools/nextflow pull epi2me-labs/wf-human-variation

/home/project/11003581/Tools/nextflow run epi2me-labs/wf-human-variation \
	--bam '/home/project/11003581/Data/Ash/P3L-lab-analysis/sorted_bams/p3292l_hac_sorted.bam' \
	--ref '/home/project/11003581/Ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna' \
	--sample_name 'p3l' \
	--cnv \
	-profile singularity