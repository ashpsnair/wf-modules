#!/bin/bash

#PBS -l select=1:ncpus=32
#PBS -l walltime=6:00:00
#PBS -P 11003581
#PBS -N annovar
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR


perl /home/project/11003581/Tools/annovar/table_annovar.pl /home/users/nus/ash.ps/scratch/YS-analysis/mutect2-vcfs/T01_vs_N01.mutect2.filtered.vcf \
    /home/project/11003581/Tools/annovar/humandb/ \
    -buildver hg38 \
    -out /home/users/nus/ash.ps/scratch/YS-analysis/vcf-to-maf/P01 \
    -protocol refGeneWithVer \
    -operation g \
    -remove -polish -vcfinput -nastring .


/home/users/nus/ash.ps/scratch/YS-analysis/vcf-to-maf/annovar/P01.hg19_multianno.vcf


setwd('/home/users/nus/ash.ps/scratch/YS-analysis/vcf-to-maf/')
library(maftools)

annovarToMaf('/home/users/nus/ash.ps/scratch/YS-analysis/vcf-to-maf/annovar/P01.hg19_multianno.txt', Center = NULL, refBuild = 'hg19', tsbCol = NULL,
table = 'refGene',basename = '/home/users/nus/ash.ps/scratch/YS-analysis/vcf-to-maf/P01.output_maf')
