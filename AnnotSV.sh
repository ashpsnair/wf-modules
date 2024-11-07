#!/bin/bash

#PBS -l select=1:ncpus=64:mem=256g
#PBS -l walltime=15:00:00
#PBS -P 11003581
#PBS -N trial-annotsv
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load gcc
module load python/3.12.1-gcc11

#sample details
vcf_file="/home/project/11003581/Data/Ash/P3L-lab-analysis/sniffles/output_dir/p3l-uniq.vcf"
output_dir="/home/project/11003581/Data/Ash/P3L-lab-analysis/annotsv"

/home/project/11003581/Tools/AnnotSV/bin/AnnotSV -SVinputFile $vcf_file \
    -bcftools /app/apps/bcftools/1.15.1/bin/bcftools \
    -bedtools /app/apps/bedtools/2.30.0/bin/bedtools \
    -annotationMode full -outputDir $output_dir




