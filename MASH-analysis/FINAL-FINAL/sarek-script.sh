#!/bin/bash

#PBS -l select=3:ncpus=64:mem=256g
#PBS -l walltime=24:00:00
#PBS -P 11003581
#PBS -N mash-TN-b3
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load nextflow/24.10.5
module load miniforge3

export NXF_CONDA_CACHEDIR=/home/project/11003581/Tools/nf-conda-cache
export NXF_CONDA_USE_MAMBA=true

nextflow run nf-core/sarek -r 3.5.1 \
   -profile conda \
   --input /home/users/nus/enambis/scratch/NCCS-MASH/b3/samplesheet.csv \
   --outdir /home/users/nus/enambis/scratch/NCCS-MASH/b3/ \
   --genome GATK.GRCh38 \
   --tools mutect2,manta,ascat \
   --joint_mutect2 \
   --igenomes_base /home/project/11003581/Ref/sarek-refs/  

##### samplesheet.csv (Example)
patient,sex,status,sample,lane,fastq_1,fastq_2
B007,XX,1,WHT392,01,/home/project/11003581/Data/NCCS-MASH/tumor-fastqs/b1/WHT392_R1.fastq.gz,/home/project/11003581/Data/NCCS-MASH/tumor-fastqs/b1/WHT392_R2.fastq.gz
B007,XX,1,WHT393,01,/home/project/11003581/Data/NCCS-MASH/tumor-fastqs/b1/WHT393_R1.fastq.gz,/home/project/11003581/Data/NCCS-MASH/tumor-fastqs/b1/WHT393_R2.fastq.gz