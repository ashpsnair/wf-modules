cd /home/project/11003581/Tools/
wget https://download.oracle.com/java/23/latest/jdk-23_linux-x64_bin.tar.gz

tar -xzf jdk-23_linux-x64_bin.tar.gz

export JAVA_HOME=/home/project/11003581/Tools/jdk-23.0.1
export PATH=$JAVA_HOME/bin:$PATH



#!/bin/bash

#PBS -l select=1:ncpus=64:mem=256g
#PBS -l walltime=6:00:00
#PBS -P 11003581
#PBS -N YS8-mutect2
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load gcc
module load python/3.12.1-gcc11
export JAVA_HOME=/home/project/11003581/Tools/jdk-23.0.1
export PATH=$JAVA_HOME/bin:$PATH

gatk Mutect2 \
  -R /home/project/11003581/Ref/Homo_sapiens/GATK/hg38/Homo_sapiens_assembly38.fasta \
  -I /home/project/11003581/Data/YS-analysis/crams/N08/N08.recal.cram \
  -I /home/project/11003581/Data/YS-analysis/crams/T08/T08.recal.cram \
  --normal-sample patient08_N08 \
  --output /home/users/nus/ash.ps/scratch/YS-analysis/gatk-mutect2/T08_vs_N08.vcf.gz \