'''
module load miniforge3
conda create --prefix /home/project/11003581/conda-envs/mimseq python=3.7


conda activate /home/project/11003581/conda-envs/mimseq


conda install -c bioconda mimseq

#check version
conda activate /home/project/11003581/conda-envs/mimseq
mimseq --version

#install cutadapt within the env

'''
####################################################################
#Code for demultiplexing the samples using cutadapt (folder-wise)
####################################################################
### Bracode File Format for cutadapt

#barcodes.fa
>WT  
ATCGTC

>Dox5  
AGCTAC

>Dox20  
GCATAC


#script

#!/bin/bash

#PBS -l select=1:ngpus=1:ncpus=64
#PBS -l walltime=15:00:00
#PBS -P 11003581
#PBS -N demultiplex
#PBS -j oe

module load miniforge3
conda activate /home/project/11003581/conda-envs/mimseq

mkdir -p /home/project/11003581/Data/JQQ-analysis/demultiplexed_fastqs/

cd /home/project/11003581/Data/JQQ/01.RawData/Rep1/

cutadapt --no-indels -q 30,30 --trimmed-only -j 32 -m 10 \
    -a file:/home/project/11003581/Data/JQQ-analysis/barcodes.fa \
    -o /home/project/11003581/Data/JQQ-analysis/demultiplexed_fastqs/{name}.1.fq.gz \
    -p /home/project/11003581/Data/JQQ-analysis/demultiplexed_fastqs/{name}.2.fq.gz \
    Rep1_DKDL240034859-1A_22TKKWLT3_L3_1.fq.gz \
    Rep1_DKDL240034859-1A_22TKKWLT3_L3_2.fq.gz \
    1> /home/project/11003581/Data/JQQ-analysis/demultiplexed_fastqs/log.txt



####################################################################
#Pawan's Code
####################################################################

#!/bin/bash

#PBS -l select=1:ngpus=1:ncpus=64
#PBS -l walltime=15:00:00
#PBS -P 11003581
#PBS -N demultiplex
#PBS -j oe

module load miniforge3
conda activate /home/project/11003581/conda-envs/mimseq

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

for i in /home/project/11003581/Data/JQQ/01.RawData/Rep3_2/Rep3_2_DKDL240034905-1A_22TKKWLT3_L3_1.fq.gz

do fn="Rep1"

cutadapt --no-indels -q 30,30 --trimmed-only -j 10 -a file:/home/project/11003581/Data/JQQ-analysis/barcodes.fa \
    -m 10 -o $fn'_{name}_trim.fastq.gz' $i 1> $fn'_log.txt'

done


### Automate pawan's code

#!/bin/bash

#PBS -l select=1:ngpus=1:ncpus=64
#PBS -l walltime=15:00:00
#PBS -P 11003581
#PBS -N auto-demultiplex
#PBS -j oe

module load miniforge3
conda activate /home/project/11003581/conda-envs/mimseq

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

# Define the directory containing the .fq.gz files
directory="/home/project/11003581/Data/JQQ/01.RawData/*/*.fq.gz"

# Loop through each .fq.gz file in the directory
for i in $directory; do
    # Extract the base filename by removing .fq.gz
    fn=$(echo "$i" | awk -F'/' '{print $NF}' | awk -F'.fq.gz' '{print $1}')

    # Run cutadapt with the specified parameters
    cutadapt --no-indels -q 30,30 --trimmed-only -j 16 -a file:/home/project/11003581/Data/JQQ-analysis/barcodes.fa \
        -m 10 -o ${fn}'_{name}_trim.fastq.gz' $i 1> ${fn}'_log.txt'
done




