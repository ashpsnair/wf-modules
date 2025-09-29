'''
module load miniforge3
conda create --prefix /home/project/11003581/conda-envs/mimseq python=3.7


conda activate /home/project/11003581/conda-envs/mimseq


conda install -c bioconda mimseq

#check version
conda activate /home/project/11003581/conda-envs/mimseq
mimseq --version

#install cutadapt within the env

#install usearch locally (v11.0.667_i86linux32)

cd /home/project/11003581/Tools/
wget https://www.drive5.com/downloads/usearch11.0.667_i86linux32.gz
gunzip usearch11.0.667_i86linux32.gz
mv usearch11.0.667_i86linux32 usearch
nano ~/.bashrc

Then add this line in the end of the file: export PATH=$PATH:/home/project/11003581/Tools

source ~/.bashrc



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

cd /home/project/11003581/Data/JQQ/01.RawData/Rep1/

cutadapt --no-indels -q 30,30 --trimmed-only -j 32 -m 10 \
    -a file:/home/project/11003581/Data/JQQ-analysis/approach2/r1/barcodes.fa \
    -o /home/project/11003581/Data/JQQ-analysis/approach2/r1/{name}.1.fq.gz \
    -p /home/project/11003581/Data/JQQ-analysis/approach2/r1/{name}.2.fq.gz \
    Rep1_DKDL240034859-1A_22TKKWLT3_L3_1.fq.gz \
    Rep1_DKDL240034859-1A_22TKKWLT3_L3_2.fq.gz \
    1> /home/project/11003581/Data/JQQ-analysis/approach2/r1/log.txt



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




##############
############## adding prefix to a file
for file in *.gz; do mv "$file" "rep2_$file"; done



##############
################# making sample.txt with conditions
##############
for file in */*.gz; do
    # Extract the filename
    filename=$(basename "$file")
    
    # Determine the condition
    if [[ "$filename" == *"WT"* ]]; then
        condition="WT"
    elif [[ "$filename" == *"DOX5d"* ]]; then
        condition="DOX5d"
    elif [[ "$filename" == *"DOX20d"* ]]; then
        condition="DOX20d"
    else
        condition="Unknown"
    fi
    
    # Write to the output file
    echo -e "$file\t$condition" >> sample.txt
done



######## demultiplexing with entire adapter


#!/bin/bash

#PBS -l select=1:ngpus=1:ncpus=64
#PBS -l walltime=15:00:00
#PBS -P 11003581
#PBS -N demultiplex
#PBS -j oe

module load miniforge3
conda activate /home/project/11003581/conda-envs/mimseq

for i in /home/project/11003581/Data/JQQ/01.RawData/Rep1

do fn="tRNA1"

cutadpat --no-indels -q 30,30 --trimmed-only -j 10 -a file:barcodes.fa -m 10 -o $fn'_{name}_trim.fq.gz' $i 1> $fn'_log.txt'

done
