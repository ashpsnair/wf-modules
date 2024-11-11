bwa_path= /home/project/11003581/Tools/bwa/bwa 

reference= /home/project/11003581/Ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna 

patient,sex,status,sample,lane,fastq_1,fastq_2
1,XY,0,N01,L1,/home/project/11003581/Data/YS1/01.RawData/N01/N01_DKDN240040899-1A_HVNMJDSXC_L1_1.fq.gz,/home/project/11003581/Data/YS1/01.RawData/N01/N01_DKDN240040899-1A_HVNMJDSXC_L1_2.fq.gz
1,XY,0,N01,L2,/home/project/11003581/Data/YS1/01.RawData/N01/N01_DKDN240040899-1A_HVNMJDSXC_L2_1.fq.gz,/home/project/11003581/Data/YS1/01.RawData/N01/N01_DKDN240040899-1A_HVNMJDSXC_L2_2.fq.gz
1,XY,0,N01,L3,/home/project/11003581/Data/YS1/01.RawData/N01/N01_DKDN240040899-1A_HVNMJDSXC_L3_1.fq.gz,/home/project/11003581/Data/YS1/01.RawData/N01/N01_DKDN240040899-1A_HVNMJDSXC_L3_2.fq.gz
1,XY,0,N01,L4,/home/project/11003581/Data/YS1/01.RawData/N01/N01_DKDN240040899-1A_HVNMJDSXC_L4_1.fq.gz,/home/project/11003581/Data/YS1/01.RawData/N01/N01_DKDN240040899-1A_HVNMJDSXC_L4_2.fq.gz
1,XY,1,T01,L1,/home/project/11003581/Data/YS1/01.RawData/T01/T01_DKDN240040885-1A_HVNMJDSXC_L1_1.fq.gz,/home/project/11003581/Data/YS1/01.RawData/T01/T01_DKDN240040885-1A_HVNMJDSXC_L1_2.fq.gz
1,XY,1,T01,L2,/home/project/11003581/Data/YS1/01.RawData/T01/T01_DKDN240040885-1A_HVNMJDSXC_L2_1.fq.gz,/home/project/11003581/Data/YS1/01.RawData/T01/T01_DKDN240040885-1A_HVNMJDSXC_L2_2.fq.gz
1,XY,1,T01,L3,/home/project/11003581/Data/YS1/01.RawData/T01/T01_DKDN240040885-1A_HVNMJDSXC_L3_1.fq.gz,/home/project/11003581/Data/YS1/01.RawData/T01/T01_DKDN240040885-1A_HVNMJDSXC_L3_2.fq.gz
1,XY,1,T01,L4,/home/project/11003581/Data/YS1/01.RawData/T01/T01_DKDN240040885-1A_HVNMJDSXC_L4_1.fq.gz,/home/project/11003581/Data/YS1/01.RawData/T01/T01_DKDN240040885-1A_HVNMJDSXC_L4_2.fq.gz



#!/bin/bash

#PBS -l select=1:ncpus=64:mem=256g
#PBS -l walltime=15:00:00
#PBS -P 11003581
#PBS -N trial-gatk
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

# Define paths
bwa_path="/home/project/11003581/Tools/bwa/bwa"
reference="/home/project/11003581/Ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
output_dir="/home/project/11003581/Data/YS-gatk"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Define input CSV file
input_csv="/home/project/11003581/Data/YS-gatk/samples.csv"  # Ensure this file is in the working directory

# Read CSV and process each sample
while IFS=',' read -r patient sex status sample lane fastq_1 fastq_2; do
    # Skip the header line
    if [[ "$patient" == "patient" ]]; then
        continue
    fi

    # Determine if the sample is tumor or normal based on status
    if [[ "$status" == "1" ]]; then
        sample_type="tumor"
    else
        sample_type="normal"
    fi

    # Concatenate FastQ files for each lane into merged files for tumor and normal samples
    merged_fastq_1="$output_dir/${sample}_R1.fastq"
    merged_fastq_2="$output_dir/${sample}_R2.fastq"

    # Concatenate FastQ files for the current lane
    echo "Concatenating FastQ files for sample $sample ($sample_type), lane $lane..."
    cat "$fastq_1" >> "$merged_fastq_1"
    cat "$fastq_2" >> "$merged_fastq_2"

done < "$input_csv"


######### Running BWA MEM
/home/project/11003581/Tools/bwa/bwa mem /home/project/11003581/Ref/hg38.fa.gz N01_R1.fastq  N01_R2.fastq > N01.sam