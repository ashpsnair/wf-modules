#!/bin/bash
#PBS -l select=1:ncpus=64:ngpus=1
#PBS -l walltime=45:00:00
#PBS -P 11003581
#PBS -N trial-nf-nanoseq
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR


module load java/17.0.6-jdk
module load singularity

/home/project/11003581/Tools/nextflow run nf-core/nanoseq \
    --input /home/project/11003581/Data/Ash/nf-core-p3l/samplesheet.csv \
    --protocol DNA \
    --call_variants \
    --skip_demultiplexing \
    --skip_quantification \
    -profile singularity



'''

group,replicate,barcode,input_file,fasta,gtf
WT_MCF10A,1,,/home/project/11003581/Data/Ash/P3L-lab-analysis/sorted_bams/WT_MCF10A_hac_sorted.bam,/home/project/11003581/Ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna,
P3L,1,,/home/project/11003581/Data/Ash/P3L-lab-analysis/sorted_bams/p3292l_hac_sorted.bam,/home/project/11003581/Ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna,


'''