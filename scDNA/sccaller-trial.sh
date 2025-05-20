#!/bin/bash

#PBS -l select=3:ncpus=64:mem=256g
#PBS -l walltime=24:00:00
#PBS -P 11003581
#PBS -N run-somatic-scDNA
#PBS -j oe

source /home/project/11003581/Tools/miniforge3/bin/activate
conda activate /home/project/11003581/conda-envs/sccaller

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

python /home/users/nus/ash.ps/scratch/scDNA/analysis/logs/sccaller_v2.0.0.py \
    --bam /home/users/nus/ash.ps/scratch/scDNA/DNA/Secondary-Analysis-DNA/2204-A2/2204-A2.sorted.bqsr.dedup.bam \
    --bulk /home/users/nus/ash.ps/scratch/scDNA/Bulk-data/Secondary-Analysis/2204/2204.sorted.bqsr.dedup.bam \
    --fasta /home/project/11003581/Ref/Homo_sapiens/Ensembl/Homo_sapiens.GRCh38.dna_rm.primary_assembly.fa \
    --output /home/users/nus/ash.ps/scratch/scDNA/analysis/somatic/somatic-2204-A2.vcf \
    --snp_type hsnp \
    --snp_in /home/users/nus/ash.ps/scratch/scDNA/analysis/somatic/2204.g.vcf \
    --cpu_num 128 \
    --engine samtools