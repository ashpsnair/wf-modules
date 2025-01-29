#!/bin/bash

#PBS -l select=4:ngpus=2
#PBS -l walltime=18:00:00
#PBS -P 11003581
#PBS -N epi2me-basecalling-shallow
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR


module load java/17.0.6-jdk
module load singularity

export NXF_SINGULARITY_CACHEDIR=/home/project/11003581/Tools/singularity-cache/
export SINGULARITY_CACHEDIR=/home/project/11003581/Tools/singularity-cache/

/home/project/11003581/Tools/nextflow run epi2me-labs/wf-basecalling \
    -profile singularity \
    --basecaller_cfg 'dna_r10.4.1_e8.2_400bps_hac@v5.0.0' \
    --dorado_ext 'fast5' \
    --input '/home/project/11003581/Data/HC/ONT/P3292L/TJPROJ6/TGS/haiwai/haiwai/HW_ONT_qc/X401SC23084120-Z01-F001/data_release/X401SC23084120-Z01-F001/raw_data/P3292L/20231210_1504_6G_PAS94472_b08ff6e9/fast5_pass/' \
    --ref '/home/project/11003581/Ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna' \
    --out_dir '/home/users/nus/ash.ps/scratch/P3L-nf/shallow/' \
    --remora_cfg 'dna_r10.4.1_e8.2_400bps_hac@v4.1.0_5mCG_5hmCG@v2' \
    --basecaller_basemod_threads 128 \
    --ubam_map_threads 64 \
    --ubam_sort_threads 32 \
    --ubam_bam2fq_threads 32 \
    --merge_threads 64 \
    --stats_threads 64


export NXF_SINGULARITY_CACHEDIR=/home/project/11003581/Tools/singularity-cache/
export SINGULARITY_CACHEDIR=/home/project/11003581/Tools/singularity-cache/


/home/project/11003581/Tools/nextflow run epi2me-labs/wf-human-variation \
    --sample_name 'P3L-shallow' \
    --bam 'wf-human-variation-demo/demo.bam' \
    --ref 'wf-human-variation-demo/demo.fasta' \
    --out_dir \
    --output_gene_summary \
    --sv \
    --cnv \
    --ref 'wf-basecalling-demo/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta' \
    --use_qdnaseq \
    --threads 128 \
    --ubam_map_threads 128 \
    --ubam_sort_threads 128 \
    --ubam_bam2fq_threads 128


