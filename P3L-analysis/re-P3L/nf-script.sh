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
    --sample_name 'P3L-shallow' \
    -profile singularity \
    --basecaller_cfg 'dna_r10.4.1_e8.2_400bps_hac@v5.0.0' \
    --dorado_ext 'fast5' \
    --input '/home/project/11003581/Data/HC/ONT/P3292L/TJPROJ6/TGS/haiwai/haiwai/HW_ONT_qc/X401SC23084120-Z01-F001/data_release/X401SC23084120-Z01-F001/raw_data/P3292L/20231210_1504_6G_PAS94472_b08ff6e9/fast5_pass/' \
    --ref '/home/project/11003581/Ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna' \
    --out_dir '/home/users/nus/ash.ps/scratch/P3L-nf/shallow/basecalling/' \
    --remora_cfg 'dna_r10.4.1_e8.2_400bps_hac@v4.1.0_5mCG_5hmCG@v2' \
    --output_fmt 'bam' \
    --basecaller_basemod_threads 64 \
    --ubam_map_threads 32 \
    --ubam_sort_threads 16 \
    --ubam_bam2fq_threads 16 \
    --merge_threads 64 \
    --stats_threads 64


'''

## ERROR 
FATAL:   While making image from oci registry: error fetching image to cache: while building SIF from layers: conveyor failed to get: no descriptor found for reference "sha256.d0d0abdeaff2ccc5f2dd3a25ba8ceb5cea12edb2e8e7901b89f29024b063290f"

#RESOLVED by the following 
ls -l $NXF_SINGULARITY_CACHEDIR
ls -l $SINGULARITY_CACHEDIR

Clear temporary files:
bash
rm -rf $SINGULARITY_CACHEDIR/oci-tmp/*
Manually download the image:
bash
cd $NXF_SINGULARITY_CACHEDIR
singularity pull docker://quay.io/biocontainers/samtools:1.15.1--h1170115_0

'''


#!/bin/bash

#PBS -l select=1:ncpus=128
#PBS -l walltime=18:00:00
#PBS -P 11003581
#PBS -N epi2me-sv-cnv-shallow
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load java/17.0.6-jdk
module load singularity

export NXF_SINGULARITY_CACHEDIR=/home/project/11003581/Tools/singularity-cache/
export SINGULARITY_CACHEDIR=/home/project/11003581/Tools/singularity-cache/


/home/project/11003581/Tools/nextflow run epi2me-labs/wf-human-variation \
    --sample_name 'P3L-shallow' \
    --bam '/home/users/nus/ash.ps/scratch/P3L-nf/shallow/basecalling/SAMPLE.pass.cram' \
    --out_dir '/home/users/nus/ash.ps/scratch/P3L-nf/shallow/variant_calling' \
    --output_gene_summary 1 \
    --sv \
    --cnv \
    --ref '/home/project/11003581/Ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna' \
    --use_qdnaseq \
    --threads 128 \
    --ubam_map_threads 128 \
    --ubam_sort_threads 128 \
    --ubam_bam2fq_threads 128 \
    --profile singularity


################################# Analysing deep seq #####################
# Define variables for the tar file and destination directory
fast5_dir="/home/users/nus/ash.ps/scratch/P3L-nf/deep/fast5/"
fast5_tar="/home/project/11003581/Data/HC/LongReadSeq_B2/X401SC23084120-Z01-F002_Data_Release_20240322/Data-X401SC23084120-Z01-F002/MCF10A_P3L_2/1101_1E_PAS74408_f5116582/MCF10A_P3L_2_fast5_pass.tar"

# Create the destination directory if it doesn't exist
mkdir -p "$fast5_dir"

# Extract the tar file into the destination directory
tar -xvf "$fast5_tar" -C "$fast5_dir"



########## Running basecalling

#!/bin/bash

#PBS -l select=4:ngpus=2
#PBS -l walltime=5:00:00
#PBS -P 11003581
#PBS -N epi2me-basecalling-deep
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR


module load java/17.0.6-jdk
module load singularity

export NXF_SINGULARITY_CACHEDIR=/home/project/11003581/Tools/singularity-cache/
export SINGULARITY_CACHEDIR=/home/project/11003581/Tools/singularity-cache/

/home/project/11003581/Tools/nextflow run epi2me-labs/wf-basecalling \
    --sample_name 'P3L-shallow' \
    -profile singularity \
    --basecaller_cfg 'dna_r10.4.1_e8.2_400bps_hac@v5.0.0' \
    --dorado_ext 'fast5' \
    --input '/home/users/nus/ash.ps/scratch/P3L-nf/deep/fast5/MCF10A_P3L_2_fast5_pass/' \
    --ref '/home/project/11003581/Ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna' \
    --out_dir '/home/users/nus/ash.ps/scratch/P3L-nf/deep/basecalling/' \
    --remora_cfg 'dna_r10.4.1_e8.2_400bps_hac@v4.1.0_5mCG_5hmCG@v2' \
    --output_fmt 'bam' \
    --basecaller_basemod_threads 64 \
    --ubam_map_threads 32 \
    --ubam_sort_threads 16 \
    --ubam_bam2fq_threads 16 \
    --merge_threads 64 \
    --stats_threads 64