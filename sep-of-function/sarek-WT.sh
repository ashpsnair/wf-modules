#samplesheet

#b2
patient,sex,status,sample,lane,fastq_1,fastq_2
WT,XX,1,WT_P30_4,2,/home/project/11003581/Data/sep-function/b1-fastqs/WT/WT_P30_4/WT_P30_1_EDHG220013938-3a_HKFHKDSX3_L2_1.fq.gz,/home/project/11003581/Data/sep-function/b1-fastqs/WT/WT_P30_4/WT_P30_1_EDHG220013938-3a_HKFHKDSX3_L2_2.fq.gz
WT,XX,1,WT_P30_4,4,/home/project/11003581/Data/sep-function/b1-fastqs/WT/WT_P30_4/WT_P30_1_EDHG220013938-3a_HKFHKDSX3_L4_1.fq.gz,/home/project/11003581/Data/sep-function/b1-fastqs/WT/WT_P30_4/WT_P30_1_EDHG220013938-3a_HKFHKDSX3_L4_2.fq.gz
WT,XX,1,WT_P30_5,1,/home/project/11003581/Data/sep-function/b1-fastqs/WT/WT_P30_5/WT_P30_2_EDHG220013939-3a_HKFNWDSX3_L1_1.fq.gz,/home/project/11003581/Data/sep-function/b1-fastqs/WT/WT_P30_5/WT_P30_2_EDHG220013939-3a_HKFNWDSX3_L1_2.fq.gz
WT,XX,1,WT_P30_6,1,/home/project/11003581/Data/sep-function/b1-fastqs/WT/WT_P30_6/WT_P30_3_EDHG220013940-3a_HKFNWDSX3_L1_1.fq.gz,/home/project/11003581/Data/sep-function/b1-fastqs/WT/WT_P30_6/WT_P30_3_EDHG220013940-3a_HKFNWDSX3_L1_2.fq.gz
WT,XX,1,WT_P30_7,1,/home/project/11003581/Data/sep-function/b1-fastqs/WT/WT_P30_7/WT_P30_4_EDHG220013941-3a_HKFNWDSX3_L1_1.fq.gz,/home/project/11003581/Data/sep-function/b1-fastqs/WT/WT_P30_7/WT_P30_4_EDHG220013941-3a_HKFNWDSX3_L1_2.fq.gz
WT,XX,1,WT_P30_8,1,/home/project/11003581/Data/sep-function/b1-fastqs/WT/WT_P30_8/WT_P30_5_EDHG220013942-3a_HKFNWDSX3_L1_1.fq.gz,/home/project/11003581/Data/sep-function/b1-fastqs/WT/WT_P30_8/WT_P30_5_EDHG220013942-3a_HKFNWDSX3_L1_2.fq.gz
WT,XX,1,WT_P30_9,4,/home/project/11003581/Data/sep-function/b1-fastqs/WT/WT_P30_9/WT_P30_6_EDHG220013943-3a_HKFMGDSX3_L4_1.fq.gz,/home/project/11003581/Data/sep-function/b1-fastqs/WT/WT_P30_9/WT_P30_6_EDHG220013943-3a_HKFMGDSX3_L4_2.fq.gz
WT,XX,1,WT_P30_9,1,/home/project/11003581/Data/sep-function/b1-fastqs/WT/WT_P30_9/WT_P30_6_EDHG220013943-3a_HKFNWDSX3_L1_1.fq.gz,/home/project/11003581/Data/sep-function/b1-fastqs/WT/WT_P30_9/WT_P30_6_EDHG220013943-3a_HKFNWDSX3_L1_2.fq.gz
WT,XX,1,WT_P30_10,1,/home/project/11003581/Data/sep-function/b1-fastqs/WT/WT_P30_10/WT_P30_7_EDHG220013944-3a_HKFNWDSX3_L1_1.fq.gz,/home/project/11003581/Data/sep-function/b1-fastqs/WT/WT_P30_10/WT_P30_7_EDHG220013944-3a_HKFNWDSX3_L1_2.fq.gz


#b4
patient,sex,status,sample,lane,fastq_1,fastq_2
WT,XX,1,WT_P30_1,1,/home/project/11003581/Data/sep-function/b1-fastqs/WT/WT_P30_1/WT_P30_1_EDHG210027025-1a_HVCK3DSX2_L4_1.fq.gz,/home/project/11003581/Data/sep-function/b1-fastqs/WT/WT_P30_1/WT_P30_1_EDHG210027025-1a_HVCK3DSX2_L4_2.fq.gz
WT,XX,1,WT_P30_2,1,
WT,XX,1,WT_P30_3,1,
WT,XX,1,WT_P30_B,1,

#!/bin/bash

#PBS -l select=1:ncpus=128:mem=256g
#PBS -l walltime=24:00:00
#PBS -P 11003581
#PBS -N run-WT
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load nextflow/24.10.5
module load miniforge3

export NXF_CONDA_CACHEDIR=/home/project/11003581/Tools/nf-conda-cache
export NXF_CONDA_USE_MAMBA=true

WORKDIR=/home/users/nus/ash.ps/scratch/sep-fun-analysis/WT/b1

nextflow run nf-core/sarek -r 3.5.1 \
   -profile conda \
   --input $WORKDIR/samplesheet.csv \
   --outdir $WORKDIR
