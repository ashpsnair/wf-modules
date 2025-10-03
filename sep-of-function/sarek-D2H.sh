#D2H samplesheet

#b2
patient,sex,status,sample,lane,fastq_1,fastq_2
D2H,XX,1,D2H_P30_4,1,/home/project/11003581/Data/sep-function/b1-fastqs/D2H/D2H_P30_4/D2H_P30_1_EDHG220013966-2a_HL5Y2DSX3_L4_1.fq.gz,/home/project/11003581/Data/sep-function/b1-fastqs/D2H/D2H_P30_4/D2H_P30_1_EDHG220013966-2a_HL5Y2DSX3_L4_2.fq.gz
D2H,XX,1,D2H_P30_5,1,/home/project/11003581/Data/sep-function/b1-fastqs/D2H/D2H_P30_5/D2H_P30_2_EDHG220013967-2a_HL5Y2DSX3_L4_1.fq.gz,/home/project/11003581/Data/sep-function/b1-fastqs/D2H/D2H_P30_5/D2H_P30_2_EDHG220013967-2a_HL5Y2DSX3_L4_2.fq.gz
D2H,XX,1,D2H_P30_6,1,/home/project/11003581/Data/sep-function/b1-fastqs/D2H/D2H_P30_6/D2H_P30_3_EDHG220013968-2a_HL5Y2DSX3_L4_1.fq.gz,/home/project/11003581/Data/sep-function/b1-fastqs/D2H/D2H_P30_6/D2H_P30_3_EDHG220013968-2a_HL5Y2DSX3_L4_2.fq.gz
D2H,XX,1,D2H_P30_7,1,/home/project/11003581/Data/sep-function/b1-fastqs/D2H/D2H_P30_7/D2H_P30_4_EDHG220013969-2a_HL5Y2DSX3_L4_1.fq.gz,/home/project/11003581/Data/sep-function/b1-fastqs/D2H/D2H_P30_7/D2H_P30_4_EDHG220013969-2a_HL5Y2DSX3_L4_2.fq.gz
D2H,XX,1,D2H_P30_8,1,/home/project/11003581/Data/sep-function/b1-fastqs/D2H/D2H_P30_8/D2H_P30_5_EDHG220013970-2a_HL5Y2DSX3_L4_1.fq.gz,/home/project/11003581/Data/sep-function/b1-fastqs/D2H/D2H_P30_8/D2H_P30_5_EDHG220013970-2a_HL5Y2DSX3_L4_2.fq.gz
D2H,XX,1,D2H_P30_9,1,/home/project/11003581/Data/sep-function/b1-fastqs/D2H/D2H_P30_9/D2H_P30_6_EDHG220013971-2a_HL5Y2DSX3_L4_1.fq.gz,/home/project/11003581/Data/sep-function/b1-fastqs/D2H/D2H_P30_9/D2H_P30_6_EDHG220013971-2a_HL5Y2DSX3_L4_2.fq.gz
D2H,XX,1,D2H_P30_10,1,/home/project/11003581/Data/sep-function/b1-fastqs/D2H/D2H_P30_10/D2H_P30_7_EDHG220013972-2a_HL5Y2DSX3_L1_1.fq.gz,/home/project/11003581/Data/sep-function/b1-fastqs/D2H/D2H_P30_10/D2H_P30_7_EDHG220013972-2a_HL5Y2DSX3_L1_2.fq.gz


#b4
patient,sex,status,sample,lane,fastq_1,fastq_2
D2H,XX,1,D2H_P30_1,1,
D2H,XX,1,D2H_P30_2,1,
D2H,XX,1,D2H_P30_3,1,
D2H,XX,1,D2H_P30_B,1,

##script
/home/users/nus/ash.ps/scratch/sep-fun-analysis/D2H/b1

#!/bin/bash

#PBS -l select=1:ncpus=128:mem=256g
#PBS -l walltime=24:00:00
#PBS -P 11003581
#PBS -N run-D2H
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load nextflow/24.10.5
module load miniforge3

export NXF_CONDA_CACHEDIR=/home/project/11003581/Tools/nf-conda-cache
export NXF_CONDA_USE_MAMBA=true

WORKDIR=/home/users/nus/ash.ps/scratch/sep-fun-analysis/D2H/b1

nextflow run nf-core/sarek -r 3.5.1 \
   -profile conda \
   --input $WORKDIR/samplesheet.csv \
   --outdir $WORKDIR \