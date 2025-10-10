#samplesheet
patient,sex,status,sample,cram,crai
BRCA1_R133C,XX,0,WT_P30_B,/home/users/nus/ash.ps/scratch/sep-fun-analysis/CONTROL/preprocessing/recalibrated/WT_P30_B/WT_P30_B.recal.cram,/home/users/nus/ash.ps/scratch/sep-fun-analysis/CONTROL/preprocessing/recalibrated/WT_P30_B/WT_P30_B.recal.cram.crai
BRCA1_R133C,XX,1,B_68_01,/home/users/nus/ash.ps/scratch/sep-fun-analysis/BRCA1_R133C/preprocessing/recalibrated/B_68_01/B_68_01.recal.cram,/home/users/nus/ash.ps/scratch/sep-fun-analysis/BRCA1_R133C/preprocessing/recalibrated/B_68_01/B_68_01.recal.cram.crai
BRCA1_R133C,XX,1,B_68_02,/home/users/nus/ash.ps/scratch/sep-fun-analysis/BRCA1_R133C/preprocessing/recalibrated/B_68_02/B_68_02.recal.cram,/home/users/nus/ash.ps/scratch/sep-fun-analysis/BRCA1_R133C/preprocessing/recalibrated/B_68_02/B_68_02.recal.cram.crai
BRCA1_R133C,XX,1,B_68_03,/home/users/nus/ash.ps/scratch/sep-fun-analysis/BRCA1_R133C/preprocessing/recalibrated/B_68_03/B_68_03.recal.cram,/home/users/nus/ash.ps/scratch/sep-fun-analysis/BRCA1_R133C/preprocessing/recalibrated/B_68_03/B_68_03.recal.cram.crai
BRCA1_R133C,XX,1,B_68_04,/home/users/nus/ash.ps/scratch/sep-fun-analysis/BRCA1_R133C/preprocessing/recalibrated/B_68_04/B_68_04.recal.cram,/home/users/nus/ash.ps/scratch/sep-fun-analysis/BRCA1_R133C/preprocessing/recalibrated/B_68_04/B_68_04.recal.cram.crai
BRCA1_R133C,XX,1,B_68_05,/home/users/nus/ash.ps/scratch/sep-fun-analysis/BRCA1_R133C/preprocessing/recalibrated/B_68_05/B_68_05.recal.cram,/home/users/nus/ash.ps/scratch/sep-fun-analysis/BRCA1_R133C/preprocessing/recalibrated/B_68_05/B_68_05.recal.cram.crai
BRCA1_R133C,XX,1,B_68_06,/home/users/nus/ash.ps/scratch/sep-fun-analysis/BRCA1_R133C/preprocessing/recalibrated/B_68_06/B_68_06.recal.cram,/home/users/nus/ash.ps/scratch/sep-fun-analysis/BRCA1_R133C/preprocessing/recalibrated/B_68_06/B_68_06.recal.cram.crai
BRCA1_R133C,XX,1,B_68_07,/home/users/nus/ash.ps/scratch/sep-fun-analysis/BRCA1_R133C/preprocessing/recalibrated/B_68_07/B_68_07.recal.cram,/home/users/nus/ash.ps/scratch/sep-fun-analysis/BRCA1_R133C/preprocessing/recalibrated/B_68_07/B_68_07.recal.cram.crai
BRCA1_R133C,XX,1,B_68_08,/home/users/nus/ash.ps/scratch/sep-fun-analysis/BRCA1_R133C/preprocessing/recalibrated/B_68_08/B_68_08.recal.cram,/home/users/nus/ash.ps/scratch/sep-fun-analysis/BRCA1_R133C/preprocessing/recalibrated/B_68_08/B_68_08.recal.cram.crai
BRCA1_R133C,XX,1,B_68_09,/home/users/nus/ash.ps/scratch/sep-fun-analysis/BRCA1_R133C/preprocessing/recalibrated/B_68_09/B_68_09.recal.cram,/home/users/nus/ash.ps/scratch/sep-fun-analysis/BRCA1_R133C/preprocessing/recalibrated/B_68_09/B_68_09.recal.cram.crai
BRCA1_R133C,XX,1,B_68_10,/home/users/nus/ash.ps/scratch/sep-fun-analysis/BRCA1_R133C/preprocessing/recalibrated/B_68_10/B_68_10.recal.cram,/home/users/nus/ash.ps/scratch/sep-fun-analysis/BRCA1_R133C/preprocessing/recalibrated/B_68_10/B_68_10.recal.cram.crai


#!/bin/bash

#PBS -l select=3:ncpus=64:mem=128g
#PBS -l walltime=24:00:00
#PBS -P 11003581
#PBS -N B1_R133C-vc
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load nextflow/24.10.5
module load miniforge3

export NXF_CONDA_CACHEDIR=/home/project/11003581/Tools/nf-conda-cache
export NXF_CONDA_USE_MAMBA=true

WORKDIR=/home/users/nus/ash.ps/scratch/sep-fun-analysis/variant_calling/B1_R133C/

nextflow run nf-core/sarek -r 3.5.1 \
   -profile mamba \
   --input $WORKDIR/samplesheet.csv \
   --outdir $WORKDIR \
   --step variant_calling \
   --tools mutect2,strelka,ascat,manta,snpeff,vep,msisensorpro \
   --joint_mutect2 \
   --pon /home/project/11003581/Ref/pons/somatic-hg38_1000g_pon.hg38.vcf.gz \
   --pon_tbi /home/project/11003581/Ref/pons/somatic-hg38_1000g_pon.hg38.vcf.gz.tbi



