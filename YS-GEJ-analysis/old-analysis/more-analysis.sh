#!/bin/bash

#PBS -l select=2:ncpus=64:mem=256g
#PBS -l walltime=12:00:00
#PBS -P 11003581
#PBS -N YS-extra-b1
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR


module load java/17.0.6-jdk
module load singularity

/home/project/11003581/Tools/nextflow run nf-core/sarek -r 3.4.4 \
   -profile singularity \
   --step variant_calling \
   --input /home/users/nus/ash.ps/scratch/YS-rerun-nf/rerun-extra-analysis/b1/samplesheet.csv \
   --outdir /home/users/nus/ash.ps/scratch/YS-rerun-nf/rerun-extra-analysis/b1 \
   --tools manta,ascat,msisensorpro \
   --max_cpus 128 \
   --max_memory '256.GB' \



###CRAM files
#batch1
patient,sex,status,sample,cram,crai
patient01,XY,0,N01,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS1-4//preprocessing/recalibrated/N01/N01.recal.cram,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS1-4//preprocessing/recalibrated/N01/N01.recal.cram.crai
patient01,XY,1,T01,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS1-4//preprocessing/recalibrated/T01/T01.recal.cram,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS1-4//preprocessing/recalibrated/T01/T01.recal.cram.crai
patient02,XY,0,N02,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS1-4//preprocessing/recalibrated/N02/N02.recal.cram,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS1-4//preprocessing/recalibrated/N02/N02.recal.cram.crai
patient02,XY,1,T02,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS1-4//preprocessing/recalibrated/T02/T02.recal.cram,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS1-4//preprocessing/recalibrated/T02/T02.recal.cram.crai
patient03,XX,0,N03,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS1-4//preprocessing/recalibrated/N03/N03.recal.cram,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS1-4//preprocessing/recalibrated/N03/N03.recal.cram.crai
patient03,XX,1,T03,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS1-4//preprocessing/recalibrated/T03/T03.recal.cram,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS1-4//preprocessing/recalibrated/T03/T03.recal.cram.crai
patient04,XY,0,N04,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS1-4//preprocessing/recalibrated/N04/N04.recal.cram,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS1-4//preprocessing/recalibrated/N04/N04.recal.cram.crai
patient04,XY,1,T04,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS1-4//preprocessing/recalibrated/T04/T04.recal.cram,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS1-4//preprocessing/recalibrated/T04/T04.recal.cram.crai
patient05,XY,0,N05,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS5-8//preprocessing/recalibrated/N05/N05.recal.cram,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS5-8//preprocessing/recalibrated/N05/N05.recal.cram.crai
patient05,XY,1,T05,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS5-8//preprocessing/recalibrated/T05/T05.recal.cram,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS5-8//preprocessing/recalibrated/T05/T05.recal.cram.crai


#batc2
patient,sex,status,sample,cram,crai
patient06,XX,0,N06,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS5-8//preprocessing/recalibrated/N06/N06.recal.cram,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS5-8//preprocessing/recalibrated/N06/N06.recal.cram.crai
patient06,XX,1,T06,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS5-8//preprocessing/recalibrated/T06/T06.recal.cram,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS5-8//preprocessing/recalibrated/T06/T06.recal.cram.crai
patient07,XY,0,N07,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS5-8//preprocessing/recalibrated/N07/N07.recal.cram,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS5-8//preprocessing/recalibrated/N07/N07.recal.cram.crai
patient07,XY,1,T07,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS5-8//preprocessing/recalibrated/T07/T07.recal.cram,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS5-8//preprocessing/recalibrated/T07/T07.recal.cram.crai
patient08,XY,0,N08,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS5-8//preprocessing/recalibrated/N08/N08.recal.cram,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS5-8//preprocessing/recalibrated/N08/N08.recal.cram.crai
patient08,XY,1,T08,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS5-8//preprocessing/recalibrated/T08/T08.recal.cram,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS5-8//preprocessing/recalibrated/T08/T08.recal.cram.crai
patient09,XY,0,N09,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS9-11//preprocessing/recalibrated/N09/N09.recal.cram,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS9-11//preprocessing/recalibrated/N09/N09.recal.cram.crai
patient09,XY,1,T09,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS9-11//preprocessing/recalibrated/T09/T09.recal.cram,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS9-11//preprocessing/recalibrated/T09/T09.recal.cram.crai
patient10,XY,0,N10,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS9-11//preprocessing/recalibrated/N10/N10.recal.cram,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS9-11//preprocessing/recalibrated/N10/N10.recal.cram.crai
patient10,XY,1,T10,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS9-11//preprocessing/recalibrated/T10/T10.recal.cram,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS9-11//preprocessing/recalibrated/T10/T10.recal.cram.crai


#Batch3
patient,sex,status,sample,cram,crai
patient11,XY,0,N11,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS9-11//preprocessing/recalibrated/N11/N11.recal.cram,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS9-11//preprocessing/recalibrated/N11/N11.recal.cram.crai
patient11,XY,1,T11,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS9-11//preprocessing/recalibrated/T11/T11.recal.cram,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS9-11//preprocessing/recalibrated/T11/T11.recal.cram.crai
patient12,XY,0,N12,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS12-14//preprocessing/recalibrated/N12/N12.recal.cram,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS12-14//preprocessing/recalibrated/N12/N12.recal.cram.crai
patient12,XY,1,T12,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS12-14//preprocessing/recalibrated/T12/T12.recal.cram,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS12-14//preprocessing/recalibrated/T12/T12.recal.cram.crai
patient13,XY,0,N13,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS12-14//preprocessing/recalibrated/N13/N13.recal.cram,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS12-14//preprocessing/recalibrated/N13/N13.recal.cram.crai
patient13,XY,1,T13,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS12-14//preprocessing/recalibrated/T13/T13.recal.cram,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS12-14//preprocessing/recalibrated/T13/T13.recal.cram.crai
patient14,XY,0,N14,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS12-14//preprocessing/recalibrated/N14/N14.recal.cram,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS12-14//preprocessing/recalibrated/N14/N14.recal.cram.crai
patient14,XY,1,T14,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS12-14//preprocessing/recalibrated/T14/T14.recal.cram,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS12-14//preprocessing/recalibrated/T14/T14.recal.cram.crai
