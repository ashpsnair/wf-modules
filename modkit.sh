#!/bin/bash

#PBS -l select=1:ngpus=1
#PBS -l walltime=10:00:00
#PBS -P personal-ash.ps
#PBS -N trial-clairs
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

output_dir="/home/users/nus/ash.ps/scratch/"
ref_fasta="/home/project/11003581/Ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

samplename="MCF10A_Dox_D3"
input_bam="/home/svu/ash.ps/MCF10A_Dox_D3/MCF10A_Dox_D3_calls_m6A_pseU_Q10.bam"

modkit pileup $INPUT_BAM $OUTPUT_DIR/{$SAMPLENAME}_pileup.bed \
  --ref $REFERENCE \
  --preset traditional

modkit sample-probs input.bam 

modkit summary $INPUT_BAM --no-sampling 