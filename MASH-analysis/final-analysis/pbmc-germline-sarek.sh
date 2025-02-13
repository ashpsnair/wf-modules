#!/bin/bash

#PBS -l select=3:ncpus=64:mem=256g
#PBS -l walltime=24:00:00
#PBS -P 11003581
#PBS -N mash-pbmc-germline-sarek
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load java/17.0.6-jdk
module load singularity

export NXF_SINGULARITY_CACHEDIR=/home/project/11003581/Tools/singularity-cache/
export SINGULARITY_CACHEDIR=/home/project/11003581/Tools/singularity-cache/

/home/project/11003581/Tools/nextflow run nf-core/sarek -r 3.4.4 \
   -name pbmc-germline-call \
   -profile singularity \
   --input /home/users/nus/ash.ps/scratch/NCCS-MASH/analysis/PBMC/samplesheet.csv \
   --outdir /home/users/nus/ash.ps/scratch/NCCS-MASH/analysis/PBMC/ \
   --genome GATK.GRCh38 \
   --tools freebayes,deepvariant,strelka,haplotypecaller,manta,tiddit,cnvkit \
   --max_cpus 128 \
   --max_memory '256.GB' \
   --email ash.ps@nus.edu.sg
   