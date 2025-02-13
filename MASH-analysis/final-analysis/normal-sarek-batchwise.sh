############ SAREK nf- code

#!/bin/bash

#PBS -l select=3:ncpus=64:mem=256g
#PBS -l walltime=24:00:00
#PBS -P 11003581
#PBS -N sarek-mash-b1
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load java/17.0.6-jdk
module load singularity

export NXF_SINGULARITY_CACHEDIR=/home/project/11003581/Tools/singularity-cache/
export SINGULARITY_CACHEDIR=/home/project/11003581/Tools/singularity-cache/

/home/project/11003581/Tools/nextflow run nf-core/sarek -r 3.4.4 \
   -profile singularity \
   --input /home/users/nus/ash.ps/scratch/NCCS-MASH/analysis/b1/samplesheet.csv \
   --outdir /home/users/nus/ash.ps/scratch/NCCS-MASH/analysis/b1/ \
   --genome GATK.GRCh38 \
   --tools mutect2,freebayes,ascat,manta,tiddit,cnvkit,controlfreec \
   --joint_mutect2 \
   --max_cpus 128 \
   --max_memory '256.GB' \
   --pon /home/project/11003581/Ref/pons/1000g_pon.hg38.vcf.gz \
   --pon_tbi /home/project/11003581/Ref/pons/1000g_pon.hg38.vcf.gz.tbi \
   --email ash.ps@nus.edu.sg \
   -name sarek-mash-b1


