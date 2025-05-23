#!/bin/bash

#PBS -l select=3:ncpus=64:mem=256g
#PBS -l walltime=24:00:00
#PBS -P 11003581
#PBS -N mash-tn-b1
#PBS -j oe

# Change to the directory where the job was submitted
cd $PBS_O_WORKDIR

module load java/11.0.15-openjdk
module load singularity

export NXF_SINGULARITY_CACHEDIR=/home/project/11003581/Tools/singularity-cache/
export SINGULARITY_CACHEDIR=/home/project/11003581/Tools/singularity-cache/

/home/project/11003581/Tools/nextflow run nf-core/sarek -r 3.5.1 \
   -profile singularity \
   -name sarek-TN-b1 \
   -work-dir /scratch/users/nus/ash.ps/NCCS-MASH/analysis/tumor/b1/work \
   --input /home/users/nus/ash.ps/scratch/NCCS-MASH/analysis/tumor/b1/samplesheet.csv \
   --outdir /home/users/nus/ash.ps/scratch/NCCS-MASH/analysis/tumor/b1/ \
   --tools mutect2,manta,strelka \
   --aligner bwa-mem2 \
   --joint_mutect2 \
   --pon /home/project/11003581/Ref/pons/1000g_pon.hg38.vcf.gz \
   --pon_tbi /home/project/11003581/Ref/pons/1000g_pon.hg38.vcf.gz.tbi \
   --igenomes_base /home/project/11003581/Ref/Homo_sapiens/GATK/hg38/


########## re running #########


#!/bin/bash

#PBS -l select=3:ncpus=64:mem=256g
#PBS -l walltime=24:00:00
#PBS -P 11003581
#PBS -N sarek-mash-TO
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load java/17.0.6-jdk
module load singularity
module load nextflow/24.10.5

export NXF_SINGULARITY_CACHEDIR=/home/project/11003581/Tools/singularity-cache/
export SINGULARITY_CACHEDIR=/home/project/11003581/Tools/singularity-cache/
export SINGULARITY_BIND="/home/users/nus/ash.ps,/home/project/11003581"

nextflow run nf-core/sarek -r 3.5.1 \
   -profile singularity \
   --input /home/users/nus/ash.ps/scratch/NCCS-MASH/analysis/tumor/analysis/samplesheet.csv \
   --outdir /home/users/nus/ash.ps/scratch/NCCS-MASH/analysis/tumor/analysis \
   --genome GATK.GRCh38 \
   --tools mutect2 \
   --joint_mutect2 \
   --max_cpus 128 \
   --max_memory '256.GB' \
   --pon /home/project/11003581/Ref/pons/1000g_pon.hg38.vcf.gz \
   --pon_tbi /home/project/11003581/Ref/pons/1000g_pon.hg38.vcf.gz.tbi \
   --email ash.ps@nus.edu.sg \
   -name sarek-mash-tumor-only




#!/bin/bash

#PBS -l select=3:ncpus=64:mem=256g
#PBS -l walltime=24:00:00
#PBS -P 11003581
#PBS -N sarek-mash-TO
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load java/17.0.6-jdk
module load singularity

export NXF_SINGULARITY_CACHEDIR=/home/project/11003581/Tools/singularity-cache/
export SINGULARITY_CACHEDIR=/home/project/11003581/Tools/singularity-cache/

/home/project/11003581/Tools/nextflow run nf-core/sarek -r 3.5.1 \
   -profile singularity \
   --input /home/users/nus/ash.ps/scratch/NCCS-MASH/analysis/tumor/new-analysis/samplesheet.csv \
   --outdir /home/users/nus/ash.ps/scratch/NCCS-MASH/analysis/tumor/new-analysis \
   --genome GATK.GRCh38 \
   --tools mutect2 \
   --joint_mutect2 \
   --max_cpus 128 \
   --max_memory '256.GB' \
   --pon /home/project/11003581/Ref/pons/1000g_pon.hg38.vcf.gz \
   --pon_tbi /home/project/11003581/Ref/pons/1000g_pon.hg38.vcf.gz.tbi \
   --email ash.ps@nus.edu.sg \
   -name sarek-mash-tumor-only