################################ PD-DM ################################

#Gzipping batch1

#!/bin/bash

#PBS -l select=3:ncpus=64:mem=64g:ngpus=1:mem=128g
#PBS -l walltime=24:00:00
#PBS -P 11003581
#PBS -N run-gzip-PD-DM-b1
#PBS -j oe

cd /home/users/nus/ash.ps/scratch/T2DM/PD-DM/fastq/

gzip SRR29413850_pass_*.fastq SRR29729001_pass_*.fastq SRR29728999_pass_*.fastq SRR29728998_pass_*.fastq SRR29728997_pass_*.fastq SRR29728996_pass_*.fastq

############ SAREK

#!/bin/bash

#PBS -l select=3:ncpus=64:mem=256g
#PBS -l walltime=24:00:00
#PBS -P 11003581
#PBS -N sarek-PD-DM-b1
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR


module load java/17.0.6-jdk
module load singularity

/home/project/11003581/Tools/nextflow run nf-core/sarek -r 3.4.4 \
   -profile singularity \
   --input /home/users/nus/ash.ps/scratch/T2DM/analysis-PD-DM/samplesheet.csv \
   --outdir /home/users/nus/ash.ps/scratch/T2DM/analysis-PD-DM \
   --tools mutect2,manta \
   --max_cpus 128 \
   --max_memory '256.GB' \
   --pon /home/project/11003581/Ref/pons/1000g_pon.hg38.vcf.gz \
   --pon_tbi /home/project/11003581/Ref/pons/1000g_pon.hg38.vcf.gz.tbi


##amplesheet.csv trial 1
patient,status,sample,lane,fastq_1,fastq_2
SRR29413850,1,SRR29413850,lane_1,/home/users/nus/ash.ps/scratch/T2DM/PD-DM/fastq/SRR29413850_1.fastq.gz,/home/users/nus/ash.ps/scratch/T2DM/PD-DM/fastq/SRR29413850_2.fastq.gz
SRR29729001,1,SRR29729001,lane_1,/home/users/nus/ash.ps/scratch/T2DM/PD-DM/fastq/SRR29729001_1.fastq.gz,/home/users/nus/ash.ps/scratch/T2DM/PD-DM/fastq/SRR29729001_2.fastq.gz

#samplesheet.csv batch1
patient,status,sample,lane,fastq_1,fastq_2
SRR29728999,1,SRR29728999,lane_1,/home/users/nus/ash.ps/scratch/T2DM/PD-DM/fastq/SRR29728999_1.fastq.gz,/home/users/nus/ash.ps/scratch/T2DM/PD-DM/fastq/SRR29728999_2.fastq.gz
SRR29728998,1,SRR29728998,lane_1,/home/users/nus/ash.ps/scratch/T2DM/PD-DM/fastq/SRR29728998_1.fastq.gz,/home/users/nus/ash.ps/scratch/T2DM/PD-DM/fastq/SRR29728998_2.fastq.gz
SRR29728997,1,SRR29728997,lane_1,/home/users/nus/ash.ps/scratch/T2DM/PD-DM/fastq/SRR29728997_1.fastq.gz,/home/users/nus/ash.ps/scratch/T2DM/PD-DM/fastq/SRR29728997_2.fastq.gz

#samplesheet.csv batch2
patient,status,sample,lane,fastq_1,fastq_2
SRR29728996,1,SRR29728996,lane_1,/home/users/nus/ash.ps/scratch/T2DM/PD-DM/fastq/SRR29728996_1.fastq.gz,/home/users/nus/ash.ps/scratch/T2DM/PD-DM/fastq/SRR29728996_2.fastq.gz
SRR29413852,1,SRR29413852,lane_1,/home/users/nus/ash.ps/scratch/T2DM/PD-DM/fastq/SRR29413852_1.fastq.gz,/home/users/nus/ash.ps/scratch/T2DM/PD-DM/fastq/SRR29413852_2.fastq.gz
SRR29413851,1,SRR29413851,lane_1,/home/users/nus/ash.ps/scratch/T2DM/PD-DM/fastq/SRR29413851_1.fastq.gz,/home/users/nus/ash.ps/scratch/T2DM/PD-DM/fastq/SRR29413851_2.fastq.gz
SRR29729002,1,SRR29729002,lane_1,/home/users/nus/ash.ps/scratch/T2DM/PD-DM/fastq/SRR29729002_1.fastq.gz,/home/users/nus/ash.ps/scratch/T2DM/PD-DM/fastq/SRR29729002_2.fastq.gz


################################ PD ################################
#Gzipping batch1

#!/bin/bash

#PBS -l select=3:ncpus=64:mem=64g:mem=128g
#PBS -l walltime=24:00:00
#PBS -P 11003581
#PBS -N run-gzip-PD-b1
#PBS -j oe

cd /home/users/nus/ash.ps/scratch/T2DM/PD-DM/fastq/

gzip SRR29413853_pass_*.fastq SRR29413854_pass_*.fastq SRR29728994_pass_*.fastq SRR29728993_pass_*.fastq SRR29728992_pass_*.fastq SRR29728991_pass_*.fastq SRR29728990_pass_*.fastq SRR29728989_pass_*.fastq

############ SAREK

#!/bin/bash

#PBS -l select=3:ncpus=64:mem=256g
#PBS -l walltime=24:00:00
#PBS -P 11003581
#PBS -N sarek-PD-b1
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR


module load java/17.0.6-jdk
module load singularity

/home/project/11003581/Tools/nextflow run nf-core/sarek -r 3.4.4 \
   -profile singularity \
   --input /home/users/nus/ash.ps/scratch/T2DM/analysis-PD/samplesheet.csv \
   --outdir /home/users/nus/ash.ps/scratch/T2DM/analysis-PD \
   --tools mutect2,manta \
   --max_cpus 128 \
   --max_memory '256.GB' \
   --pon /home/project/11003581/Ref/pons/1000g_pon.hg38.vcf.gz \
   --pon_tbi /home/project/11003581/Ref/pons/1000g_pon.hg38.vcf.gz.tbi

#  PD trial
patient,status,sample,lane,fastq_1,fastq_2
SRR29413853,1,SRR29413853,lane_1,/home/users/nus/ash.ps/scratch/T2DM/PD/fastq/SRR29413853_pass_1.fastq.gz,/home/users/nus/ash.ps/scratch/T2DM/PD/fastq/SRR29413853_pass_2.fastq.gz
SRR29413854,1,SRR29413854,lane_1,/home/users/nus/ash.ps/scratch/T2DM/PD/fastq/SRR29413854_pass_1.fastq.gz,/home/users/nus/ash.ps/scratch/T2DM/PD/fastq/SRR29413854_pass_2.fastq.gz



#PD batch1
patient,status,sample,lane,fastq_1,fastq_2
SRR29413858,1,SRR29413858,lane_1,/home/users/nus/ash.ps/scratch/T2DM/PD/fastq/SRR29413858_1.fastq.gz,/home/users/nus/ash.ps/scratch/T2DM/PD/fastq/SRR29413858_2.fastq.gz
SRR29728988,1,SRR29728988,lane_1,/home/users/nus/ash.ps/scratch/T2DM/PD/fastq/SRR29728988_1.fastq.gz,/home/users/nus/ash.ps/scratch/T2DM/PD/fastq/SRR29728988_2.fastq.gz
SRR29728989,1,SRR29728989,lane_1,/home/users/nus/ash.ps/scratch/T2DM/PD/fastq/SRR29728989_1.fastq.gz,/home/users/nus/ash.ps/scratch/T2DM/PD/fastq/SRR29728989_2.fastq.gz
SRR29728990,1,SRR29728990,lane_1,/home/users/nus/ash.ps/scratch/T2DM/PD/fastq/SRR29728990_1.fastq.gz,/home/users/nus/ash.ps/scratch/T2DM/PD/fastq/SRR29728990_2.fastq.gz
SRR29728991,1,SRR29728991,lane_1,/home/users/nus/ash.ps/scratch/T2DM/PD/fastq/SRR29728991_1.fastq.gz,/home/users/nus/ash.ps/scratch/T2DM/PD/fastq/SRR29728991_2.fastq.gz


#PD batch2
patient,status,sample,lane,fastq_1,fastq_2
SRR29728992,1,SRR29728992,lane_1,/home/users/nus/ash.ps/scratch/T2DM/PD/fastq/SRR29728992_1.fastq.gz,/home/users/nus/ash.ps/scratch/T2DM/PD/fastq/SRR29728992_2.fastq.gz
SRR29728993,1,SRR29728993,lane_1,/home/users/nus/ash.ps/scratch/T2DM/PD/fastq/SRR29728993_1.fastq.gz,/home/users/nus/ash.ps/scratch/T2DM/PD/fastq/SRR29728993_2.fastq.gz
SRR29728994,1,SRR29728994,lane_1,/home/users/nus/ash.ps/scratch/T2DM/PD/fastq/SRR29728994_1.fastq.gz,/home/users/nus/ash.ps/scratch/T2DM/PD/fastq/SRR29728994_2.fastq.gz
SRR29729003,1,SRR29729003,lane_1,/home/users/nus/ash.ps/scratch/T2DM/PD/fastq/SRR29729003_1.fastq.gz,/home/users/nus/ash.ps/scratch/T2DM/PD/fastq/SRR29729003_2.fastq.gz
SRR29729004,1,SRR29729004,lane_1,/home/users/nus/ash.ps/scratch/T2DM/PD/fastq/SRR29729004_1.fastq.gz,/home/users/nus/ash.ps/scratch/T2DM/PD/fastq/SRR29729004_2.fastq.gz