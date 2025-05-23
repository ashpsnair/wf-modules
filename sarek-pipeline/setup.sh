###code working in enambis account
#!/bin/bash

#PBS -l select=3:ncpus=64:mem=256g
#PBS -l walltime=24:00:00
#PBS -P 11003581
#PBS -N sarek-trial
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load nextflow/24.10.5
module load miniforge3

export NXF_CONDA_CACHEDIR=/home/project/11003581/Tools/nf-conda-cache
export NXF_CONDA_USE_MAMBA=true

nextflow run nf-core/sarek -r 3.5.1 \
   -profile conda \
   --input /home/users/nus/enambis/scratch/wgs-1/samplesheet.csv \
   --outdir /home/users/nus/enambis/scratch/wgs-1/ \
   --genome GATK.GRCh38 \
   --tools mutect2 \
   --joint_mutect2 \
   --igenomes_base /home/project/11003581/Ref/sarek-refs/ \
   -resume


## samplesheet.csv that ran
patient,sex,status,sample,lane,fastq_1,fastq_2
patient01,XY,0,N01,lane_1,/home/project/11003581/Data/YS-analysis/raw_data/N01_1.fq.gz,/home/project/11003581/Data/YS-analysis/raw_data/N01_2.fq.gz
patient01,XY,1,T01,lane_1,/home/project/11003581/Data/YS-analysis/raw_data/T01_1.fq.gz,/home/project/11003581/Data/YS-analysis/raw_data/T01_2.fq.gz
patient02,XY,1,T02,lane_1,/home/project/11003581/Data/YS-analysis/raw_data/T02_1.fq.gz,/home/project/11003581/Data/YS-analysis/raw_data/T02_2.fq.gz
patient02,XY,0,N02,lane_1,/home/project/11003581/Data/YS-analysis/raw_data/N02_1.fq.gz,/home/project/11003581/Data/YS-analysis/raw_data/N02_2.fq.gz








#modifications to nextflow.config
vim /home/project/11003581/Tools/assets/nf-core/sarek/nextflow.config 

singularity {
  enabled    = true
  autoMounts = true
  cacheDir   = '/home/project/11003581/Tools/singularity-cache'
  runOptions = '--bind /scratch,/home/project/11003581,/home/users/nus,/scratch/users/nus/enambis/wgs,/home/users/nus/enambis/scratch/wgs --writable-tmpfs'
}

process {
  beforeScript = '''
      module load singularity/3.10.0   # or apptainer
      export PATH=$PATH:/app/apps/singularity/3.10.0/bin
  '''
}
