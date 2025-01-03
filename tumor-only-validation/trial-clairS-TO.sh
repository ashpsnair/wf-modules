'''
####### getting the image and configuring it(steps i did)

module load singularity
# singularity pull docker pre-built image
singularity pull docker://hkubal/clairs:latest

#create environment
module load singularity
module load miniforge3
conda config --add channels defaults
conda create -n singularity-env -c conda-forge singularity -y

'''

#!/bin/bash

#PBS -l select=1:ngpus=1
#PBS -l walltime=10:00:00
#PBS -P personal-ash.ps
#PBS -N trial-clairs
#PBS -j oe

module load singularity
module load miniforge3

conda activate singularity-env


# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

ref_fasta="/home/project/11003581/Ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
output_dir="/home/project/11003581/Data/Ash/P3L-lab-analysis/clairs_output"
samplename="p3292l"
input_dir="/home/project/11003581/Data/Ash/P3L-lab-analysis/sorted_bams/"
p3l_bam="/home/project/11003581/Data/Ash/P3L-lab-analysis/sorted_bams/p3292l_hac_sorted.bam" 
normal_bam="/home/project/11003581/Data/Ash/P3L-lab-analysis/sorted_bams/p3292l_hac_sorted.bam" 

# give entire path of SIF image
SIF_IMAGE="/home/project/11003581/Tools/sifs/clairs_latest.sif"

########
mkdir -p $output_dir

# run the sandbox like this afterward
singularity exec \
  -B ${input_dir},${output_dir} \
  $SIF_IMAGE \
  /opt/bin/run_clairs \
  --bind '/home/project/11003581/Ref/':'/ref/' \
  --normal_bam_fn $normal_bam \
  --tumor_bam_fn $p3l_bam \
  --ref_fn /ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
  --output_dir $output_dir \
  --platform ont_r10_dorado_hac_4khz \
  --conda_prefix /scratch/GPFS/app/apps/miniforge3/23.10 \
  --threads 10




