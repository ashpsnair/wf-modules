

# 2. Activate Miniforge
source /home/project/11003581/Tools/miniforge3/bin/activate

# 3. Create the SCcaller env under your conda-envs folder
conda create -y --prefix /home/project/11003581/conda-envs/sccaller python=2.7

# 4. Activate that env
conda activate /home/project/11003581/conda-envs/sccaller



# Clone repos
git clone https://github.com/biosinodx/SCcaller.git
git clone https://github.com/XiaoDongLab/SCcaller-pipeline.git

# Install the SCcaller entrypoint
cd SCcaller
chmod +x sccaller_v2.0.0.py
ln -s $(pwd)/sccaller_v2.0.0.py $HOME/.local/bin/sccaller


##########################################################################################################################
##sccaller_pipeline_1.sh
###########################################################################################################################

#!/bin/bash

#PBS -l select=3:ncpus=64:mem=256g
#PBS -l walltime=24:00:00
#PBS -P 11003581
#PBS -N run-scDNA
#PBS -j oe

source /home/project/11003581/Tools/miniforge3/bin/activate
conda activate /home/project/11003581/conda-envs/sccaller

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

python /home/users/nus/ash.ps/scratch/scDNA/analysis/logs/sccaller_v2.0.0.py \
    --bam /home/users/nus/ash.ps/scratch/scDNA/DNA/Secondary-Analysis-DNA/2204-A2/2204-A2.sorted.bqsr.dedup.bam \
    --fasta /home/project/11003581/Ref/Homo_sapiens/GATK/hg38/Homo_sapiens_assembly38.fasta \
    --output /home/users/nus/ash.ps/scratch/scDNA/analysis/2204-A2.vcf \
    --bulk /home/users/nus/ash.ps/scratch/scDNA/Bulk-data/Secondary-Analysis/2204/2204.sorted.bqsr.dedup.bam \
    --snp_type dbsnp \
    --snp_in /home/project/11003581/Ref/Homo_sapiens/GATK/hg38/dbsnp_144.hg38.vcf.gz  \
    --cpu_num 128 \
    --engine samtools 



######## Calling somatic SNVs and INDELs not present in bulk DNA

#!/bin/bash

#PBS -l select=3:ncpus=64:mem=256g
#PBS -l walltime=24:00:00
#PBS -P 11003581
#PBS -N run-somatic-scDNA
#PBS -j oe

source /home/project/11003581/Tools/miniforge3/bin/activate
conda activate /home/project/11003581/conda-envs/sccaller

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

python /home/users/nus/ash.ps/scratch/scDNA/analysis/logs/sccaller_v2.0.0.py \
    --bam /home/users/nus/ash.ps/scratch/scDNA/DNA/Secondary-Analysis-DNA/2204-A2/2204-A2.sorted.bqsr.dedup.bam \
    --bulk /home/users/nus/ash.ps/scratch/scDNA/Bulk-data/Secondary-Analysis/2204/2204.sorted.bqsr.dedup.bam \
    --fasta /home/project/11003581/Ref/Homo_sapiens/GATK/hg38/Homo_sapiens_assembly38.fasta \
    --output /home/users/nus/ash.ps/scratch/scDNA/analysis/somatic/somatic-2204-A2.vcf \
    --snp_type hsnp \
    --snp_in /home/users/nus/ash.ps/scratch/scDNA/Bulk-data/Secondary-Analysis/2204/2204.g.vcf.gz \
    --cpu_num 128 \
    --engine samtools