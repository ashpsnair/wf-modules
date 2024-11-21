#source: https://gist.github.com/ckandoth/4bccadcacd58aad055ed369a78bf2e7c

#Create and activate a conda environment with VEP, its dependencies, and other related tools:
module load miniforge3
conda update -y -n base -c defaults conda
conda config --set solver libmamba
conda create --prefix /home/project/11003581/conda-envs/ -y 
conda activate vep
conda install -y -c conda-forge -c bioconda -c defaults ensembl-vep==113.0 htslib==1.20 bcftools==1.20 samtools==1.20 ucsc-liftover==447

#Download VEP's offline cache for GRCh38, and the reference FASTA:
mkdir /home/project/11003581/Refs/vep
wget https://ftp.ensembl.org/pub/release-112/variation/indexed_vep_cache/homo_sapiens_vep_112_GRCh38.tar.gz
tar -zvf /home/project/11003581/Refs/vep/homo_sapiens_vep_112_GRCh38.tar.gz -C /home/project/11003581/Refs/vep/


#vcf to maf
#source: https://github.com/mskcc/vcf2maf

export VCF2MAF_URL=`curl -sL https://api.github.com/repos/mskcc/vcf2maf/releases | grep -m1 tarball_url | cut -d\" -f4`
curl -L -o mskcc-vcf2maf.tar.gz $VCF2MAF_URL; tar -zxf mskcc-vcf2maf.tar.gz; cd mskcc-vcf2maf-*
perl vcf2maf.pl --man
perl maf2maf.pl --man



### Running VEP
#!/bin/bash

#PBS -l select=1:ncpus=32
#PBS -l walltime=6:00:00
#PBS -P 11003581
#PBS -N vcf-to-maf
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load samtools 

perl /home/project/11003581/Tools/vcf2maf-1.6.22/vcf2maf.pl \
    --tabix-exec=/home/project/11003581/Tools/bin/tabix \
    --input-vcf=/home/users/nus/ash.ps/scratch/YS-analysis/mutect2-vcfs/T01_vs_N01.mutect2.filtered.vcf \
    --output-maf=/home/users/nus/ash.ps/scratch/YS-analysis/vcf-to-maf/P01.maf \
    --ref-fasta=/home/project/11003581/Ref/vep/homo_sapiens_vep_113_GRCh38.tar.gz