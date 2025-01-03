'''
#source: https://gist.github.com/ckandoth/4bccadcacd58aad055ed369a78bf2e7c

#Create and activate a conda environment with VEP, its dependencies, and other related tools:
module load miniforge3
conda create --prefix /home/project/11003581/conda-envs/vep pip -y
conda activate /home/project/11003581/conda-envs/vep
conda install -y -c conda-forge -c bioconda -c defaults ensembl-vep==113.0 htslib==1.20 bcftools==1.20 samtools==1.20 ucsc-liftover==447

vep_install -a cf -s homo_sapiens -y GRCh38 -c /home/project/11003581/Refs/vep --CONVERT

#Download VEPs offline cache for GRCh38, and the reference FASTA:
mkdir /home/project/11003581/Refs/vep
cd /home/project/11003581/Refs/vep
wget https://ftp.ensembl.org/pub/release-113/variation/indexed_vep_cache/homo_sapiens_vep_113_GRCh38.tar.gz
tar -zvf /home/project/11003581/Refs/vep/homo_sapiens_vep_113_GRCh38.tar.gz -C /home/project/11003581/Refs/vep/
wget https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna_index/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna_index/Homo_sapiens.GRCh38.dna.toplevel.fa.gz.fai
wget https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna_index/Homo_sapiens.GRCh38.dna.toplevel.fa.gz.gzi


#Download file from https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna_index/
# reference FASTA which we must bgzip instead of gzip
gzip -d Homo_sapiens.GRCh38.dna.toplevel.fa.gz
/home/project/11003581/Tools/bin/bgzip Homo_sapiens.GRCh38.dna.toplevel.fa

#vcf to maf
#source: https://github.com/mskcc/vcf2maf

export VCF2MAF_URL=`curl -sL https://api.github.com/repos/mskcc/vcf2maf/releases | grep -m1 tarball_url | cut -d\" -f4`
curl -L -o mskcc-vcf2maf.tar.gz $VCF2MAF_URL; tar -zxf mskcc-vcf2maf.tar.gz; cd mskcc-vcf2maf-*
perl vcf2maf.pl --man
perl maf2maf.pl --man

'''

### Running VEP
#!/bin/bash

#PBS -l select=1:ncpus=32
#PBS -l walltime=6:00:00
#PBS -P 11003581
#PBS -N vcf-to-maf
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load miniforge3
conda activate /home/project/11003581/conda-envs/vep

perl /home/project/11003581/Tools/vcf2maf-1.6.22/vcf2maf.pl \
    --tabix-exec=/home/project/11003581/Tools/bin/tabix \
    --input-vcf=/home/users/nus/ash.ps/scratch/YS-analysis/mutect2-vcfs/T01_vs_N01.mutect2.filtered.vcf \
    --output-maf=/home/users/nus/ash.ps/scratch/YS-analysis/vcf-to-maf/P01.maf \
    --ref-fasta=/home/project/11003581/Ref/vep/homo_sapiens_vep_113_GRCh38.tar.gz \
    --vep-path=/home/project/11003581/conda-envs/vep/bin \
    --vep-data=/home/project/11003581/Ref/vep/


######### RUnning loop for every vcf files ##########


#!/bin/bash

#PBS -l select=1:ncpus=32
#PBS -l walltime=6:00:00
#PBS -P 11003581
#PBS -N vcf-to-maf
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load miniforge3
conda activate /home/project/11003581/conda-envs/vep

# Define input and output directories
input_dir="/home/users/nus/ash.ps/scratch/MASH-vcfs/vcf-inputs/"
output_dir="/home/users/nus/ash.ps/scratch/MASH-analysis/mafs"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop through all gzipped VCF files in the input directory
for vcf_gz_file in "$input_dir"/*.vcf; do
    # Extract the base name of the VCF file (without path and extension)
    base_name=$(basename "$vcf_gz_file" .vcf)

    # Construct the output MAF file path
    maf_file="$output_dir/${base_name}.maf"

    # Run vcf2maf for each decompressed VCF file
    perl /home/project/11003581/Tools/vcf2maf-1.6.22/vcf2maf.pl \
        --tabix-exec=/home/project/11003581/Tools/bin/tabix \
        --input-vcf="$input_dir/$base_name.vcf" \
        --output-maf="$maf_file" \
        --ref-fasta=/home/project/11003581/Ref/vep/homo_sapiens_vep_113_GRCh38.tar.gz \
        --vep-path=/home/project/11003581/conda-envs/vep/bin \
        --vep-data=/home/project/11003581/Ref/vep/

done

