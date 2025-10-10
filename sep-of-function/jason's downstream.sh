################### Merge strelka snvs and indels 
#!/bin/bash

#PBS -l select=1:ncpus=64
#PBS -l walltime=1:00:00
#PBS -P 11003581
#PBS -N intersect
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load bcftools/1.15.1
################### Merge strelka snvs and indels 
#!/bin/bash

#PBS -l select=1:ncpus=64
#PBS -l walltime=1:00:00
#PBS -P 11003581
#PBS -N intersect
#PBS -j oe

cd $PBS_O_WORKDIR
module load bcftools/1.15.1

############### Merge strelka ###############

SNV_DIR="/home/users/nus/ash.ps/scratch/sep-fun-analysis/Jasons/strelka/snv/"
INDEL_DIR="/home/users/nus/ash.ps/scratch/sep-fun-analysis/Jasons/strelka/indel/"
OUTPUT_DIR="/home/users/nus/ash.ps/scratch/sep-fun-analysis/Jasons/strelka/merged/"
mkdir -p "$OUTPUT_DIR"

for snv_file in "$SNV_DIR"*_somatic_strelka2_snp_merged.vcf.gz; do
    sample_name=$(basename "$snv_file" _somatic_strelka2_snp_merged.vcf.gz)
    indel_file="$INDEL_DIR${sample_name}_somatic_strelka2_indel_merged.vcf.gz"

    if [ -f "$indel_file" ]; then
        echo "Concatenating SNV + INDEL for ${sample_name}"
        bcftools concat -a -O z -o "${OUTPUT_DIR}${sample_name}.merged.vcf.gz" "$snv_file" "$indel_file"
        bcftools index "${OUTPUT_DIR}${sample_name}.merged.vcf.gz"
        echo "Merged and indexed: ${sample_name}"
    else
        echo "Warning: No matching indel file found for ${sample_name}"
    fi
done

echo "Merging complete. Results stored in ${OUTPUT_DIR}"

########### Merge mutect2 #############
SNV_DIR="/home/users/nus/ash.ps/scratch/sep-fun-analysis/Jasons/mutect2/snv/"
INDEL_DIR="/home/users/nus/ash.ps/scratch/sep-fun-analysis/Jasons/mutect2/indel/"
OUTPUT_DIR="/home/users/nus/ash.ps/scratch/sep-fun-analysis/Jasons/mutect2/merged/"
mkdir -p "$OUTPUT_DIR"

for snv_file in "$SNV_DIR"*_somatic_mutect2_snp_merged.vcf.gz; do
    sample_name=$(basename "$snv_file" _somatic_mutect2_snp_merged.vcf.gz)
    indel_file="$INDEL_DIR${sample_name}_somatic_mutect2_indel_merged.vcf.gz"

    if [ -f "$indel_file" ]; then
        echo "Concatenating SNV + INDEL for ${sample_name}"
        bcftools concat -a -O z -o "${OUTPUT_DIR}${sample_name}.merged.vcf.gz" "$snv_file" "$indel_file"
        bcftools index "${OUTPUT_DIR}${sample_name}.merged.vcf.gz"
        echo "Merged and indexed: ${sample_name}"
    else
        echo "Warning: No matching indel file found for ${sample_name}"
    fi
done

echo "Merging complete. Results stored in ${OUTPUT_DIR}"



################### intersect strelka and mutect2
#!/bin/bash

#PBS -l select=1:ncpus=64
#PBS -l walltime=1:00:00
#PBS -P 11003581
#PBS -N intersect
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load bcftools/1.15.1

# Define directories
mutect2_dir="/home/users/nus/ash.ps/scratch/sep-fun-analysis/Jasons/mutect2/merged/"
strelka_dir="/home/users/nus/ash.ps/scratch/sep-fun-analysis/Jasons/strelka/merged"

