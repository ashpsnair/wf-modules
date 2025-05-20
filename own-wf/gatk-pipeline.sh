#!/bin/bash

#PBS -l select=3:ncpus=64:mem=256g
#PBS -l walltime=24:00:00
#PBS -P 11003581
#PBS -N run-workflow
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load java/17.0.6-jdk
module load python/3.12.1-gcc11


# Define full path to executables
bwa_mem2=/home/project/11003581/Tools/bwa-mem2-2.2.1_x64-linux/bwa-mem2
samtools=/app/apps/samtools/1.15.1/bin/samtools

# Define variables
ref_fasta=/home/project/11003581/Ref/Homo_sapiens/GATK/hg38/Homo_sapiens_assembly38.fasta
t_f1=/home/users/nus/ash.ps/scratch/NCCS-MASH/analysis/tumor/raw-fastqs/concat-fastqs/b1/WHT463_R1.fastq.gz
t_f2=/home/users/nus/ash.ps/scratch/NCCS-MASH/analysis/tumor/raw-fastqs/concat-fastqs/b1/WHT463_R2.fastq.gz
n_bam=/home/users/nus/ash.ps/scratch/NCCS-MASH/analysis/normal/b3/preprocessing/recalibrated/WHT462/WHT462.recal.cram
output_dir=/home/users/nus/ash.ps/scratch/NCCS-MASH/FINAL/somatic

# Known sites for BQSR (update these paths as needed)
dbsnp_vcf=/home/project/11003581/Ref/Homo_sapiens/GATK/hg38/dbsnp_146.hg38.vcf.gz
mills_vcf=/home/project/11003581/Ref/Homo_sapiens/GATK/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

# Resources for Mutect2 (update these paths as needed)
gnomad_vcf=/home/project/11003581/Ref/Homo_sapiens/GATK/hg38/af-only-gnomad.hg38.vcf.gz
pon_vcf=/home/project/11003581/Ref/pons/1000g_pon.hg38.vcf.gz

# Output file names
tumor_sorted_bam=${output_dir}/tumor.sorted.bam
tumor_marked_bam=${output_dir}/tumor.marked.bam
tumor_metrics=${output_dir}/tumor.metrics.txt
tumor_recal_table=${output_dir}/tumor.recal.table
tumor_final_bam=${output_dir}/tumor.final.bam

# Get normal sample name from CRAM header
normal_sample_name=$($samtools view -H "$n_bam" | grep '^@RG' | head -1 | sed -n 's/.*SM:\([^\t]*\).*/\1/p')

mkdir -p "$output_dir"

########## Alignment with BWA-MEM2 (Tumor only) ########
"$bwa_mem2" mem -t 16 "$ref_fasta" "$t_f1" "$t_f2" | \
    $samtools sort -@4 -o "$tumor_sorted_bam" -

if [[ ! -s "$tumor_sorted_bam" ]]; then
    echo "ERROR: $tumor_sorted_bam not created. Check BWA-MEM2 and samtools."
    exit 1
fi

########## Mark Duplicates (Tumor only) ########
gatk MarkDuplicates \
  -I "$tumor_sorted_bam" \
  -O "$tumor_marked_bam" \
  -M "$tumor_metrics"

if [[ ! -s "$tumor_marked_bam" ]]; then
    echo "ERROR: $tumor_marked_bam not created. Check MarkDuplicates step."
    exit 1
fi

########## Base Quality Score Recalibration (Tumor only) ########
gatk BaseRecalibrator \
  -R "$ref_fasta" \
  -I "$tumor_marked_bam" \
  --known-sites "$dbsnp_vcf" \
  --known-sites "$mills_vcf" \
  -O "$tumor_recal_table"

gatk ApplyBQSR \
  -R "$ref_fasta" \
  -I "$tumor_marked_bam" \
  --bqsr-recal-file "$tumor_recal_table" \
  -O "$tumor_final_bam"

if [[ ! -s "$tumor_final_bam" ]]; then
    echo "ERROR: $tumor_final_bam not created. Check ApplyBQSR step."
    exit 1
fi

########## Somatic Variant Calling with Mutect2 ########
gatk Mutect2 \
  -R "$ref_fasta" \
  -I "$tumor_final_bam" \
  -I "$n_bam" \
  -normal "$normal_sample_name" \
  --germline-resource "$gnomad_vcf" \
  --panel-of-normals "$pon_vcf" \
  -O "${output_dir}/somatic.unfiltered.vcf.gz"

########## Filter Variants ########
gatk FilterMutectCalls \
  -R "$ref_fasta" \
  -V "${output_dir}/somatic.unfiltered.vcf.gz" \
  -O "${output_dir}/somatic.filtered.vcf.gz"
