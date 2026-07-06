#!/bin/bash

#PBS -l select=3:ncpus=64:mem=256g
#PBS -l walltime=10:00:00
#PBS -P 11003581
#PBS -N run-WT-pon
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load java/17.0.6-jdk
module load python/3.12.1-gcc11
module load samtools/1.15.1

# Define full path to executables
bwa_mem2=/home/project/11003581/Tools/bwa-mem2-2.2.1_x64-linux/bwa-mem2
samtools=/app/apps/samtools/1.15.1/bin/samtools

set -euo pipefail

# --- REQUIRED INPUTS ---
REF=/home/project/11003581/Ref/sarek-refs/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta
CRAM_LIST="/home/users/nus/ash.ps/scratch/sep-fun-analysis/WT-pon/wt_crams.txt"        
OUTDIR=/home/users/nus/ash.ps/scratch/sep-fun-analysis/WT-pon
PON_OUT="/home/users/nus/ash.ps/scratch/sep-fun-analysis/WT-pon/WT_clones.PON.vcf.gz"

mkdir -p "$OUTDIR"

# 1) Mutect2 on each WT clone CRAM (normal-only discovery)
while IFS= read -r cram; do
  sample=$(basename "$cram" .cram)
  gatk Mutect2 \
    -R "$REF" \
    -I "$cram" \
    --max-mnp-distance 0 \
    --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
    -O "$OUTDIR/${sample}.unfiltered.vcf.gz"
done < "$CRAM_LIST"

# 2) Create the Panel of Normals (PON)
gatk CreateSomaticPanelOfNormals \
  $(for v in "$OUTDIR"/*.unfiltered.vcf.gz; do echo -V "$v"; done) \
  -O "$PON_OUT"

# 3) (Optional) index the PON
gatk IndexFeatureFile -I "$PON_OUT"

echo "Done. PON: $PON_OUT"







############## samples list

/home/users/nus/ash.ps/scratch/sep-fun-analysis/WT/b1/preprocessing/recalibrated/WT_P30_10/WT_P30_10.recal.cram
/home/users/nus/ash.ps/scratch/sep-fun-analysis/WT/b1/preprocessing/recalibrated/WT_P30_4/WT_P30_4.recal.cram
/home/users/nus/ash.ps/scratch/sep-fun-analysis/WT/b1/preprocessing/recalibrated/WT_P30_5/WT_P30_5.recal.cram
/home/users/nus/ash.ps/scratch/sep-fun-analysis/WT/b1/preprocessing/recalibrated/WT_P30_6/WT_P30_6.recal.cram
/home/users/nus/ash.ps/scratch/sep-fun-analysis/WT/b1/preprocessing/recalibrated/WT_P30_7/WT_P30_7.recal.cram
/home/users/nus/ash.ps/scratch/sep-fun-analysis/WT/b1/preprocessing/recalibrated/WT_P30_8/WT_P30_8.recal.cram
/home/users/nus/ash.ps/scratch/sep-fun-analysis/WT/b1/preprocessing/recalibrated/WT_P30_9/WT_P30_9.recal.cram
