##################### --------------------------------------------- #####################
# RAD51_S181P
##################### --------------------------------------------- #####################

#!/bin/bash

#PBS -l select=1:ncpus=128:mem=256g
#PBS -l walltime=04:00:00
#PBS -P 11003581
#PBS -N slice-bam-rad51
#PBS -j oe

cd $PBS_O_WORKDIR
module load samtools

REFERENCE="/home/project/11003581/Ref/sarek-refs/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"
BASE_INPUT_DIR="/home/users/nus/ash.ps/scratch/sep-fun-analysis/RAD51_S181P/preprocessing/recalibrated"
OUTPUT_DIR="/home/users/nus/ash.ps/scratch/sep-fun-analysis/RAD51_S181P/preprocessing/chr15-bams"
mkdir -p "$OUTPUT_DIR"

find "$BASE_INPUT_DIR" -type f -name "*.cram" | while read -r cramfile; do
  filename=$(basename "$cramfile" .cram)
  outbam="$OUTPUT_DIR/${filename}_chr15.bam"

  # Detect whether contigs are '15' or 'chr15'
  if samtools idxstats "$cramfile" | cut -f1 | grep -qx "chr15"; then
    CHR="chr15"
  elif samtools idxstats "$cramfile" | cut -f1 | grep -qx "15"; then
    CHR="15"
  else
    echo "Could not find chromosome 15 name in $cramfile"; continue
  fi

  echo "Processing $cramfile -> $outbam ($CHR)"
  # Extract region, keep header, sort just to be safe, then index
  samtools view -@ 16 -T "$REFERENCE" -b -h "$cramfile" "$CHR" \
    | samtools sort -@ 16 -o "$outbam" -
  samtools index -@ 16 "$outbam"
done


##################### --------------------------------------------- #####################
# BRCA2_P3292L
##################### --------------------------------------------- #####################

#!/bin/bash

#PBS -l select=1:ncpus=128:mem=256g
#PBS -l walltime=04:00:00
#PBS -P 11003581
#PBS -N slice-bam-brca2
#PBS -j oe

cd $PBS_O_WORKDIR
module load samtools

REFERENCE="/home/project/11003581/Ref/sarek-refs/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"
BASE_INPUT_DIR="/home/users/nus/ash.ps/scratch/sep-fun-analysis/BRCA2_P3292L/preprocessing/recalibrated"
OUTPUT_DIR="/home/users/nus/ash.ps/scratch/sep-fun-analysis/BRCA2_P3292L/preprocessing/chr13-bams"
mkdir -p "$OUTPUT_DIR"

find "$BASE_INPUT_DIR" -type f -name "*.cram" | while read -r cramfile; do
  filename=$(basename "$cramfile" .cram)
  outbam="$OUTPUT_DIR/${filename}_chr13.bam"

  # Detect whether contigs are '13' or 'chr13'
  if samtools idxstats "$cramfile" | cut -f1 | grep -qx "chr13"; then
    CHR="chr13"
  elif samtools idxstats "$cramfile" | cut -f1 | grep -qx "13"; then
    CHR="13"
  else
    echo "Could not find chromosome 13 name in $cramfile"; continue
  fi

  echo "Processing $cramfile -> $outbam ($CHR)"
  # Extract region, keep header, sort just to be safe, then index
  samtools view -@ 16 -T "$REFERENCE" -b -h "$cramfile" "$CHR" \
    | samtools sort -@ 32 -o "$outbam" -
  samtools index -@ 16 "$outbam"
done


##################### --------------------------------------------- #####################
# BRCA1_R133C
##################### --------------------------------------------- #####################

#!/bin/bash

#PBS -l select=1:ncpus=128:mem=256g
#PBS -l walltime=04:00:00
#PBS -P 11003581
#PBS -N slice-bam-brca1
#PBS -j oe

cd $PBS_O_WORKDIR
module load samtools

REFERENCE="/home/project/11003581/Ref/sarek-refs/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"
BASE_INPUT_DIR="/home/users/nus/ash.ps/scratch/sep-fun-analysis/BRCA1_R133C/preprocessing/recalibrated"
OUTPUT_DIR="/home/users/nus/ash.ps/scratch/sep-fun-analysis/BRCA1_R133C/preprocessing/chr17-bams"
mkdir -p "$OUTPUT_DIR"

find "$BASE_INPUT_DIR" -type f -name "*.cram" | while read -r cramfile; do
  filename=$(basename "$cramfile" .cram)
  outbam="$OUTPUT_DIR/${filename}_chr17.bam"

  # Detect whether contigs are '17' or 'chr17'
  if samtools idxstats "$cramfile" | cut -f1 | grep -qx "chr17"; then
    CHR="chr17"
  elif samtools idxstats "$cramfile" | cut -f1 | grep -qx "17"; then
    CHR="17"
  else
    echo "Could not find chromosome 17 name in $cramfile"; continue
  fi

  echo "Processing $cramfile -> $outbam ($CHR)"
  # Extract region, keep header, sort just to be safe, then index
  samtools view -@ 16 -T "$REFERENCE" -b -h "$cramfile" "$CHR" \
    | samtools sort -@ 32 -o "$outbam" -
  samtools index -@ 16 "$outbam"
done



##################### --------------------------------------------- #####################
# H2A_KO
##################### --------------------------------------------- #####################

#!/bin/bash

#PBS -l select=1:ncpus=128:mem=256g
#PBS -l walltime=04:00:00
#PBS -P 11003581
#PBS -N slice-bam-brca1
#PBS -j oe

cd $PBS_O_WORKDIR
module load samtools

REFERENCE="/home/project/11003581/Ref/sarek-refs/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"
BASE_INPUT_DIR="/home/users/nus/ash.ps/scratch/sep-fun-analysis/H2A_KO/preprocessing/recalibrated"
OUTPUT_DIR="/home/users/nus/ash.ps/scratch/sep-fun-analysis/H2A_KO/preprocessing/chr11-bams"
mkdir -p "$OUTPUT_DIR"

find "$BASE_INPUT_DIR" -type f -name "*.cram" | while read -r cramfile; do
  filename=$(basename "$cramfile" .cram)
  outbam="$OUTPUT_DIR/${filename}_chr11.bam"

  # Detect whether contigs are '11' or 'chr11'
  if samtools idxstats "$cramfile" | cut -f1 | grep -qx "chr11"; then
    CHR="chr11"
  elif samtools idxstats "$cramfile" | cut -f1 | grep -qx "11"; then
    CHR="11"
  else
    echo "Could not find chromosome 11 name in $cramfile"; continue
  fi

  echo "Processing $cramfile -> $outbam ($CHR)"
  # Extract region, keep header, sort just to be safe, then index
  samtools view -@ 16 -T "$REFERENCE" -b -h "$cramfile" "$CHR" \
    | samtools sort -@ 32 -o "$outbam" -
  samtools index -@ 16 "$outbam"
done
