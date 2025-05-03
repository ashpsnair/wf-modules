

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

###########################################################################################################################

#!/bin/bash
#PBS -l select=3:ncpus=64:mem=256g
#PBS -l walltime=24:00:00
#PBS -P 11003581
#PBS -N mash-tn-b1
#PBS -j oe

set -euo pipefail

# Change to the directory where the job was submitted
cd $PBS_O_WORKDIR

### CONFIGURATION ###
# Miniforge & conda/env paths
MINIFORGE=/home/project/11003581/Tools/miniforge3/bin/activate
CONDA_ENV=/home/project/11003581/conda-envs/sccaller

# SCcaller script
SCCALLER_SCRIPT=/home/project/11003581/Tools/sccaller/SCcaller/sccaller_v2.0.0.py

# Reference & known‐variant VCFs
REF_FA=/home/project/11003581/Ref/Homo_sapiens/GATK/hg38/Homo_sapiens_assembly38.fasta
DBSNP_VCF=/home/project/11003581/Ref/Homo_sapiens/GATK/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz

# Input FASTQ parent directory and output base
FASTQ_ROOT=/home/users/nus/ash.ps/scratch/scDNA/DNA/Clean-FQ
OUT_ROOT=/home/project/11003581/results/sccaller

# Threads
THREADS=128

### PREP ###
mkdir -p "$OUT_ROOT"
source "$MINIFORGE"
conda activate "$CONDA_ENV"

# (only needed once if ref not yet indexed)
if [ ! -f "${REF_FA}.fai" ]; then
  bwa index "$REF_FA"
  samtools faidx "$REF_FA"
fi

### PROCESS EACH CELL ###
for sample_dir in "$FASTQ_ROOT"/*/ ; do
  sample=$(basename "$sample_dir")
  echo ">>> Processing $sample"

  # locate paired FASTQs (assumes *_1.fq.gz and *_2.fq.gz)
  R1=$(ls "$sample_dir"/*_1.fq.gz 2>/dev/null | head -n1)
  R2=$(ls "$sample_dir"/*_2.fq.gz 2>/dev/null | head -n1)
  if [[ -z "$R1" || -z "$R2" ]]; then
    echo "  ⚠️  FASTQ pairs not found in $sample_dir – skipping"
    continue
  fi

  # 1) Alignment → sorted BAM
  BAM="$OUT_ROOT/${sample}.bam"
  echo "  • Aligning → $BAM"
  bwa mem -M -t "$THREADS" "$REF_FA" "$R1" "$R2" \
    | samtools view -bSh - \
    | samtools sort -@ "$THREADS" -o "$BAM"

  # 2) Index
  echo "  • Indexing BAM"
  samtools index "$BAM"

  # 3) SCcaller: SNVs + INDELs
  VCF="$OUT_ROOT/${sample}.vcf"
  echo "  • Calling variants → $VCF"
  python "$SCCALLER_SCRIPT" \
    --bam "$BAM" \
    --fasta "$REF_FA" \
    --output "$VCF" \
    --snp_type dbsnp \
    --snp_in "$DBSNP_VCF" \
    --cpu_num "$THREADS" \
    --engine samtools

  echo "  ✓ Done $sample"
done

echo "All samples processed."
