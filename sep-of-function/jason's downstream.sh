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

# ====== CONFIG (edit paths if needed) ======
# 1) Your existing inputs for per-caller SNV/INDEL:
STRELKA_SNV_DIR="/home/users/nus/ash.ps/scratch/sep-fun-analysis/Jasons/strelka/snv/"
STRELKA_INDEL_DIR="/home/users/nus/ash.ps/scratch/sep-fun-analysis/Jasons/strelka/indel/"
STRELKA_MERGED_DIR="/home/users/nus/ash.ps/scratch/sep-fun-analysis/Jasons/strelka/merged/"

M2_SNV_DIR="/home/users/nus/ash.ps/scratch/sep-fun-analysis/Jasons/mutect2/snv/"
M2_INDEL_DIR="/home/users/nus/ash.ps/scratch/sep-fun-analysis/Jasons/mutect2/indel/"
M2_MERGED_DIR="/home/users/nus/ash.ps/scratch/sep-fun-analysis/Jasons/mutect2/merged/"

# 2) Reference FASTA for normalization (must have .fai index)
REF="/home/project/11003581/refs/hg38.fa"

# 3) Outputs (and temp work)
WORKDIR="/home/users/nus/ash.ps/scratch/sep-fun-analysis/Jasons/conc_work/"
SUMMARY_TSV="/home/users/nus/ash.ps/scratch/sep-fun-analysis/Jasons/concordance_summary.tsv"

THREADS=64
PASS_ONLY=true   # set to true to count only PASS variants
# ==========================================

mkdir -p "$STRELKA_MERGED_DIR" "$M2_MERGED_DIR" "$WORKDIR"

# ---------- helpers ----------
ensure_bgzip_index () {
  local f="$1"
  if ! bcftools view -h "$f" >/dev/null 2>&1; then
    local tmp="${f%.vcf.gz}.rebgz.vcf.gz"
    bcftools view -O z -o "$tmp" "$f"
    mv "$tmp" "$f"
  fi
  [ -f "${f}.tbi" ] || bcftools index -t "$f"
}

normalize_split () {
  # Left-align against REF and split multiallelics
  local in="$1" out="$2"
  ensure_bgzip_index "$in"
  bcftools norm -f "$REF" -m -both -O z -o "$out" "$in"
  bcftools index -t "$out"
}

count_records () {
  # Count records; optionally only PASS
  local in="$1"
  if $PASS_ONLY; then
    bcftools view -f PASS -H "$in" | wc -l
  else
    bcftools view -H "$in" | wc -l
  fi
}

# ---------- 1) CONCAT SNV + INDEL (same sample) for STRELKA ----------
echo "[Strelka] Concatenating SNV + INDEL per sample ..."
shopt -s nullglob
for snv_file in "${STRELKA_SNV_DIR}"*_somatic_strelka2_snp_merged.vcf.gz; do
  sample_name="$(basename "$snv_file" _somatic_strelka2_snp_merged.vcf.gz)"
  indel_file="${STRELKA_INDEL_DIR}${sample_name}_somatic_strelka2_indel_merged.vcf.gz"
  out_file="${STRELKA_MERGED_DIR}${sample_name}.merged.vcf.gz"
  if [ -f "$indel_file" ]; then
    bcftools concat -a -O z -o "$out_file" "$snv_file" "$indel_file"
    bcftools index -t "$out_file"
    echo "[Strelka] OK: $sample_name"
  else
    echo "[Strelka] WARN: no matching INDEL for $sample_name"
  fi
done

# ---------- 2) CONCAT SNV + INDEL (same sample) for MUTECT2 ----------
echo "[Mutect2] Concatenating SNV + INDEL per sample ..."
for snv_file in "${M2_SNV_DIR}"*_somatic_mutect2_snp_merged.vcf.gz; do
  sample_name="$(basename "$snv_file" _somatic_mutect2_snp_merged.vcf.gz)"
  indel_file="${M2_INDEL_DIR}${sample_name}_somatic_mutect2_indel_merged.vcf.gz"
  out_file="${M2_MERGED_DIR}${sample_name}.merged.vcf.gz"
  if [ -f "$indel_file" ]; then
    bcftools concat -a -O z -o "$out_file" "$snv_file" "$indel_file"
    bcftools index -t "$out_file"
    echo "[Mutect2] OK: $sample_name"
  else
    echo "[Mutect2] WARN: no matching INDEL for $sample_name"
  fi
done

# ---------- 3) PER-SAMPLE CONCORDANCE (intersection) + COUNTS ----------
echo -e "sample\tmutect_total\tstrelka_total\tconcordant\tmutect_only\tstrelka_only" > "$SUMMARY_TSV"

for st_vcf in "${STRELKA_MERGED_DIR}"*.merged.vcf.gz; do
  sample="$(basename "$st_vcf" .merged.vcf.gz)"
  m2_vcf="${M2_MERGED_DIR}${sample}.merged.vcf.gz"
  if [ ! -f "$m2_vcf" ]; then
    echo "[Concordance] WARN: missing Mutect2 for $sample; skipping"
    continue
  fi

  # Normalize both before comparing
  st_norm="${WORKDIR}${sample}.strelka.norm.vcf.gz"
  m2_norm="${WORKDIR}${sample}.mutect2.norm.vcf.gz"
  normalize_split "$st_vcf" "$st_norm"
  normalize_split "$m2_vcf" "$m2_norm"

  # Counts per caller
  st_total=$(count_records "$st_norm")
  m2_total=$(count_records "$m2_norm")

  # Intersection & uniques
  tmpdir="$(mktemp -d -p "$WORKDIR" isec_${sample}_XXXX)"
  # -n=2 ensures output of shared sites; directory contains:
  # 0000 = st_only, 0001 = m2_only, 0002 = shared (A∩B; representation from st_norm)
  if $PASS_ONLY; then
    bcftools isec -f PASS -p "$tmpdir" -Oz "$st_norm" "$m2_norm"
  else
    bcftools isec -p "$tmpdir" -Oz "$st_norm" "$m2_norm"
  fi

  st_only_vcf="$tmpdir/0000.vcf.gz"
  m2_only_vcf="$tmpdir/0001.vcf.gz"
  both_vcf="$tmpdir/0002.vcf.gz"

  # Some of these may be empty/not created; guard counts
  [[ -f "$both_vcf" ]] && both_cnt=$(bcftools view -H "$both_vcf" | wc -l) || both_cnt=0
  [[ -f "$m2_only_vcf" ]] && m2_only_cnt=$(bcftools view -H "$m2_only_vcf" | wc -l) || m2_only_cnt=0
  [[ -f "$st_only_vcf" ]] && st_only_cnt=$(bcftools view -H "$st_only_vcf" | wc -l) || st_only_cnt=0

  echo -e "${sample}\t${m2_total}\t${st_total}\t${both_cnt}\t${m2_only_cnt}\t${st_only_cnt}" >> "$SUMMARY_TSV"

  rm -rf "$tmpdir"
  echo "[Concordance] ${sample} -> m2:${m2_total} st:${st_total} both:${both_cnt} m2_only:${m2_only_cnt} st_only:${st_only_cnt}"
done

echo "[DONE] Summary written to: $SUMMARY_TSV"
