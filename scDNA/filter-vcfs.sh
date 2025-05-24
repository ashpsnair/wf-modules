#!/bin/bash
# ── PBS DIRECTIVES ───────────────────────────────────────────────────────────────
#PBS -N scdna-filter
#PBS -l select=2:ncpus=64:mem=128gb
#PBS -l walltime=10:00:00
#PBS -P 11003581
#PBS -j oe

cd "$PBS_O_WORKDIR"
module load bcftools/1.15.1

BGZIP=/home/project/11003581/Tools/bin/bgzip
TABIX=/home/project/11003581/Tools/bin/tabix

SRC_ROOT=/home/users/nus/ash.ps/scratch/scDNA/analysis
DEST_ROOT=$SRC_ROOT/filtered-vcfs
mkdir -p "$DEST_ROOT"/{2204,3401}

for vcf in \
    $SRC_ROOT/2204-somatic/b*/2204*/*.vcf{,.gz} \
    $SRC_ROOT/3401-somatic/b*/3401*/*.vcf{,.gz}
do
  [[ -e "$vcf" ]] || continue

  sample=$(basename "${vcf%.vcf*}")              # strip .vcf or .vcf.gz
  dest_dir="$DEST_ROOT"/$( [[ $vcf == *"/2204-"* ]] && echo 2204 || echo 3401 )
  out_vcf="$dest_dir/${sample}_filtered.vcf.gz"

  echo "▶ filtering $(basename "$vcf") → $(basename "$out_vcf")"

  # Clean stream: remove REF=N and FILTER≠PASS
  ( [[ $vcf == *.gz ]] && zcat "$vcf" || cat "$vcf" ) \
  | awk 'BEGIN {OFS="\t"}
         /^#/ {print; next}
         $4 != "N" && $7 == "PASS" {print}' \
  | "$BGZIP" -c > "$out_vcf"

  "$TABIX" -p vcf "$out_vcf"
done

echo "✔ Done. Files saved in $DEST_ROOT/{2204,3401}"
