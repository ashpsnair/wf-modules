#!/bin/bash
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

REPORT_2204="$DEST_ROOT/2204_variant_counts.txt"
REPORT_3401="$DEST_ROOT/3401_variant_counts.txt"
echo -e "Cell\tBefore\tAfter" > "$REPORT_2204"
echo -e "Cell\tBefore\tAfter" > "$REPORT_3401"

for vcf in \
    $SRC_ROOT/2204-somatic/b*/2204*/*.vcf{,.gz} \
    $SRC_ROOT/3401-somatic/b*/3401*/*.vcf{,.gz}
do
  [[ -e "$vcf" ]] || continue

  sample=$(basename "${vcf%.vcf*}")
  is_2204=false
  if [[ $vcf == *"/2204-"* ]]; then
    dest_dir=$DEST_ROOT/2204
    report_file=$REPORT_2204
    is_2204=true
  else
    dest_dir=$DEST_ROOT/3401
    report_file=$REPORT_3401
  fi

  out_vcf="$dest_dir/${sample}_filtered.vcf.gz"
  echo "▶ filtering $(basename "$vcf") → $(basename "$out_vcf")"

  # count variants before filtering (excluding headers)
  before=$( ( [[ $vcf == *.gz ]] && zcat "$vcf" || cat "$vcf" ) | grep -vc '^#' )

  # filter and count after removing REF = N
  filtered=$( ( [[ $vcf == *.gz ]] && zcat "$vcf" || cat "$vcf" ) \
    | awk 'BEGIN{c=0} /^#/ {next} $4 != "N" {c++} END{print c}' )

  # apply actual filtering and compress
  ( [[ $vcf == *.gz ]] && zcat "$vcf" || cat "$vcf" ) \
    | awk 'BEGIN {OFS="\t"} /^#/ {print; next} $4 != "N" {print}' \
    | "$BGZIP" -c > "$out_vcf"

  "$TABIX" -p vcf "$out_vcf"

  echo -e "${sample}\t${before}\t${filtered}" >> "$report_file"
done

echo "✔ Done. Filtered VCFs + reports saved in $DEST_ROOT/{2204,3401}"
