#!/bin/bash
#PBS -l select=1:ncpus=64:mem=128g
#PBS -l walltime=04:00:00
#PBS -P 11003581
#PBS -N batch2-annotate-filter
#PBS -j oe

set -euo pipefail
cd "$PBS_O_WORKDIR"

module load bcftools/1.15.1
module load r/4.2.0

########################
# PATHS (batch2)
########################
BASE_DIR="/home/users/nus/ash.ps/scratch/FINAL-sepfunc/batch2"
RAW_DIR="${BASE_DIR}/unfilt-vcfs"               # has *.concordant.vcf.gz
FILT_DIR="${BASE_DIR}/filtered-vcfs"            # output for FILTER=PASS
ANN_DIR="${BASE_DIR}/annotated"                 # per-sample annovar outputs
POP_VCF_DIR="${BASE_DIR}/pop-filter-vcfs"       # freq-filtered VCFs
POP_TXT_DIR="${BASE_DIR}/pop-filter-multianno"  # freq-filtered multianno TXT
MAF_DIR="${BASE_DIR}/mafs"

mkdir -p "$FILT_DIR" "$ANN_DIR" "$POP_VCF_DIR" "$POP_TXT_DIR" "$MAF_DIR"

########################
# 1) FILTER RAW VCFs (FILTER=PASS)
########################
echo "== Filtering VCFs from: $RAW_DIR -> $FILT_DIR"

# Expecting filenames like: B_68_01_1.concordant.vcf.gz
find "$RAW_DIR" -type f -name "*.concordant.vcf.gz" | while read -r file; do
  bn=$(basename "$file" .concordant.vcf.gz)   # sample base name
  out="${FILT_DIR}/${bn}.PASS.vcf.gz"
  echo "  Filtering $bn ..."
  # Keep only FILTER=PASS; output bgzip + index
  bcftools view -i 'FILTER="PASS"' -Oz -o "$out" "$file"
  bcftools index -f "$out"
done

echo "== Filtering complete."

########################
# 2) ANNOVAR (per sample)
########################
echo "== Running ANNOVAR -> $ANN_DIR"

HUMANDB="/home/project/11003581/Tools/annovar/humandb"
ANNOVAR_PL="/home/project/11003581/Tools/annovar/table_annovar.pl"

# Protocols & operations (same as your earlier run)
PROT="refGene,cosmic70,exac03,gnomad_genome,ALL.sites.2015_08,SAS.sites.2015_08,EAS.sites.2015_08,esp6500siv2_all,avsnp151,icgc28,dbnsfp30a,intervar_20180118,clinvar_20240611"
OPER="g,f,f,f,f,f,f,f,f,f,f,f,f"

# Loop over filtered PASS VCFs
find "$FILT_DIR" -type f -name "*.PASS.vcf.gz" | while read -r vcf_gz; do
  samp=$(basename "$vcf_gz" .PASS.vcf.gz)
  outdir="${ANN_DIR}/${samp}"
  mkdir -p "$outdir"
  echo "  ANNOVAR: $samp"

  # ANNOVAR prefers uncompressed VCF input when -vcfinput + -polish; either is ok but we can stream
  # Decompress to temp (safer for polish), clean after
  tmp_vcf="${outdir}/${samp}.PASS.tmp.vcf"
  bcftools view -Ov -o "$tmp_vcf" "$vcf_gz"

  perl "$ANNOVAR_PL" "$tmp_vcf" \
    "$HUMANDB" \
    -buildver hg38 \
    -out "${outdir}/${samp}" \
    -protocol "$PROT" \
    -operation "$OPER" \
    -remove -vcfinput -polish -nastring .

  rm -f "$tmp_vcf"
done

echo "== ANNOVAR complete."

########################
# 3) POP FREQ FILTER on ANNOVAR .vcf
########################
echo "== Pop-filter VCFs -> $POP_VCF_DIR"

# Find all *.hg38_multianno.vcf in per-sample subfolders
find "$ANN_DIR" -type f -name "*.hg38_multianno.vcf" | while read -r vcf_file; do
  # e.g., annotated/<sample>/<sample>.hg38_multianno.vcf -> <sample>_pop_filt.vcf
  sample=$(basename "$vcf_file" .hg38_multianno.vcf)
  out="${POP_VCF_DIR}/${sample}_pop_filt.vcf"

  bcftools view "$vcf_file" | awk '
  BEGIN {FS="\t"; OFS="\t"}
  /^#/ {print; next}
  {
      split($8, info, ";");
      exac_all="."; exac_eas="."; exac_sas="."; gnomad_genome_all="."; gnomad_genome_eas="."; all_sites_2015_08="."; esp6500siv2_all=".";
      for (i in info) {
          split(info[i], pair, "=");
          if (pair[1] == "ExAC_ALL") exac_all = pair[2];
          if (pair[1] == "ExAC_EAS") exac_eas = pair[2];
          if (pair[1] == "ExAC_SAS") exac_sas = pair[2];
          if (pair[1] == "gnomAD_genome_ALL") gnomad_genome_all = pair[2];
          if (pair[1] == "gnomAD_genome_EAS") gnomad_genome_eas = pair[2];
          if (pair[1] == "ALL.sites.2015_08") all_sites_2015_08 = pair[2];
          if (pair[1] == "esp6500siv2_all") esp6500siv2_all = pair[2];
      }
      if ((exac_all == "." || exac_all <= 0.01) &&
          (exac_eas == "." || exac_eas <= 0.01) &&
          (exac_sas == "." || exac_sas <= 0.01) &&
          (gnomad_genome_all == "." || gnomad_genome_all <= 0.01) &&
          (gnomad_genome_eas == "." || gnomad_genome_eas <= 0.01) &&
          (all_sites_2015_08 == "." || all_sites_2015_08 <= 0.01) &&
          (esp6500siv2_all == "." || esp6500siv2_all <= 0.01))
          print
  }' > "$out"

  echo "  VCF pop-filtered: $vcf_file -> $out"
done

echo "== Pop-filter VCF complete."

########################
# 4) POP FREQ FILTER on ANNOVAR .txt (multianno)
########################
echo "== Pop-filter multianno TXT -> $POP_TXT_DIR"

find "$ANN_DIR" -type f -name "*.hg38_multianno.txt" | while read -r txt_file; do
  # e.g., annotated/<sample>/<sample>.hg38_multianno.txt -> <sample>_pop_filt.txt
  base=$(basename "$txt_file" .hg38_multianno.txt)
  out="${POP_TXT_DIR}/${base}.hg38_multianno_pop_filt.txt"

  awk 'BEGIN {FS=OFS="\t"}
  NR==1 {print; next}
  {
      # NOTE: column indexes hinge on your annovar version/db order; using your original mapping:
      exac_all=$12; exac_eas=$15; exac_sas=$19; gnomad_genome_all=$20; gnomad_genome_eas=$24; all_sites_2015_08=$28; esp6500siv2_all=$29;
      if ((exac_all == "." || exac_all <= 0.01) &&
          (exac_eas == "." || exac_eas <= 0.01) &&
          (exac_sas == "." || exac_sas <= 0.01) &&
          (gnomad_genome_all == "." || gnomad_genome_all <= 0.01) &&
          (gnomad_genome_eas == "." || gnomad_genome_eas <= 0.01) &&
          (all_sites_2015_08 == "." || all_sites_2015_08 <= 0.01) &&
          (esp6500siv2_all == "." || esp6500siv2_all <= 0.01))
          print
  }' "$txt_file" > "$out"

  echo "  TXT pop-filtered: $txt_file -> $out"
done

echo "== Pop-filter TXT complete."

########################
# 5) Build MAFs (from pop-filtered multianno .txt)
########################
echo "== Creating MAFs -> $MAF_DIR"

Rscript - <<'RS'
library(maftools)

input_dir <- "/home/users/nus/ash.ps/scratch/FINAL-sepfunc/batch2/pop-filter-multianno"
output_dir <- "/home/users/nus/ash.ps/scratch/FINAL-sepfunc/batch2/mafs"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

create_maf <- function(files, output_name){
  if (length(files) > 0) {
    maf_df <- annovarToMaf(
      files = files,
      Center = NULL,
      refBuild = "hg38",
      tsbCol = NULL,
      ens2hugo = TRUE,
      basename = NULL,
      sep = "\t",
      MAFobj = FALSE,
      sampleAnno = NULL
    )
    out <- file.path(output_dir, paste0(output_name, ".maf"))
    write.table(maf_df, file = out, sep = "\t", quote = FALSE, row.names = FALSE)
    cat("Created", output_name, "MAF with", length(files), "files\n")
  } else {
    cat("No files for", output_name, "\n")
  }
}

all_files <- list.files(path = input_dir, pattern = "\\.hg38_multianno_pop_filt\\.txt$", full.names = TRUE)

# One combined MAF across all samples
create_maf(all_files, "combined")

# (Optional) If you want custom groupings for batch2, adjust the regex below.
# Here we just create an example split by prefix B_* vs others:
create_maf(grep("/B_", all_files, value = TRUE), "B_prefix")

cat("MAF generation complete.\n")
RS

echo "== All done =="
