#!/bin/bash
#PBS -l select=1:ncpus=64:mem=128g
#PBS -l walltime=04:00:00
#PBS -P 11003581
#PBS -N annovar_annotate_only
#PBS -j oe

set -euo pipefail
cd "$PBS_O_WORKDIR"

module load bcftools/1.15.1
module load r/4.2.0

########################
# PATHS
########################
BASE_DIR="/home/users/nus/ash.ps/scratch/FINAL-sepfunc/batch2"
RAW_DIR="${BASE_DIR}/unfilt-vcfs"     # has *.concordant.vcf.gz
ANN_DIR="${BASE_DIR}/annotated"       # per-sample annovar outputs
MAF_DIR="${BASE_DIR}/mafs_per_sample" # per-sample maf outputs

mkdir -p "$ANN_DIR" "$MAF_DIR"

########################
# ANNOVAR CONFIG
########################
HUMANDB="/home/project/11003581/Tools/annovar/humandb"
ANNOVAR_PL="/home/project/11003581/Tools/annovar/table_annovar.pl"

# Keep your same protocols/operations
PROT="refGene,cosmic70,exac03,gnomad_genome,ALL.sites.2015_08,SAS.sites.2015_08,EAS.sites.2015_08,esp6500siv2_all,avsnp151,icgc28,dbnsfp30a,intervar_20180118,clinvar_20240611"
OPER="g,f,f,f,f,f,f,f,f,f,f,f,f"

########################
# 1) ANNOTATE (NO FILTERING)
########################
echo "== Running ANNOVAR (no filtering) =="
echo "Input:  $RAW_DIR"
echo "Output: $ANN_DIR"

find "$RAW_DIR" -type f -name "*.concordant.vcf.gz" | sort | while read -r vcf_gz; do
  samp=$(basename "$vcf_gz" .concordant.vcf.gz)
  outdir="${ANN_DIR}/${samp}"
  mkdir -p "$outdir"

  echo "  ANNOVAR: $samp"

  # Convert gz VCF -> plain VCF for table_annovar -vcfinput -polish
  tmp_vcf="${outdir}/${samp}.tmp.vcf"
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

echo "== ANNOVAR complete =="

########################
# 2) PER-SAMPLE MAFs (from ANNOVAR multianno TXT)
########################
echo "== Creating per-sample MAFs -> $MAF_DIR =="

Rscript - <<'RS'
library(maftools)

ann_dir <- "/home/users/nus/ash.ps/scratch/FINAL-sepfunc/batch2/annotated"
maf_dir <- "/home/users/nus/ash.ps/scratch/FINAL-sepfunc/batch2/mafs_per_sample"
dir.create(maf_dir, recursive = TRUE, showWarnings = FALSE)

# Grab all multianno txt files (one per sample produced by table_annovar)
files <- list.files(
  path = ann_dir,
  pattern = "\\.hg38_multianno\\.txt$",
  full.names = TRUE,
  recursive = TRUE
)

if (length(files) == 0) {
  stop("No *.hg38_multianno.txt found under: ", ann_dir)
}

for (f in files) {
  # Example input:
  # annotated/<sample>/<sample>.hg38_multianno.txt
  base <- sub("\\.hg38_multianno\\.txt$", "", basename(f))
  out  <- file.path(maf_dir, paste0(base, ".maf"))

  cat("Creating MAF for:", base, "\n")

  maf_df <- annovarToMaf(
    files = f,
    Center = NULL,
    refBuild = "hg38",
    tsbCol = NULL,
    ens2hugo = TRUE,
    basename = NULL,
    sep = "\t",
    MAFobj = FALSE,
    sampleAnno = NULL
  )

  write.table(maf_df, file = out, sep = "\t", quote = FALSE, row.names = FALSE)
}

cat("Done. Created", length(files), "MAF files in:", maf_dir, "\n")
RS

echo "== All done =="
