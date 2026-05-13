#!/bin/bash

#PBS -l select=1:ncpus=64:mem=128g
#PBS -l walltime=02:00:00
#PBS -P 11003581
#PBS -N annovar_maf_pipeline
#PBS -j oe

set -euo pipefail

cd $PBS_O_WORKDIR

module load bcftools/1.15.1

base_dir="/home/users/nus/ash.ps/scratch/YS-organoid-downstream"

############################################################
# STEP 1: ANNOVAR ANNOTATION
############################################################

input_dir="${base_dir}/final-common"
output_dir="${base_dir}/annotated"

mkdir -p "$output_dir"

echo "=== ANNOTATION STEP ==="

for vcf_gz in "$input_dir"/*_common.vcf.gz; do
    
    samplename=$(basename "$vcf_gz" _common.vcf.gz)
    
    echo "Processing: $samplename"
    
    mkdir -p "$output_dir/$samplename"
    
    # Convert to plain VCF (ANNOVAR requirement)
    vcf_file="$output_dir/$samplename/${samplename}.vcf"
    bcftools view "$vcf_gz" -Ov -o "$vcf_file"
    
    perl /home/project/11003581/Tools/annovar/table_annovar.pl "$vcf_file" \
        /home/project/11003581/Tools/annovar/humandb/ \
        -buildver hg38 \
        -out "$output_dir/$samplename/$samplename" \
        -protocol refGene,cosmic70,exac03,gnomad_genome,avsnp151,dbnsfp30a,clinvar_20240611 \
        -operation g,f,f,f,f,f,f \
        -remove -vcfinput -polish -nastring .
done

############################################################
# STEP 2A: POPULATION FILTER (VCF)
############################################################

INPUT_DIR="${base_dir}/annotated"
OUTPUT_DIR="${base_dir}/pop-filter-vcfs"

mkdir -p "$OUTPUT_DIR"

echo "=== POPULATION FILTER (VCF) ==="

find "$INPUT_DIR" -name "*.hg38_multianno.vcf" | while read -r vcf_file; do
    
    sample_name=$(basename "$vcf_file" .hg38_multianno.vcf)
    output_file="$OUTPUT_DIR/${sample_name}_pop_filt.vcf"
    
    bcftools view "$vcf_file" | awk '
    BEGIN {FS="\t"; OFS="\t"}
    /^#/ {print; next}
    {
        split($8, info, ";");
        exac="."; gnomad=".";
        
        for (i in info) {
            split(info[i], pair, "=");
            if (pair[1] == "ExAC_ALL") exac = pair[2];
            if (pair[1] == "gnomAD_genome_ALL") gnomad = pair[2];
        }
        
        if ((exac == "." || exac <= 0.01) &&
            (gnomad == "." || gnomad <= 0.01))
            print
    }' > "$output_file"
    
    echo "VCF filtered: $sample_name"
done


############################################################
# STEP 2 B: POPULATION FILTER (MULTIANNO)
############################################################

INPUT_DIR="${base_dir}/annotated"
OUTPUT_DIR="${base_dir}/pop-filter-multianno"

mkdir -p "$OUTPUT_DIR"

echo "=== POPULATION FILTER ==="

find "$INPUT_DIR" -name "*.hg38_multianno.txt" | while read -r txt_file; do
    
    output_file="$OUTPUT_DIR/$(basename "${txt_file%.txt}_pop_filt.txt")"
    
    awk 'BEGIN {FS=OFS="\t"}
    NR==1 {print; next}
    {
        exac_all=$12; gnomad_all=$20;
        if ((exac_all == "." || exac_all <= 0.01) &&
            (gnomad_all == "." || gnomad_all <= 0.01))
            print
    }' "$txt_file" > "$output_file"
    
    echo "Filtered: $txt_file"
done


############################################################
# STEP 3: CREATE MAF FILES
############################################################

module load r/4.2.0

echo "=== MAF GENERATION ==="

Rscript - << 'EOF'

library(maftools)

input_dir <- "/home/users/nus/ash.ps/scratch/YS-organoid-downstream/pop-filter-multianno"
output_dir <- "/home/users/nus/ash.ps/scratch/YS-organoid-downstream/mafs"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(input_dir, pattern = "\\.txt$", full.names = TRUE)

############################################
# PER-SAMPLE MAF
############################################

cat("Creating per-sample MAFs...\n")

maf_list <- list()

for (file in files) {
  
  sample_name <- gsub("_common.hg38_multianno_pop_filt.txt", "", basename(file))
  
  cat("Processing:", sample_name, "\n")
  
  maf <- annovarToMaf(
    file,
    refBuild = "hg38",
    ens2hugo = TRUE,
    MAFobj = FALSE
  )
  
  maf$Tumor_Sample_Barcode <- sample_name
  
  write.table(
    maf,
    file = file.path(output_dir, paste0(sample_name, ".maf")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  maf_list[[sample_name]] <- maf
}

############################################
# COMBINED MAF
############################################

cat("Creating combined MAF...\n")

combined_maf <- do.call(rbind, maf_list)

write.table(
  combined_maf,
  file = file.path(output_dir, "combined.maf"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("MAF generation complete\n")

EOF

echo "=== PIPELINE COMPLETE ==="