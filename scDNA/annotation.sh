#!/bin/bash

#PBS -l select=3:ncpus=64:mem=256g
#PBS -l walltime=06:00:00
#PBS -P 11003581
#PBS -N scdna-annotate
#PBS -j oe
#PBS -M ash.ps@nus.edu.sg

cd $PBS_O_WORKDIR
module load bcftools/1.15.1

base_dir="/home/users/nus/ash.ps/scratch/scDNA/analysis"
annot_dir="${base_dir}/annotation"

mkdir -p "${annot_dir}"/{annotated,pop-filter-vcfs,pop-filter-multianno,mafs}

###############################################
### annotation
###############################################

for cat in 2204 3401; do
  input_dir="${base_dir}/filtered-vcfs/$cat"
  output_dir="${annot_dir}/annotated/$cat"
  mkdir -p "$output_dir"

  for vcf_file in "$input_dir"/*.vcf.gz; do
    samplename=$(basename "$vcf_file" .vcf.gz)
    echo "Annotating: $samplename"
    mkdir -p "$output_dir/$samplename"

    perl /home/project/11003581/Tools/annovar/table_annovar.pl "$vcf_file" \
      /home/project/11003581/Tools/annovar/humandb/ \
      -buildver hg38 \
      -out "$output_dir/$samplename/$samplename" \
      -protocol refGene,cosmic70,exac03,gnomad_genome,ALL.sites.2015_08,SAS.sites.2015_08,EAS.sites.2015_08,esp6500siv2_all,avsnp151,icgc28,dbnsfp30a,intervar_20180118,clinvar_20240611 \
      -operation g,f,f,f,f,f,f,f,f,f,f,f,f \
      -remove -vcfinput -polish -nastring .
  done
done

###############################################
### pop filtration vcf files
###############################################

for cat in 2204 3401; do
  INPUT_DIR="${annot_dir}/annotated/$cat"
  OUTPUT_DIR="${annot_dir}/pop-filter-vcfs/$cat"
  mkdir -p "$OUTPUT_DIR"

  find "$INPUT_DIR" -name "*.hg38_multianno.vcf" | while read -r vcf_file; do
    output_file="$OUTPUT_DIR/$(basename "${vcf_file%.vcf}_pop_filt.vcf")"
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
      }' > "$output_file"
    echo "Processed: $vcf_file → $output_file"
  done
done

###############################################
### pop filtration on multianno files
###############################################

for cat in 2204 3401; do
  INPUT_DIR="${annot_dir}/annotated/$cat"
  OUTPUT_DIR="${annot_dir}/pop-filter-multianno/$cat"
  mkdir -p "$OUTPUT_DIR"

  find "$INPUT_DIR" -name "*.hg38_multianno.txt" | while read -r txt_file; do
    output_file="$OUTPUT_DIR/$(basename "${txt_file%.txt}_pop_filt.txt")"
    awk 'BEGIN {FS=OFS="\t"}
      NR==1 {print; next}
      {
        exac_all=$12; exac_eas=$15; exac_sas=$19; gnomad_genome_all=$20; gnomad_genome_eas=$24; all_sites_2015_08=$28; esp6500siv2_all=$29;
        if ((exac_all == "." || exac_all <= 0.01) &&
            (exac_eas == "." || exac_eas <= 0.01) &&
            (exac_sas == "." || exac_sas <= 0.01) &&
            (gnomad_genome_all == "." || gnomad_genome_all <= 0.01) &&
            (gnomad_genome_eas == "." || gnomad_genome_eas <= 0.01) &&
            (all_sites_2015_08 == "." || all_sites_2015_08 <= 0.01) &&
            (esp6500siv2_all == "." || esp6500siv2_all <= 0.01))
          print
      }' "$txt_file" > "$output_file"
    echo "Processed: $txt_file → $output_file"
  done
done

###############################################
### create maf files
###############################################

module load r/4.2.0

Rscript -e "
library(maftools)

base <- '/home/users/nus/ash.ps/scratch/scDNA/analysis/annotation'
input_dir <- file.path(base, 'pop-filter-multianno')
output_dir <- file.path(base, 'mafs')
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

create_maf <- function(files, output_name) {
  if (length(files) > 0) {
    maf <- annovarToMaf(
      files,
      refBuild = 'hg38',
      sep = '\t'
    )
    write.table(maf, file = file.path(output_dir, paste0(output_name, '.maf')), sep = '\t', quote = FALSE, row.names = FALSE)
    cat('Created', output_name, 'MAF with', length(files), 'files\\n')
  }
}

# all 2204 and 3401 files
all_files <- list.files(path = input_dir, pattern = '\\\\_pop_filt\\\\.txt$', recursive = TRUE, full.names = TRUE)
create_maf(all_files, 'combined_maf')

cat('MAF generation complete.\\n')
"




####################### BATCWISE- example ########################
#!/bin/bash

#PBS -N anno-3401-G-H
#PBS -l select=3:ncpus=64:mem=256g
#PBS -l walltime=22:00:00
#PBS -P 11003581
#PBS -j oe
#PBS -M ash.ps@nus.edu.sg

module load bcftools/1.15.1
module load r/4.2.0

base_dir="/home/users/nus/ash.ps/scratch/scDNA/analysis"
annot_dir="${base_dir}/annotation"
input_dir="${base_dir}/filtered-vcfs/3401"
anno_dir="${annot_dir}/annotated/3401"
popvcf_dir="${annot_dir}/pop-filter-vcfs/3401"
poptxt_dir="${annot_dir}/pop-filter-multianno/3401"
maf_dir="${annot_dir}/mafs"

mkdir -p "$anno_dir" "$popvcf_dir" "$poptxt_dir" "$maf_dir"

# Annotation
for vcf_file in "$input_dir"/somatic-3401-[G-H]*_filtered.vcf.gz; do
  [[ -e "$vcf_file" ]] || continue
  sample=$(basename "$vcf_file" .vcf.gz)
  mkdir -p "$anno_dir/$sample"
  perl /home/project/11003581/Tools/annovar/table_annovar.pl "$vcf_file" \
    /home/project/11003581/Tools/annovar/humandb/ \
    -buildver hg38 \
    -out "$anno_dir/$sample/$sample" \
    -protocol refGene,cosmic70,exac03,gnomad_genome,ALL.sites.2015_08,SAS.sites.2015_08,EAS.sites.2015_08,esp6500siv2_all,avsnp151,icgc28,dbnsfp30a,intervar_20180118,clinvar_20240611 \
    -operation g,f,f,f,f,f,f,f,f,f,f,f,f \
    -remove -vcfinput -polish -nastring .
done

# VCF Pop-filter
find "$anno_dir" -name "*.hg38_multianno.vcf" | while read -r vcf_file; do
  sample=$(basename "${vcf_file%.hg38_multianno.vcf}")
  output_file="$popvcf_dir/${sample}_pop_filt.vcf"
  bcftools view "$vcf_file" | awk '
    BEGIN {FS=OFS="\t"}
    /^#/ {print; next}
    {
      split($8, info, ";");
      exac_all="."; exac_eas="."; exac_sas="."; gnomad_all="."; gnomad_eas="."; all_sites="."; esp=".";
      for (i in info) {
        split(info[i], kv, "=");
        if (kv[1]=="ExAC_ALL") exac_all=kv[2];
        if (kv[1]=="ExAC_EAS") exac_eas=kv[2];
        if (kv[1]=="ExAC_SAS") exac_sas=kv[2];
        if (kv[1]=="gnomAD_genome_ALL") gnomad_all=kv[2];
        if (kv[1]=="gnomAD_genome_EAS") gnomad_eas=kv[2];
        if (kv[1]=="ALL.sites.2015_08") all_sites=kv[2];
        if (kv[1]=="esp6500siv2_all") esp=kv[2];
      }
      if ((exac_all=="." || exac_all<=0.01) &&
          (exac_eas=="." || exac_eas<=0.01) &&
          (exac_sas=="." || exac_sas<=0.01) &&
          (gnomad_all=="." || gnomad_all<=0.01) &&
          (gnomad_eas=="." || gnomad_eas<=0.01) &&
          (all_sites=="." || all_sites<=0.01) &&
          (esp=="." || esp<=0.01))
        print
    }' > "$output_file"
done

# TXT Pop-filter
find "$anno_dir" -name "*.hg38_multianno.txt" | while read -r txt_file; do
  sample=$(basename "${txt_file%.hg38_multianno.txt}")
  output_file="$poptxt_dir/${sample}_pop_filt.txt"
  awk 'BEGIN {FS=OFS="\t"}
    NR==1 {print; next}
    {
      exac_all=$12; exac_eas=$15; exac_sas=$19; gnomad_all=$20; gnomad_eas=$24; all_sites=$28; esp=$29;
      if ((exac_all == "." || exac_all <= 0.01) &&
          (exac_eas == "." || exac_eas <= 0.01) &&
          (exac_sas == "." || exac_sas <= 0.01) &&
          (gnomad_all == "." || gnomad_all <= 0.01) &&
          (gnomad_eas == "." || gnomad_eas <= 0.01) &&
          (all_sites == "." || all_sites <= 0.01) &&
          (esp == "." || esp <= 0.01))
        print
    }' "$txt_file" > "$output_file"
done

# MAF generation
Rscript -e "
library(maftools)
files <- list.files('${poptxt_dir}', pattern = '_pop_filt\\\\.txt$', full.names = TRUE)
if (length(files) > 0) {
  maf <- annovarToMaf(files, refBuild = 'hg38', sep = '\\t')
  write.table(maf, file = file.path('${maf_dir}', '3401_G_H.maf'), sep = '\\t', quote = FALSE, row.names = FALSE)
  cat('Created 3401_G_H.maf with', length(files), 'samples\\n')
} else {
  cat('No files found for MAF generation in batch 3401-G-H\\n')
}
"
