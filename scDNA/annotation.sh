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


