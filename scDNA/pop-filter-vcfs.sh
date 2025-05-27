#!/bin/bash

#PBS -l select=3:ncpus=64:mem=256g
#PBS -l walltime=04:00:00
#PBS -P 11003581
#PBS -N scdna-pop-filter-vcfs
#PBS -j oe
#PBS -M ash.ps@nus.edu.sg
#PBS -m abe

cd $PBS_O_WORKDIR

module load bcftools/1.15.1

base_dir="/home/users/nus/ash.ps/scratch/scDNA/analysis/annotation"
INPUT_DIR="${base_dir}/annotated"
OUTPUT_DIR="${base_dir}/pop-filter-vcfs"
SUMMARY_FILE="${OUTPUT_DIR}/variant_counts.txt"

mkdir -p "$OUTPUT_DIR"
echo -e "File\tVariants_Before\tVariants_After" > "$SUMMARY_FILE"

find "$INPUT_DIR" -type f -name "*.hg38_multianno.vcf" | while read -r vcf_file; do
    base_name=$(basename "${vcf_file%.vcf}")
    output_file="$OUTPUT_DIR/${base_name}_pop_filt.vcf"

    echo "Processing: $vcf_file"

    # Count variants before filtering (excluding header lines)
    variants_before=$(grep -vc "^#" "$vcf_file")

    # Sanitize header and fix FILTER values inline
    awk '
      BEGIN {FS=OFS="\t"}
      # Fix FILTER header line
      /^##FILTER=<ID=No.variants<4/ {
        print "##FILTER=<ID=LowVarSupport,Description=\"Low variant support\">"
        next
      }
      # Fix INFO tag header lines with invalid characters (optional)
      /^##INFO=<ID=fathmm-MKL_coding_score/ {
        print "##INFO=<ID=fathmm_MKL_coding_score,Number=.,Type=String,Description=\"fixed\""
        next
      }
      /^##INFO=<ID=fathmm-MKL_coding_pred/ {
        print "##INFO=<ID=fathmm_MKL_coding_pred,Number=.,Type=String,Description=\"fixed\""
        next
      }
      /^##INFO=<ID=GERP\+\+_RS/ {
        print "##INFO=<ID=GERPpp_RS,Number=.,Type=String,Description=\"fixed\""
        next
      }
      # Pass all other header lines
      /^#/ {print; next}

      # Fix FILTER column values in records
      {
        split($7, filters, ",");
        for (i in filters) {
          if (filters[i] == "No.variants<4") filters[i] = "LowVarSupport";
        }
        $7 = filters[1];
        for (i = 2; i in filters; i++) $7 = $7 "," filters[i];

        print
      }
    ' "$vcf_file" | awk '
      BEGIN {FS=OFS="\t"}
      /^#/ {print; next}
      {
        split($8, info, ";");
        exac_all=""; exac_eas=""; exac_sas=""; gnomad_all=""; gnomad_eas=""; all_sites=""; esp="";
        for (i in info) {
          split(info[i], pair, "=");
          if (pair[1] == "ExAC_ALL") exac_all = pair[2];
          if (pair[1] == "ExAC_EAS") exac_eas = pair[2];
          if (pair[1] == "ExAC_SAS") exac_sas = pair[2];
          if (pair[1] == "gnomAD_genome_ALL") gnomad_all = pair[2];
          if (pair[1] == "gnomAD_genome_EAS") gnomad_eas = pair[2];
          if (pair[1] == "ALL.sites.2015_08") all_sites = pair[2];
          if (pair[1] == "esp6500siv2_all") esp = pair[2];
        }
        if ((exac_all != "." && exac_all <= 0.01) &&
            (exac_eas != "." && exac_eas <= 0.01) &&
            (exac_sas != "." && exac_sas <= 0.01) &&
            (gnomad_all != "." && gnomad_all <= 0.01) &&
            (gnomad_eas != "." && gnomad_eas <= 0.01) &&
            (all_sites != "." && all_sites <= 0.01) &&
            (esp != "." && esp <= 0.01))
          print
      }
    ' > "$output_file"

    # Count filtered variants
    variants_after=$(grep -vc "^#" "$output_file")

    # Write to summary
    echo -e "${base_name}\t${variants_before}\t${variants_after}" >> "$SUMMARY_FILE"
    echo "→ Finished: $output_file (Before: $variants_before | After: $variants_after)"
done