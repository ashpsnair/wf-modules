#!/bin/bash

#PBS -l select=3:ncpus=64:mem=256g
#PBS -l walltime=04:00:00
#PBS -P 11003581
#PBS -N scdna-pop-filter-vcfs
#PBS -j oe
#PBS -M ash.ps@nus.edu.sg

cd $PBS_O_WORKDIR

#!/bin/bash

INPUT_DIR="/home/users/nus/ash.ps/scratch/scDNA/analysis/filtered-vcfs/2204"
OUTPUT_DIR="/home/users/nus/ash.ps/scratch/scDNA/analysis/annotation/pop-filter-vcfs/2204"

mkdir -p "$OUTPUT_DIR"

find "$INPUT_DIR" -name "*.vcf.gz" | while read -r vcf_file; do
  base_name=$(basename "${vcf_file%.vcf.gz}")
  output_file="$OUTPUT_DIR/${base_name}_pop_filt.vcf"

  zcat "$vcf_file" | awk '
    BEGIN {FS="\t"; OFS="\t"}
    /^##FILTER=<ID=No.variants<4/ {
      print "##FILTER=<ID=LowVarSupport,Description=\"Low variant support\">"
      next
    }
    /^#/ { print; next }

    {
      # Clean up malformed FILTER values
      split($7, filters, ",");
      for (i in filters) {
        if (filters[i] == "No.variants<4") filters[i] = "LowVarSupport";
      }
      $7 = filters[1];
      for (i = 2; i in filters; i++) $7 = $7 "," filters[i];

      # Parse INFO
      split($8, info, ";");
      for (i in info) {
        split(info[i], kv, "=");
        tag[kv[1]] = kv[2];
      }

      exac_all = tag["ExAC_ALL"]; if (exac_all == "") exac_all = ".";
      exac_eas = tag["ExAC_EAS"]; if (exac_eas == "") exac_eas = ".";
      exac_sas = tag["ExAC_SAS"]; if (exac_sas == "") exac_sas = ".";
      gnomad_all = tag["gnomAD_genome_ALL"]; if (gnomad_all == "") gnomad_all = ".";
      gnomad_eas = tag["gnomAD_genome_EAS"]; if (gnomad_eas == "") gnomad_eas = ".";
      all_sites = tag["ALL.sites.2015_08"]; if (all_sites == "") all_sites = ".";
      esp = tag["esp6500siv2_all"]; if (esp == "") esp = ".";

      if ((exac_all == "." || exac_all <= 0.01) &&
          (exac_eas == "." || exac_eas <= 0.01) &&
          (exac_sas == "." || exac_sas <= 0.01) &&
          (gnomad_all == "." || gnomad_all <= 0.01) &&
          (gnomad_eas == "." || gnomad_eas <= 0.01) &&
          (all_sites == "." || all_sites <= 0.01) &&
          (esp == "." || esp <= 0.01))
        print

      delete tag
    }
  ' > "$output_file"

  echo "Processed: $vcf_file → $output_file"
done
