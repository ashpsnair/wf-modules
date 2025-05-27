#!/bin/bash

#PBS -l select=3:ncpus=64:mem=256g
#PBS -l walltime=04:00:00
#PBS -P 11003581
#PBS -N scdna-pop-filter-vcfs
#PBS -j oe
#PBS -M ash.ps@nus.edu.sg
#PBS -m abe

cd "$PBS_O_WORKDIR"

module load bcftools/1.15.1

base_dir="/home/users/nus/ash.ps/scratch/scDNA/analysis"
annot_dir="${base_dir}/annotation"
stats_file="${base_dir}/stats/pop_filter_summary.tsv"

echo -e "Cohort\tSample\tBefore_Filter\tAfter_Filter" > "$stats_file"

for cat in 2204 3401; do
  INPUT_DIR="${annot_dir}/annotated/$cat"
  OUTPUT_DIR="${annot_dir}/pop-filter-multianno/$cat"
  mkdir -p "$OUTPUT_DIR"

  find "$INPUT_DIR" -name "*.hg38_multianno.txt" | while read -r txt_file; do
    sample_name=$(basename "$txt_file" .hg38_multianno.txt)
    output_file="$OUTPUT_DIR/${sample_name}_pop_filt.txt"

    echo "🔍 Processing: $sample_name"

    # Use awk to dynamically map column names to indices and apply filters
    awk -v OFS="\t" '
      BEGIN {
        before = 0; after = 0;
      }
      NR == 1 {
        for (i = 1; i <= NF; i++) {
          if ($i == "ExAC_ALL") exac_all = i;
          if ($i == "ExAC_EAS") exac_eas = i;
          if ($i == "ExAC_SAS") exac_sas = i;
          if ($i == "gnomAD_genome_ALL") gnomad_all = i;
          if ($i == "gnomAD_genome_EAS") gnomad_eas = i;
          if ($i == "ALL.sites.2015_08") all_2015 = i;
          if ($i == "esp6500siv2_all") esp6500 = i;
        }
        print > "'"$output_file"'";
        next;
      }
      {
        before++;
        # Convert "." to large dummy value (9999) to bypass numeric comparison
        a = ($exac_all == "." ? 0 : $exac_all);
        b = ($exac_eas == "." ? 0 : $exac_eas);
        c = ($exac_sas == "." ? 0 : $exac_sas);
        d = ($gnomad_all == "." ? 0 : $gnomad_all);
        e = ($gnomad_eas == "." ? 0 : $gnomad_eas);
        f = ($all_2015 == "." ? 0 : $all_2015);
        g = ($esp6500 == "." ? 0 : $esp6500);

        if (a <= 0.01 && b <= 0.01 && c <= 0.01 &&
            d <= 0.01 && e <= 0.01 && f <= 0.01 && g <= 0.01) {
          print > "'"$output_file"'";
          after++;
        }
      }
      END {
        printf("📊 %s → Before: %d | After: %d\n", "'"$sample_name"'", before, after) > "/dev/stderr";
        printf("%s\t%s\t%d\t%d\n", "'"$cat"'", "'"$sample_name"'", before, after) >> "'"$stats_file"'";
      }
    ' "$txt_file"

  done
done

echo -e "\n✅ All done. Summary written to: $stats_file"
