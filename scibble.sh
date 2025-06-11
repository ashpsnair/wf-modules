#!/bin/bash

germline_dir="/Users/ash/NUS Dropbox/Aiswarya PS/Prof/MASH/MASH-new/All-data/true-germline"
normal_dir="/Users/ash/NUS Dropbox/Aiswarya PS/Prof/MASH/MASH-new/All-data/normal-raw"
out_matrix="$germline_dir/germline_in_normal_matrix.tsv"
out_stats="$germline_dir/germline_in_normal_stats.tsv"

# Get list of normal samples
normal_vcfs=("$normal_dir"/*.vcf.gz)
normal_names=()
for vcf in "${normal_vcfs[@]}"; do
    normal_names+=("$(basename "$vcf" .vcf.gz)")
done

# Header for the matrix
echo -e "Variant\t${normal_names[*]}" | tr ' ' '\t' > "$out_matrix"

# For stats
echo -e "Germline_Variant\tNum_Normals_With_Variant" > "$out_stats"

# For each germline VCF
for germline_vcf in "$germline_dir"/*.vcf.gz; do
    # Extract variants with AF > 0.5 (assuming FORMAT/AF or INFO/AF, adjust as needed)
    # Try FORMAT/AF first, fallback to INFO/AF if needed
    bcftools view -i 'FORMAT/AF[0]>0.5' "$germline_vcf" 2>/dev/null | \
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' > /tmp/germline_sites.txt
    if [[ ! -s /tmp/germline_sites.txt ]]; then
        bcftools view -i 'INFO/AF>0.5' "$germline_vcf" | \
            bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' > /tmp/germline_sites.txt
    fi

    while read -r variant; do
        row="$variant"
        present_count=0
        for vcf in "${normal_vcfs[@]}"; do
            # Check if this variant is present in the normal sample
            if bcftools view -H "$vcf" | awk -v var="$variant" 'BEGIN{FS="\t"} {if($1"\t"$2"\t"$4"\t"$5==var) {exit 1}}'; then
                row="$row\t0"
            else
                row="$row\t1"
                present_count=$((present_count+1))
            fi
        done
        echo -e "$row" >> "$out_matrix"
        echo -e "$variant\t$present_count" >> "$out_stats"
    done < /tmp/germline_sites.txt
done

rm /tmp/germline_sites.txt

echo "Matrix saved to $out_matrix"
echo "Stats saved to $out_stats"
