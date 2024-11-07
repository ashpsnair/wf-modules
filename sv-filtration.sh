#!/bin/bash

#PBS -l select=1:ncpus=64:mem=256g
#PBS -l walltime=15:00:00
#PBS -P 11003581
#PBS -N trial-auto-filter
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load gcc
module load python/3.12.1-gcc11
module load bcftools/1.15.1

# Define file names
WT_VCF="/home/project/11003581/Data/Ash/P3L-lab-analysis/sniffles/sniffles_out/WT_MCF10A_sniffles.vcf"
P3L_VCF="/home/project/11003581/Data/Ash/P3L-lab-analysis/sniffles/sniffles_out/p3292l_sniffles.vcf"
OUTPUT_DIR="/home/project/11003581/Data/Ash/P3L-lab-analysis"

# Step 1: Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Step 2: Compress the VCF files using bcftools view
/app/apps/bcftools/1.15.1/bin/bcftools view -Oz -o "${WT_VCF}.gz" "$WT_VCF"
/app/apps/bcftools/1.15.1/bin/bcftools view -Oz -o "${P3L_VCF}.gz" "$P3L_VCF"


/app/apps/bcftools/1.15.1/bin/bcftools index "${WT_VCF}.gz"
/app/apps/bcftools/1.15.1/bin/bcftools index "${P3L_VCF}.gz"


# Step 3: Run bcftools isec to get unique and shared variants
/app/apps/bcftools/1.15.1/bin/bcftools isec -n=1 -p "$OUTPUT_DIR/autofilter-sv/" "${WT_VCF}.gz" "${P3L_VCF}.gz"

# Step 3: Rename the file containing unique variants from p3l
mv "$OUTPUT_DIR/autofilter-sv/0001.vcf" "$OUTPUT_DIR/autofilter-sv/p3l-uniq-filtered.vcf"

# Step 5: Filter the unique variants based on specified criteria
/app/apps/bcftools/1.15.1/bin/bcftools filter -i 'INFO/AF[0] > 0.2 && FILTER="PASS" && (INFO/SVLEN[0] > 1000 || INFO/SVTYPE="BND")' \
"$OUTPUT_DIR/autofilter-sv/p3l-uniq-filtered.vcf" -o "$OUTPUT_DIR/autofilter-sv/p3l-sv-filtered.vcf" -O v

# Step 6: Remove rows with IMPRECISE in INFO column
awk '$0 ~ /^#/ || !/IMPRECISE/' "$OUTPUT_DIR/autofilter-sv/p3l-sv-filtered.vcf" > "$OUTPUT_DIR/autofilter-sv/p3l-final-filtered.vcf"


# Step 4: Generate statistics about the variants
{
    echo "WT Variants: $(wc -l < "$WT_VCF")"
    echo "p3l Variants: $(wc -l < "$P3L_VCF")"
    echo "p3l Unique Variants: $(wc -l < "$OUTPUT_DIR/autofilter-sv/p3l-uniq-filtered.vcf")"
     echo "p3l Unique Variants after filtration: $(wc -l < "$OUTPUT_DIR/autofilter-sv/p3l-final-filtered.vcf")"
} > "$OUTPUT_DIR/autofilter-sv/variant_stats.txt"

########## Annotation ##########

/home/project/11003581/Tools/AnnotSV/bin/AnnotSV -SVinputFile $OUTPUT_DIR/autofilter-sv/p3l-final-filtered.vcf \
    -bcftools /app/apps/bcftools/1.15.1/bin/bcftools \
    -bedtools /app/apps/bedtools/2.30.0/bin/bedtools \
    -annotationMode full -outputDir $OUTPUT_DIR/annotsv
