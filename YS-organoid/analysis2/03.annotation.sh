#!/bin/bash

#PBS -l select=1:ncpus=16:mem=32g
#PBS -l walltime=06:00:00
#PBS -P 11003581
#PBS -N annovar-annotation
#PBS -j oe
#PBS -m ae

# =============================================================================
# annovar-annotation.sh
#
# PURPOSE: Annotate base-filtered VCFs using ANNOVAR. Produces both the
#          multianno TXT file (used for downstream text-based filters) and
#          an annotated VCF (used for VCF-based filters).
#
# DATABASES USED:
#   refGene          — gene-level annotation (exonic, intronic, splicing, …)
#   cosmic70         — COSMIC cancer mutations
#   exac03           — ExAC population allele frequencies
#   gnomad_genome    — gnomAD whole-genome allele frequencies
#   avsnp151         — dbSNP identifiers (build 151)
#   dbnsfp30a        — functional effect prediction scores
#   clinvar_20240611 — ClinVar clinical significance
#
# INPUT  : $BASE_FILTER_DIR/<sample>.base_filtered.vcf.gz
#
# OUTPUT : $ANNOTATED_DIR/<sample>/
#            <sample>.hg38_multianno.txt   — tab-delimited annotation table
#            <sample>.hg38_multianno.vcf   — annotated VCF (plain text)
#            <sample>.hg38_multianno.vcf.gz — bgzipped + indexed VCF
#
# NOTES  : ANNOVAR requires uncompressed VCF input; this script decompresses
#          before each run and cleans up afterward.
#          Adjust ANNOVAR_DIR and HUMANDB_DIR to match your installation.
# =============================================================================

set -euo pipefail

cd "$PBS_O_WORKDIR"

echo "=== annovar-annotation STARTED: $(date) ==="

module load bcftools/1.15.1

# ---------------------------------------------------------------------------
# PATHS
# ---------------------------------------------------------------------------

BASE_DIR="/home/users/nus/ash.ps/scratch/YS-organoid-downstream/analysis2.1"

BASE_FILTER_DIR="${BASE_DIR}/02_base_filtered"
ANNOTATED_DIR="${BASE_DIR}/03_annotated"

ANNOVAR_DIR="/home/project/11003581/Tools/annovar"
HUMANDB_DIR="${ANNOVAR_DIR}/humandb"

BGZIP="/home/project/11003581/Tools/bin/bgzip"

mkdir -p "$ANNOTATED_DIR"

SUMMARY_CSV="${BASE_DIR}/annotation_summary.csv"
echo "sample,input_variants,annotated_variants" > "$SUMMARY_CSV"

# ---------------------------------------------------------------------------
# MAIN LOOP
# ---------------------------------------------------------------------------

for vcf_gz in "${BASE_FILTER_DIR}"/*.base_filtered.vcf.gz; do

    sample=$(basename "$vcf_gz" .base_filtered.vcf.gz)
    out_dir="${ANNOTATED_DIR}/${sample}"

    echo ""
    echo "=== Processing: $sample ==="

    mkdir -p "$out_dir"

    # -----------------------------------------------------------------------
    # Decompress to plain VCF (ANNOVAR requirement)
    # -----------------------------------------------------------------------

    plain_vcf="${out_dir}/${sample}.vcf"

    bcftools view "$vcf_gz" -O v -o "$plain_vcf"

    n_input=$(grep -vc '^#' "$plain_vcf" || true)
    echo "  Input variants: $n_input"

    # -----------------------------------------------------------------------
    # Run ANNOVAR table_annovar.pl
    #
    # Protocols used (matching Cell paper methodology):
    #   refGene          : gene function annotation
    #   cosmic70         : cancer somatic mutations
    #   exac03           : population AF (ExAC)
    #   gnomad_genome    : population AF (gnomAD)
    #   avsnp151         : dbSNP rsIDs
    #   dbnsfp30a        : functional prediction (SIFT, PolyPhen, CADD, …)
    #   clinvar_20240611 : clinical significance
    #
    # Operations: g=gene-based  f=filter-based
    # -----------------------------------------------------------------------

    perl "${ANNOVAR_DIR}/table_annovar.pl" \
        "$plain_vcf" \
        "$HUMANDB_DIR" \
        -buildver hg38 \
        -out "${out_dir}/${sample}" \
        -protocol refGene,cosmic70,exac03,gnomad_genome,avsnp151,dbnsfp30a,clinvar_20240611 \
        -operation  g,f,f,f,f,f,f \
        -remove \
        -vcfinput \
        -polish \
        -nastring .

    # -----------------------------------------------------------------------
    # Compress and index the annotated VCF output
    # -----------------------------------------------------------------------

    anno_vcf="${out_dir}/${sample}.hg38_multianno.vcf"
    anno_vcf_gz="${anno_vcf}.gz"

    if [[ -f "$anno_vcf" ]]; then

        "$BGZIP" -c "$anno_vcf" > "$anno_vcf_gz"
        bcftools index -f --tbi "$anno_vcf_gz"

        n_annotated=$(bcftools view -H "$anno_vcf_gz" | wc -l)
        echo "  Annotated variants: $n_annotated"

    else
        echo "[WARN] Annotated VCF not found: $anno_vcf"
        n_annotated=0
    fi

    # -----------------------------------------------------------------------
    # Clean up uncompressed input VCF (annotation table .txt is kept)
    # -----------------------------------------------------------------------

    rm -f "$plain_vcf"

    echo "$sample,$n_input,$n_annotated" >> "$SUMMARY_CSV"
    echo "[OK] Done: $sample"

done

# ---------------------------------------------------------------------------
# DONE
# ---------------------------------------------------------------------------

echo ""
echo "=== annovar-annotation COMPLETE: $(date) ==="
echo ""
echo "Annotated outputs : $ANNOTATED_DIR"
echo "Summary CSV       : $SUMMARY_CSV"
echo ""
cat "$SUMMARY_CSV"