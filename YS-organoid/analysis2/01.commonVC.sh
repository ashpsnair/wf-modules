#!/bin/bash

#PBS -l select=2:ncpus=32:mem=64g
#PBS -l walltime=04:00:00
#PBS -P 11003581
#PBS -N commonVC_variants
#PBS -j oe
#PBS -m ae

# =============================================================================
# commonVC_variants.sh
#
# PURPOSE: Intersect Mutect2 and Strelka2 VCFs to identify variants called
#          by both callers. Includes full preprocessing: concat strelka
#          SNVs+INDELs, left-align / normalise indels, decompose
#          multi-allelic sites, dedup, and PASS-filter before isec.
#
# INPUT  : Mutect2 filtered VCFs   : $MUTECT2_DIR/*.mutect2.filtered.vcf.gz
#          Strelka2 SNV VCFs       : $STRELKA_SNV_DIR/*.strelka.somatic_snvs.vcf.gz
#          Strelka2 INDEL VCFs     : $STRELKA_INDEL_DIR/*.strelka.somatic_indels.vcf.gz
#
# OUTPUT : Per-sample intersection dirs (bcftools isec layout):
#            $INTERSECTION_DIR/<sample>/0000.vcf.gz  = Mutect2-only (PASS, norm)
#            $INTERSECTION_DIR/<sample>/0001.vcf.gz  = Strelka-only (PASS, norm)
#            $INTERSECTION_DIR/<sample>/0002.vcf.gz  = common (from Mutect2)
#            $INTERSECTION_DIR/<sample>/0003.vcf.gz  = common (from Strelka)
#          Strelka merged VCFs     : $STRELKA_MERGED_DIR/
#          Summary CSV             : $BASE_DIR/commonVC_summary.csv
#
# NOTES  : - bcftools norm -m -any decomposes multi-allelic records before
#             isec so that e.g. A>C,T is compared as A>C and A>T separately.
#           - bcftools norm -f <REF> left-aligns indels (critical for correct
#             matching between callers whose internal representations differ).
#           - bcftools norm -d exact removes exact duplicates after splitting.
#           - PASS filter applied to both callers before isec.
#           - isec uses -c all (match on REF+ALT after normalisation, ignore
#             INFO/FORMAT differences between callers).
# =============================================================================

set -euo pipefail
shopt -s nullglob

cd "$PBS_O_WORKDIR"

echo "=== commonVC_variants STARTED: $(date) ==="

module load bcftools/1.15.1

# ---------------------------------------------------------------------------
# PATHS — edit these for your project
# ---------------------------------------------------------------------------

REF="/home/project/11003581/Ref/sarek-refs/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"

BASE_DIR="/home/users/nus/ash.ps/scratch/YS-organoid-downstream/analysis2.1"

MUTECT2_DIR="${BASE_DIR}/raw-vcfs/mutect2"
STRELKA_SNV_DIR="${BASE_DIR}/raw-vcfs/strelka/snv"
STRELKA_INDEL_DIR="${BASE_DIR}/raw-vcfs/strelka/indel"

STRELKA_MERGED_DIR="${BASE_DIR}/00_strelka_merged"
INTERSECTION_DIR="${BASE_DIR}/01_commonVC"

TEMP_DIR="/tmp/commonVC_$$"

mkdir -p \
    "$STRELKA_MERGED_DIR" \
    "$INTERSECTION_DIR" \
    "$TEMP_DIR"

SUMMARY_CSV="${BASE_DIR}/commonVC_summary.csv"
echo "sample,mutect2_raw,strelka_raw,mutect2_pass,strelka_pass,mutect2_only,strelka_only,common" \
    > "$SUMMARY_CSV"

# ---------------------------------------------------------------------------
# SANITY CHECK
# ---------------------------------------------------------------------------

echo ""
echo "=== SANITY CHECK ==="
echo "--- Mutect2 VCFs ---"
find "$MUTECT2_DIR" -name "*.mutect2.filtered.vcf.gz" | sort
echo "--- Strelka SNV VCFs ---"
find "$STRELKA_SNV_DIR" -name "*.vcf.gz" | sort
echo "--- Strelka INDEL VCFs ---"
find "$STRELKA_INDEL_DIR" -name "*.vcf.gz" | sort

# ---------------------------------------------------------------------------
# STEP 1 — Concatenate Strelka SNVs + INDELs, sort, index
# ---------------------------------------------------------------------------

echo ""
echo "=== STEP 1: Concatenating Strelka SNVs + INDELs ==="

while read -r mutect2_file; do

    sample=$(basename "$mutect2_file" .mutect2.filtered.vcf.gz)

    snv_file="${STRELKA_SNV_DIR}/${sample}.strelka.somatic_snvs.vcf.gz"
    indel_file="${STRELKA_INDEL_DIR}/${sample}.strelka.somatic_indels.vcf.gz"
    merged_file="${STRELKA_MERGED_DIR}/${sample}.strelka.merged.vcf.gz"

    echo ""
    echo "--- $sample ---"

    if [[ ! -f "$snv_file" ]]; then
        echo "[ERROR] Missing SNV file: $snv_file"
        continue
    fi

    if [[ ! -f "$indel_file" ]]; then
        echo "[ERROR] Missing INDEL file: $indel_file"
        continue
    fi

    # Concat (allow overlaps between SNV and INDEL files with -a)
    bcftools concat \
        -a \
        -O z \
        -o "$TEMP_DIR/${sample}.concat.vcf.gz" \
        "$snv_file" \
        "$indel_file"

    # Sort by genomic coordinate
    bcftools sort \
        -O z \
        -o "$merged_file" \
        "$TEMP_DIR/${sample}.concat.vcf.gz"

    bcftools index -f --tbi "$merged_file"

    rm -f "$TEMP_DIR/${sample}.concat.vcf.gz"

    n=$(bcftools view -H "$merged_file" | wc -l)
    echo "[OK] Strelka merged: $n variants"

done < <(find "$MUTECT2_DIR" -name "*.mutect2.filtered.vcf.gz" | sort)

# ---------------------------------------------------------------------------
# STEP 2 — PASS filter, left-align + normalise, decompose multi-allelic,
#           dedup, then intersect
# ---------------------------------------------------------------------------

echo ""
echo "=== STEP 2: Normalise + Intersect ==="

while read -r mutect2_file; do

    sample=$(basename "$mutect2_file" .mutect2.filtered.vcf.gz)
    strelka_file="${STRELKA_MERGED_DIR}/${sample}.strelka.merged.vcf.gz"
    isec_dir="${INTERSECTION_DIR}/${sample}"

    echo ""
    echo "--- $sample ---"

    if [[ ! -f "$strelka_file" ]]; then
        echo "[ERROR] Missing merged Strelka: $strelka_file"
        continue
    fi

    mkdir -p "$isec_dir"

    # Raw counts before any filtering
    m2_raw=$(bcftools view -H "$mutect2_file" | wc -l)
    sk_raw=$(bcftools view -H "$strelka_file" | wc -l)

    # ---- PASS filter --------------------------------------------------------

    bcftools view -f PASS \
        -O z -o "$TEMP_DIR/${sample}.m2.pass.vcf.gz" \
        "$mutect2_file"
    bcftools index -f --tbi "$TEMP_DIR/${sample}.m2.pass.vcf.gz"

    bcftools view -f PASS \
        -O z -o "$TEMP_DIR/${sample}.sk.pass.vcf.gz" \
        "$strelka_file"
    bcftools index -f --tbi "$TEMP_DIR/${sample}.sk.pass.vcf.gz"

    m2_pass=$(bcftools view -H "$TEMP_DIR/${sample}.m2.pass.vcf.gz" | wc -l)
    sk_pass=$(bcftools view -H "$TEMP_DIR/${sample}.sk.pass.vcf.gz" | wc -l)

    # ---- Left-align + normalise indels (-m -any decomposes multi-allelic) ---
    # -m -any  : split multi-allelic into separate records
    # -f       : left-align indels using the reference genome
    # Running twice (split then left-align) is the safest approach

    bcftools norm \
        -m -any \
        -f "$REF" \
        -O z -o "$TEMP_DIR/${sample}.m2.norm.vcf.gz" \
        "$TEMP_DIR/${sample}.m2.pass.vcf.gz"
    bcftools index -f --tbi "$TEMP_DIR/${sample}.m2.norm.vcf.gz"

    bcftools norm \
        -m -any \
        -f "$REF" \
        -O z -o "$TEMP_DIR/${sample}.sk.norm.vcf.gz" \
        "$TEMP_DIR/${sample}.sk.pass.vcf.gz"
    bcftools index -f --tbi "$TEMP_DIR/${sample}.sk.norm.vcf.gz"

    # ---- Remove exact duplicates (can arise after splitting) ----------------

    bcftools norm \
        -d exact \
        -O z -o "$TEMP_DIR/${sample}.m2.dedup.vcf.gz" \
        "$TEMP_DIR/${sample}.m2.norm.vcf.gz"
    bcftools index -f --tbi "$TEMP_DIR/${sample}.m2.dedup.vcf.gz"

    bcftools norm \
        -d exact \
        -O z -o "$TEMP_DIR/${sample}.sk.dedup.vcf.gz" \
        "$TEMP_DIR/${sample}.sk.norm.vcf.gz"
    bcftools index -f --tbi "$TEMP_DIR/${sample}.sk.dedup.vcf.gz"

    # ---- Intersect ----------------------------------------------------------
    # -c all : match variants by REF+ALT position regardless of INFO/FORMAT
    # Output files:
    #   0000.vcf.gz = private to Mutect2
    #   0001.vcf.gz = private to Strelka
    #   0002.vcf.gz = common (Mutect2 records)
    #   0003.vcf.gz = common (Strelka records)

    bcftools isec \
        -p "$isec_dir" \
        -c all \
        -O z \
        "$TEMP_DIR/${sample}.m2.dedup.vcf.gz" \
        "$TEMP_DIR/${sample}.sk.dedup.vcf.gz"

    for f in "$isec_dir"/*.vcf.gz; do
        bcftools index -f --tbi "$f"
    done

    # ---- Counts -------------------------------------------------------------

    m2_only=$(bcftools view -H "$isec_dir/0000.vcf.gz" | wc -l)
    sk_only=$(bcftools view -H "$isec_dir/0001.vcf.gz" | wc -l)
    common=$(bcftools view  -H "$isec_dir/0002.vcf.gz" | wc -l)

    echo "  Raw       — Mutect2: $m2_raw  Strelka: $sk_raw"
    echo "  PASS      — Mutect2: $m2_pass  Strelka: $sk_pass"
    echo "  Mutect2-only: $m2_only  Strelka-only: $sk_only  Common: $common"

    echo "$sample,$m2_raw,$sk_raw,$m2_pass,$sk_pass,$m2_only,$sk_only,$common" \
        >> "$SUMMARY_CSV"

    # ---- Clean temp ---------------------------------------------------------

    rm -f "$TEMP_DIR/${sample}."*

    echo "[OK] Done: $sample"

done < <(find "$MUTECT2_DIR" -name "*.mutect2.filtered.vcf.gz" | sort)

# ---------------------------------------------------------------------------
# CLEANUP
# ---------------------------------------------------------------------------

rm -rf "$TEMP_DIR"

echo ""
echo "=== commonVC_variants COMPLETE: $(date) ==="
echo ""
echo "Strelka merged  : $STRELKA_MERGED_DIR"
echo "Intersections   : $INTERSECTION_DIR"
echo "  (common variants are in <sample>/0002.vcf.gz)"
echo "Summary CSV     : $SUMMARY_CSV"
echo ""
cat "$SUMMARY_CSV"