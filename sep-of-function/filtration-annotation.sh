#!/bin/bash
#PBS -l select=1:ncpus=128:mem=256g
#PBS -l walltime=06:00:00
#PBS -P 11003581
#PBS -N batch2-filter-annotate
#PBS -j oe

set -euo pipefail
cd "$PBS_O_WORKDIR"

# ===== Modules & tools =====
module load bcftools/1.15.1

BGZIP="/home/project/11003581/Tools/bin/bgzip"
[[ -x "$BGZIP" ]] || BGZIP=$(which bgzip)

# ===== Paths =====
BASE="/home/users/nus/ash.ps/scratch/FINAL-sepfunc/batch2"
RAW="${BASE}/unfilt-vcfs"                 # /cohort/*.concordant.vcf.gz
PASS_DIR="${BASE}/filtered-vcfs"          # step 1 output
CLEAN_DIR="${BASE}/clean-vcfs"            # step 2 input to ANNOVAR (plain VCF)
ANN_DIR="${BASE}/annotated"               # step 2 output (ANNOVAR)
POP_TXT_DIR="${BASE}/pop-filter-multianno" # step 3 output
LOG_DIR="${BASE}/logs"
TMP_DIR="${BASE}/tmp"
mkdir -p "$PASS_DIR" "$CLEAN_DIR" "$ANN_DIR" "$POP_TXT_DIR" "$LOG_DIR" "$TMP_DIR"

STAMP="$(date +%Y%m%d_%H%M%S)"
CSV_FILE="${LOG_DIR}/counts_${STAMP}.csv"
echo "sample,cohort,step,variants" > "$CSV_FILE"

# ===== ANNOVAR setup =====
ANNOVAR_DIR="/home/project/11003581/Tools/annovar"
HUMANDB="${ANNOVAR_DIR}/humandb"
ANNOVAR_PL="${ANNOVAR_DIR}/table_annovar.pl"
export PATH="${ANNOVAR_DIR}:$PATH"
chmod 755 "${ANNOVAR_DIR}"/*.pl 2>/dev/null || true
PROT="refGene,cosmic70,exac03,gnomad_genome,ALL.sites.2015_08,SAS.sites.2015_08,EAS.sites.2015_08,esp6500siv2_all,avsnp151,icgc28,dbnsfp30a,intervar_20180118,clinvar_20240611"
OPER="g,f,f,f,f,f,f,f,f,f,f,f,f"

count_vcf() { bcftools view -H "$1" | wc -l; }
count_txt() { awk 'NR>1{c++} END{print c+0}' "$1"; }

# =======================================================
# Step 1 — FILTER=PASS
# =======================================================
find "$RAW" -mindepth 1 -maxdepth 1 -type d | sort | while read -r COHORT_DIR; do
  cohort="$(basename "$COHORT_DIR")"
  out_dir="${PASS_DIR}/${cohort}"
  mkdir -p "$out_dir"

  find "$COHORT_DIR" -type f -name "*.concordant.vcf.gz" | sort | while read -r f; do
    sample="$(basename "$f" .concordant.vcf.gz)"
    out="${out_dir}/${sample}.PASS.vcf.gz"
    bcftools view -i 'FILTER="PASS"' -Oz -o "$out" "$f"
    bcftools index -f "$out" >/dev/null
    echo "${sample},${cohort},PASS,$(count_vcf "$out")" >> "$CSV_FILE"
  done
done

# =======================================================
# Step 2 — Make plain VCFs & run ANNOVAR
# =======================================================
find "$PASS_DIR" -type f -name "*.PASS.vcf.gz" | sort | while read -r vcf_gz; do
  cohort="$(basename "$(dirname "$vcf_gz")")"
  sample="$(basename "$vcf_gz" .PASS.vcf.gz)"

  out_plain_dir="${CLEAN_DIR}/${cohort}"
  mkdir -p "$out_plain_dir"
  clean_vcf="${out_plain_dir}/${sample}.min.vcf"

  bcftools view -G -Ov "$vcf_gz" > "$clean_vcf"

  outdir="${ANN_DIR}/${sample}"
  mkdir -p "$outdir"

  perl "$ANNOVAR_PL" "$clean_vcf" "$HUMANDB" \
      -buildver hg38 \
      -out "${outdir}/${sample}" \
      -protocol $PROT \
      -operation $OPER \
      -remove -vcfinput -polish -nastring .

  anno_txt="${outdir}/${sample}.hg38_multianno.txt"
  n=$(awk 'NR>1{c++} END{print c+0}' "$anno_txt")
  echo "${sample},${cohort},ANNOVAR,${n}" >> "$CSV_FILE"
done

# =======================================================
# Step 3 — Population-frequency filter
# =======================================================
find "$ANN_DIR" -type f -name "*.hg38_multianno.txt" | sort | while read -r txt_file; do
  sample="$(basename "$txt_file" .hg38_multianno.txt)"
  out="${POP_TXT_DIR}/${sample}.hg38_multianno_pop_filt.txt"

  awk 'BEGIN{FS=OFS="\t"}
    NR==1{print; next}
    {
      exac_all=$12; exac_eas=$15; exac_sas=$19; gnomad_genome_all=$20; gnomad_genome_eas=$24; all_sites_2015_08=$28; esp6500siv2_all=$29;
      if ((exac_all=="." || exac_all<=0.01) &&
          (exac_eas=="." || exac_eas<=0.01) &&
          (exac_sas=="." || exac_sas<=0.01) &&
          (gnomad_genome_all=="." || gnomad_genome_all<=0.01) &&
          (gnomad_genome_eas=="." || gnomad_genome_eas<=0.01) &&
          (all_sites_2015_08=="." || all_sites_2015_08<=0.01) &&
          (esp6500siv2_all=="." || esp6500siv2_all<=0.01)) print
    }' "$txt_file" > "$out"

  cohort=$(grep -m1 "^${sample}," "$CSV_FILE" | cut -d, -f2)
  echo "${sample},${cohort},POP_TXT,$(count_txt "$out")" >> "$CSV_FILE"
done

echo "✅ Done. Counts CSV: $CSV_FILE"
