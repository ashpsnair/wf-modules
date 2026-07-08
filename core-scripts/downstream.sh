#!/bin/bash

#PBS -l select=1:ncpus=64:mem=128g
#PBS -l walltime=04:00:00
#PBS -P 11003581
#PBS -N downstream-4HNE
#PBS -j oe
#PBS -M ash.ps@nus.edu.sg

cd $PBS_O_WORKDIR

# ================= MODULES / PATH =================
module load bcftools/1.15.1
module load r/4.2.0
export PATH="/home/project/11003581/Tools/bin:$PATH"   # bgzip, tabix
# ====================================================

set -euo pipefail

# ================= CONFIG — EDIT THESE =================
REF=/home/project/11003581/Ref/sarek-refs/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta
OUTDIR=/home/users/nus/ash.ps/scratch/YS-4HNE/downstream/4HNE
MUTECT_LIST=$OUTDIR/mutect2.txt
STRELKA_LIST=$OUTDIR/strelka.txt
TUMOR_NAME_FALLBACK=TUMOR
NORMAL_NAME_FALLBACK=NORMAL
THREADS=${NCPUS:-$(nproc)}
# =========================================================

mkdir -p "$OUTDIR"/{strelka_merged,mutect_norm,strelka_norm,isec,consensus_vcfs,qc,\
filters/LRK_filter,filters/hard_filter,logs}

RUN_LOG="$OUTDIR/logs/run_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$RUN_LOG") 2>&1
log()  { echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*"; }
step() { echo; echo "===== $* ====="; }

BINOM_R="$OUTDIR/filters/LRK_filter/.binom_threshold.R"

get_tumor_name() {
  local vcf="$1" name
  name=$(bcftools view -h "$vcf" | sed -n 's/^##tumor_sample=//p' | head -1)
  [[ -n "$name" ]] && echo "$name" || echo "$TUMOR_NAME_FALLBACK"
}
get_normal_name() {
  local vcf="$1" name
  name=$(bcftools view -h "$vcf" | sed -n 's/^##normal_sample=//p' | head -1)
  [[ -n "$name" ]] && echo "$name" || echo "$NORMAL_NAME_FALLBACK"
}
get_sample_index() {  # 0-based column index, for FORMAT/TAG[idx] indexing
  local vcf="$1" name="$2"
  bcftools query -l "$vcf" | awk -v n="$name" '$0==n{print NR-1; found=1} END{if(!found) print 0}'
}

# ------------------------------------------------------------------
# STEP 1: merge Strelka snv + indel per sample
# ------------------------------------------------------------------
step1_merge_strelka() {
  step "STEP 1: merge Strelka snv+indel"
  while IFS= read -r indel_vcf; do
    [[ -z "$indel_vcf" ]] && continue
    [[ "$indel_vcf" == *somatic_indels* ]] || continue   # skip snv lines if list has both
    snv_vcf=${indel_vcf/somatic_indels/somatic_snvs}
    sample=$(basename "$indel_vcf" .strelka.somatic_indels.vcf.gz)
    [[ -f "$indel_vcf" ]] || { log "[ERROR] missing $indel_vcf"; exit 1; }
    [[ -f "$snv_vcf"   ]] || { log "[ERROR] missing $snv_vcf"; exit 1; }

    out="$OUTDIR/strelka_merged/${sample}.strelka.merged.vcf.gz"
    log "[merge] $sample"
    bcftools concat -a --threads "$THREADS" "$snv_vcf" "$indel_vcf" -Oz -o "${out}.tmp"
    bcftools sort "${out}.tmp" -Oz -o "$out" 2> "$OUTDIR/logs/${sample}.sort.log"
    bcftools index -f -t --threads "$THREADS" "$out"
    rm -f "${out}.tmp"
  done < "$STRELKA_LIST"
  log "STEP 1 done."
}

# ------------------------------------------------------------------
# STEP 2: normalize both callers, PASS-filter, isec consensus
# ------------------------------------------------------------------
step2_consensus() {
  step "STEP 2: normalize + isec consensus"
  local total; total=$(wc -l < "$MUTECT_LIST")
  local i=0
  while IFS= read -r mutect_vcf; do
    [[ -z "$mutect_vcf" ]] && continue
    i=$((i+1))
    sample=$(basename "$mutect_vcf" .mutect2.filtered.vcf.gz)
    strelka_vcf="$OUTDIR/strelka_merged/${sample}.strelka.merged.vcf.gz"
    [[ -f "$strelka_vcf" ]] || { log "[ERROR] no merged Strelka VCF for $sample"; exit 1; }
    tumor=$(get_tumor_name "$mutect_vcf")
    log "[$i/$total] consensus: $sample (tumor: $tumor)"

    m_norm="$OUTDIR/mutect_norm/${sample}.mutect2.norm.pass.vcf.gz"
    s_norm="$OUTDIR/strelka_norm/${sample}.strelka.norm.pass.vcf.gz"

    bcftools norm -m -any -f "$REF" --threads "$THREADS" "$mutect_vcf" 2> "$OUTDIR/logs/${sample}.mnorm.log" \
      | bcftools view -f PASS --threads "$THREADS" -Oz -o "$m_norm"
    bcftools index -f -t --threads "$THREADS" "$m_norm"

    bcftools norm -m -any -f "$REF" --threads "$THREADS" "$strelka_vcf" 2> "$OUTDIR/logs/${sample}.snorm.log" \
      | bcftools view -f PASS --threads "$THREADS" -Oz -o "$s_norm"
    bcftools index -f -t --threads "$THREADS" "$s_norm"

    isec_dir="$OUTDIR/isec/${sample}"
    rm -rf "$isec_dir"; mkdir -p "$isec_dir"
    bcftools isec -p "$isec_dir" -c none "$m_norm" "$s_norm" 2> "$OUTDIR/logs/${sample}.isec.log"

    for f in 0000 0001 0002 0003; do
      [[ -f "$isec_dir/$f.vcf" ]] && bgzip -f -@ "$THREADS" "$isec_dir/$f.vcf" && bcftools index -f -t --threads "$THREADS" "$isec_dir/$f.vcf.gz"
    done

    if [[ -f "$isec_dir/0002.vcf.gz" ]]; then
      cp "$isec_dir/0002.vcf.gz" "$OUTDIR/consensus_vcfs/${sample}.consensus.vcf.gz"
    else
      log "  [warn] $sample: zero overlapping calls between callers — writing empty consensus VCF"
      bcftools view -h "$m_norm" | bgzip -@ "$THREADS" > "$OUTDIR/consensus_vcfs/${sample}.consensus.vcf.gz"
    fi
    bcftools index -f -t --threads "$THREADS" "$OUTDIR/consensus_vcfs/${sample}.consensus.vcf.gz"

    n_mutect_only=$(bcftools view -H "$isec_dir/0000.vcf.gz" 2>/dev/null | wc -l)
    n_strelka_only=$(bcftools view -H "$isec_dir/0001.vcf.gz" 2>/dev/null | wc -l)
    n_shared=$(bcftools view -H "$isec_dir/0002.vcf.gz" 2>/dev/null | wc -l)
    log "  -> mutect_only=$n_mutect_only strelka_only=$n_strelka_only shared=$n_shared"
    printf "sample\tmutect_only\tstrelka_only\tshared\n%s\t%d\t%d\t%d\n" \
      "$sample" "$n_mutect_only" "$n_strelka_only" "$n_shared" > "$OUTDIR/qc/${sample}.counts.tsv"
  done < "$MUTECT_LIST"
  log "STEP 2 done."
}

# ------------------------------------------------------------------
# STEP 3: LRK filter
#   1) remove sites present in >1 subclone (regardless of genotype)
#   2) keep only calls whose tumor alt-read count falls within the
#      left-tailed 95% CI of Binomial(DP, 0.5)
# Region files use 2-column chrom/pos format with a non-".bed"
# extension (bcftools treats ".bed" files as 0-based BED, which
# silently breaks single-position start==end regions).
# ------------------------------------------------------------------
write_binom_r() {
  cat > "$BINOM_R" << 'EOF'
args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]; outfile <- args[2]
d <- read.table(infile, sep = "\t", header = FALSE,
                 col.names = c("CHROM","POS","REF","ALT","AD","DP"),
                 stringsAsFactors = FALSE)
d$alt_count <- sapply(strsplit(d$AD, ","), function(x) as.numeric(x[length(x)]))
d$DP <- as.numeric(d$DP)
d$lower95 <- qbinom(0.05, size = d$DP, prob = 0.5)
keep <- d[!is.na(d$DP) & d$DP > 0 & d$alt_count >= d$lower95, c("CHROM","POS")]
write.table(keep, outfile, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
EOF
}

step3_LRK_filter() {
  step "STEP 3: LRK filter"
  write_binom_r
  local CONS_DIR="$OUTDIR/consensus_vcfs"
  local OUT="$OUTDIR/filters/LRK_filter"
  local TMP="$OUT/tmp"; mkdir -p "$TMP"

  : > "$TMP/all_sites.txt"
  local total_input=0
  for vcf in "$CONS_DIR"/*.consensus.vcf.gz; do
    n=$(bcftools view -H "$vcf" | wc -l)
    total_input=$((total_input + n))
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$vcf" >> "$TMP/all_sites.txt"
  done
  log "total records across all consensus VCFs (with duplicates): $total_input"

  sort "$TMP/all_sites.txt" | uniq -c | awk '$1>1{print $2"\t"$3}' > "$TMP/recurrent.regions"
  local n_unique n_recurrent
  n_unique=$(sort -u "$TMP/all_sites.txt" | wc -l)
  n_recurrent=$(wc -l < "$TMP/recurrent.regions")
  log "unique sites: $n_unique | flagged recurrent (>1 subclone): $n_recurrent"

  local total; total=$(ls "$CONS_DIR"/*.consensus.vcf.gz | wc -l)
  local i=0
  for vcf in "$CONS_DIR"/*.consensus.vcf.gz; do
    i=$((i+1))
    sample=$(basename "$vcf" .consensus.vcf.gz)
    n_before=$(bcftools view -H "$vcf" | wc -l)

    local step1="$TMP/${sample}.step1.vcf.gz"
    if [[ -s "$TMP/recurrent.regions" ]]; then
      bcftools view -T ^"$TMP/recurrent.regions" --threads "$THREADS" "$vcf" -Oz -o "$step1"
    else
      cp "$vcf" "$step1"
    fi
    bcftools index -f -t --threads "$THREADS" "$step1"
    n_after_recurrent=$(bcftools view -H "$step1" | wc -l)

    local tumor; tumor=$(get_tumor_name "$step1")
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD\t%DP]\n' -s "$tumor" "$step1" \
      > "$TMP/${sample}.ad_dp.tsv"
    Rscript "$BINOM_R" "$TMP/${sample}.ad_dp.tsv" "$TMP/${sample}.keep_sites.tsv"
    awk 'BEGIN{OFS="\t"}{print $1,$2}' "$TMP/${sample}.keep_sites.tsv" > "$TMP/${sample}.keep.regions"

    if [[ -s "$TMP/${sample}.keep.regions" ]]; then
      bcftools view -T "$TMP/${sample}.keep.regions" --threads "$THREADS" "$step1" -Oz -o "$OUT/${sample}.LRK.vcf.gz"
    else
      bcftools view -h "$step1" | bgzip -@ "$THREADS" > "$OUT/${sample}.LRK.vcf.gz"
    fi
    bcftools index -f -t --threads "$THREADS" "$OUT/${sample}.LRK.vcf.gz"
    n_final_vcf=$(bcftools view -H "$OUT/${sample}.LRK.vcf.gz" | wc -l)
    log "[$i/$total] LRK: $sample  before=$n_before  after_recurrent_removal=$n_after_recurrent  final=$n_final_vcf"
  done
  rm -rf "$TMP"
  log "STEP 3 done."
}

# ------------------------------------------------------------------
# STEP 4: hard filter, on LRK-filtered VCFs
#   - FILTER=PASS
#   - FORMAT/AF[tumor]  >= 0.2
#   - FORMAT/DP[tumor]  > 10  AND  FORMAT/DP[normal] > 10
#   - INFO/MMQ[1] (alt-allele median mapping quality) > 40
#   - alt allele supported on both strands (FORMAT/SB: alt-fwd>0, alt-rev>0)
#   - extra: remove any site present in ALL samples (defensive; step 3
#     already removes anything recurring in >1 subclone)
# NOT implemented: rank-sum-style "no bias in read position/base quality/
# mapping quality" — Mutect2 doesn't emit MQRankSum/ReadPosRankSum/
# BaseQRankSum by default; only the direct thresholds above are enforced.
# ------------------------------------------------------------------
step4_hard_filter() {
  step "STEP 4: hard filter"
  local LRK_DIR="$OUTDIR/filters/LRK_filter"
  local OUT="$OUTDIR/filters/hard_filter"

  local total; total=$(ls "$LRK_DIR"/*.LRK.vcf.gz | wc -l)
  local i=0
  for vcf in "$LRK_DIR"/*.LRK.vcf.gz; do
    i=$((i+1))
    sample=$(basename "$vcf" .LRK.vcf.gz)
    tumor=$(get_tumor_name "$vcf")
    normal=$(get_normal_name "$vcf")
    tumor_idx=$(get_sample_index "$vcf" "$tumor")
    normal_idx=$(get_sample_index "$vcf" "$normal")
    n_in=$(bcftools view -H "$vcf" | wc -l)
    log "[$i/$total] $sample (tumor=$tumor[$tumor_idx] normal=$normal[$normal_idx]) input=$n_in"

    thresh="$OUT/${sample}.thresh.vcf.gz"
    bcftools view -f PASS "$vcf" \
      | bcftools filter -e "FORMAT/AF[${tumor_idx}:0]<0.2 || FORMAT/DP[${tumor_idx}]<=10 || FORMAT/DP[${normal_idx}]<=10 || INFO/MMQ[1]<=40" \
      -Oz -o "$thresh" --threads "$THREADS"
    bcftools index -f -t --threads "$THREADS" "$thresh"
    n_thresh=$(bcftools view -H "$thresh" | wc -l)

    strand_tsv="$OUT/${sample}.strand.tsv"
    bcftools query -f '%CHROM\t%POS[\t%SB]\n' -s "$tumor" "$thresh" > "$strand_tsv" 2>/dev/null || : > "$strand_tsv"
    regions="$OUT/${sample}.strand_pass.regions"
    awk 'BEGIN{OFS="\t"} NF>=3 { split($3,s,","); if (s[3]+0>0 && s[4]+0>0) print $1,$2 }' \
        "$strand_tsv" > "$regions"

    if [[ -s "$regions" ]]; then
      bcftools view -T "$regions" --threads "$THREADS" "$thresh" -Oz -o "$OUT/${sample}.hardfilter.vcf.gz"
    else
      log "  [warn] FORMAT/SB not usable for $sample — strand-bias check skipped"
      cp "$thresh" "$OUT/${sample}.hardfilter.vcf.gz"
    fi
    bcftools index -f -t --threads "$THREADS" "$OUT/${sample}.hardfilter.vcf.gz"
    n_final=$(bcftools view -H "$OUT/${sample}.hardfilter.vcf.gz" | wc -l)
    log "  -> after thresholds=$n_thresh, after strand check=$n_final"

    rm -f "$thresh" "$thresh.tbi" "$strand_tsv" "$regions"
  done

  step "STEP 4b: remove sites shared across ALL samples"
  : > "$OUT/all_sites.txt"
  local n_samples; n_samples=$(ls "$OUT"/*.hardfilter.vcf.gz | wc -l)
  for vcf in "$OUT"/*.hardfilter.vcf.gz; do
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$vcf" >> "$OUT/all_sites.txt"
  done
  sort "$OUT/all_sites.txt" | uniq -c | awk -v n="$n_samples" '$1==n{print $2"\t"$3}' > "$OUT/shared_all.regions"
  log "sites present in all $n_samples samples: $(wc -l < "$OUT/shared_all.regions")"

  for vcf in "$OUT"/*.hardfilter.vcf.gz; do
    sample=$(basename "$vcf" .hardfilter.vcf.gz)
    final="$OUT/${sample}.final.vcf.gz"
    if [[ -s "$OUT/shared_all.regions" ]]; then
      bcftools view -T ^"$OUT/shared_all.regions" --threads "$THREADS" "$vcf" -Oz -o "$final"
    else
      cp "$vcf" "$final"
    fi
    bcftools index -f -t --threads "$THREADS" "$final"
    log "$sample final: $(bcftools view -H "$final" | wc -l) records"
  done
  rm -f "$OUT/all_sites.txt"
  log "STEP 4 done."
}

# ================= RUN =================
step1_merge_strelka
step2_consensus
step3_LRK_filter
step4_hard_filter
log "Pipeline complete. Final VCFs: $OUTDIR/filters/hard_filter/*.final.vcf.gz  (full log: $RUN_LOG)"
