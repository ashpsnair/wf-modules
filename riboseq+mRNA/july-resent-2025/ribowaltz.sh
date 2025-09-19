#!/bin/bash
#PBS -N ribowaltz_full_all
#PBS -P 11003581
#PBS -l select=1:ncpus=32:mem=64g
#PBS -l walltime=04:00:00
#PBS -j oe

module load r/4.2.0
export R_LIBS_USER=/home/project/11003581/R/library/4.2/

BAMDIR="/home/users/nus/ash.ps/scratch/JQQ/analysis-july/mamba/alignment/star/sorted/ribo_only"
GTFPATH="/home/users/nus/ash.ps/scratch/refs/Homo_sapiens.GRCh38.114.gtf"
OUTDIR="/home/users/nus/ash.ps/scratch/JQQ/analysis-july/downstream/ribowaltz"

mkdir -p "$OUTDIR"

echo "[$(date)] Starting riboWaltz pipeline for all samples..." | tee "$OUTDIR/job.log"

Rscript - <<'EOF'
library(riboWaltz)
library(ggplot2)

bamdir  <- "$BAMDIR"
gtfpath <- "$GTFPATH"
outdir  <- "$OUTDIR"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ===== Annotation =====
annotation_dt <- create_annotation(gtfpath = gtfpath,
                                   dataSource = "Ensembl 114",
                                   organism = "Homo sapiens")

# ===== Load BAMs =====
reads_list <- bamtolist(bamfolder = bamdir, annotation = annotation_dt)

# ===== Filter reads (28–32 nt) =====
filtered_list <- length_filter(reads_list,
                               length_filter_mode = "custom",
                               length_range = 28:32)

# ===== Compute P-site offsets =====
psite_offset <- psite(filtered_list, flanking = 6, extremity = "auto")

# ===== Attach P-site info =====
reads_psite_list <- psite_info(filtered_list, psite_offset)

# ===== Group by condition =====
input_samples <- list(
  "Day0"  = grep("_0d\\.",  names(reads_list), value = TRUE),
  "Day5"  = grep("_5d\\.",  names(reads_list), value = TRUE),
  "Day20" = grep("_20d\\.", names(reads_list), value = TRUE)
)
print(input_samples)

# ===== 1. Group-level Read length distribution =====
rlength_res <- rlength_distr(reads_list,
                             sample = input_samples,
                             multisamples = "average",
                             plot_style = "dodge",
                             cl = 99,
                             colour = c("#333f50", "#39827c", "gray70"))
ggsave(file.path(outdir, "read_length_distribution_grouped.pdf"),
       plot = rlength_res[["plot"]], width = 7, height = 5)

# ===== 1b. Individual sample Read length distribution =====
for (s in names(reads_list)) {
  rl <- rlength_distr(reads_list,
                      sample = list(s = s),
                      multisamples = "independent",
                      plot_style = "dodge",
                      cl = 99,
                      colour = "#333f50")
  ggsave(file.path(outdir, paste0("read_length_distribution_", s, ".pdf")),
         plot = rl[["plot"]], width = 7, height = 5)
}

# ===== 2. Heatmaps for each replicate =====
for (s in names(reads_list)) {
  hm <- rends_heat(reads_list, annotation_dt,
                   sample = s, cl = 85,
                   utr5l = 25, cdsl = 40, utr3l = 25,
                   colour = "#333f50")
  ggsave(file.path(outdir, paste0("read_extremity_", s, ".pdf")),
         plot = hm[[paste0("plot_", s)]], width = 7, height = 5)
}

# ===== 3. Group-level P-sites per region =====
psite_region_stack <- region_psite(reads_psite_list, annotation_dt,
                                   sample = input_samples,
                                   multisamples = "average",
                                   plot_style = "stack",
                                   cl = 85,
                                   colour = c("#333f50", "gray70", "#39827c"))
ggsave(file.path(outdir, "psites_per_region_stacked_grouped.pdf"),
       plot = psite_region_stack[["plot"]], width = 7, height = 5)

# ===== 3b. Individual sample P-sites per region =====
for (s in names(reads_list)) {
  pr <- region_psite(reads_psite_list, annotation_dt,
                     sample = list(s = s),
                     multisamples = "independent",
                     plot_style = "stack",
                     cl = 85,
                     colour = c("#333f50", "gray70", "#39827c"))
  ggsave(file.path(outdir, paste0("psites_per_region_", s, ".pdf")),
         plot = pr[["plot"]], width = 7, height = 5)
}

# ===== 4. Frame periodicity (grouped) =====
frames_cds <- frame_psite(reads_psite_list, annotation_dt,
                           sample = input_samples,
                           multisamples = "average",
                           plot_style = "facet",
                           region = "cds",
                           colour = c("#333f50", "#39827c"))
ggsave(file.path(outdir, "frame_psite_cds_grouped.pdf"),
       plot = frames_cds[["plot"]], width = 7, height = 5)

# ===== 4b. Frame periodicity (per sample) =====
for (s in names(reads_list)) {
  fr <- frame_psite(reads_psite_list, annotation_dt,
                    sample = list(s = s),
                    multisamples = "independent",
                    plot_style = "facet",
                    region = "cds",
                    colour = c("#333f50"))
  ggsave(file.path(outdir, paste0("frame_psite_cds_", s, ".pdf")),
         plot = fr[["plot"]], width = 7, height = 5)
}

# ===== 5. Metagene profiles (grouped) =====
meta_res <- metaprofile_psite(reads_psite_list, annotation_dt,
                              sample = input_samples,
                              multisamples = "average",
                              plot_style = "overlap",
                              utr5l = 20, cdsl = 40, utr3l = 20,
                              colour = c("#333f50", "#39827c", "gray70"))
ggsave(file.path(outdir, "metagene_profile_grouped.pdf"),
       plot = meta_res[["plot"]], width = 7, height = 5)

# ===== 5b. Metagene profiles (per sample) =====
for (s in names(reads_list)) {
  mp <- metaprofile_psite(reads_psite_list, annotation_dt,
                          sample = list(s = s),
                          multisamples = "independent",
                          plot_style = "overlap",
                          utr5l = 20, cdsl = 40, utr3l = 20,
                          colour = c("#333f50"))
  ggsave(file.path(outdir, paste0("metagene_profile_", s, ".pdf")),
         plot = mp[["plot"]], width = 7, height = 5)
}
EOF

echo "[$(date)] riboWaltz all-sample pipeline finished." | tee -a "$OUTDIR/job.log"
