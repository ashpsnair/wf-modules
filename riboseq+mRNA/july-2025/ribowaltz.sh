library(riboWaltz)
library(ggplot2)

# ====== Paths ======
bamdir  <- "/home/users/nus/ash.ps/scratch/JQQ/analysis-july/mamba/alignment/star/sorted/ribo_only"
gtfpath <- "/home/users/nus/ash.ps/scratch/refs/Homo_sapiens.GRCh38.114.gtf"
outdir  <- "/home/users/nus/ash.ps/scratch/JQQ/analysis-july/downstream/ribowaltz"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ====== Annotation ======
annotation_dt <- create_annotation(gtfpath = gtfpath,
                                   dataSource = "Ensembl 114",
                                   organism = "Homo sapiens")

# ====== Load BAMs ======
reads_list <- bamtolist(bamfolder = bamdir, annotation = annotation_dt)

# 3. Filter read lengths (28–32 nt)
filtered_list <- length_filter(reads_list,
                               length_filter_mode = "custom",
                               length_range = 28:32)

# 4. Compute P-site offsets
psite_offset <- psite(filtered_list, flanking = 6, extremity = "auto")

# 5. Attach P-site information
reads_psite_list <- psite_info(filtered_list, psite_offset)

# ====== Group samples by condition ======
# Match your BAM filenames, e.g. Ribo_rep1_20d.transcriptome.sorted.bam
input_samples <- list(
  "Day0"  = grep("_0d\\.",  names(reads_list), value = TRUE),
  "Day5"  = grep("_5d\\.",  names(reads_list), value = TRUE),
  "Day20" = grep("_20d\\.", names(reads_list), value = TRUE)
)
print(input_samples)

# ====== 1. Read length distribution ======
rlength_res <- rlength_distr(reads_list,
                             sample = input_samples,
                             multisamples = "average",
                             plot_style = "dodge",
                             cl = 99,
                             colour = c("#333f50", "#39827c", "gray70"))

ggsave(file.path(outdir, "read_length_distribution.pdf"),
       plot = rlength_res[["plot"]],
       width = 10, height = 10)

# ====== 2. Read extremity localization heatmap (per replicate) ======
# Example: plot for first replicate of Day0
heatmap_res <- rends_heat(reads_list, annotation_dt,
                          sample = input_samples$Day0[1], # change for other reps
                          cl = 85,
                          utr5l = 25, cdsl = 40, utr3l = 25,
                          colour = "#333f50")

# Loop over all groups and all replicates
for (grp in names(input_samples)) {
  for (s in input_samples[[grp]]) {
    heatmap_res <- rends_heat(reads_list, annotation_dt,
                              sample = s,
                              cl = 85,
                              utr5l = 25, cdsl = 40, utr3l = 25,
                              colour = "#333f50")
    ggsave(file.path(outdir, paste0("read_extremity_", s, ".pdf")),
           plot = heatmap_res[[paste0("plot_", s)]],
           width = 10, height = 7)
  }
}



# ====== 3. P-sites per region (stacked) ======
psite_region_stack <- region_psite(reads_psite_list, annotation_dt,
                                   sample = input_samples,
                                   multisamples = "average",
                                   plot_style = "stack",
                                   cl = 85,
                                   colour = c("#333f50", "gray70", "#39827c"))
ggsave(file.path(outdir, "psites_per_region_stacked.pdf"),
       plot = psite_region_stack[["plot"]],
       width = 10, height = 7)

# ====== 4. P-sites per region (RNA-aware/dodge) ======
psite_region_dodge <- region_psite(reads_psite_list, annotation_dt,
                                   sample = input_samples,
                                   multisamples = "average",
                                   plot_style = "dodge",
                                   cl = 85,
                                   colour = c("#333f50", "gray70", "#39827c"))
ggsave(file.path(outdir, "psites_per_region_dodge.pdf"),
       plot = psite_region_dodge[["plot"]],
       width = 10, height = 7)

# ====== 5. Trinucleotide periodicity (stratified by read length) ======
frames_strat <- frame_psite_length(reads_psite_list, annotation_dt,
                                   sample = input_samples,
                                   multisamples = "average",
                                   plot_style = "facet",
                                   region = "all",
                                   cl = 85,
                                   colour = "#333f50")
ggsave(file.path(outdir, "frame_psite_length.pdf"),
       plot = frames_strat[["plot"]],
       width = 10, height = 7)

# ====== 6. Trinucleotide periodicity (CDS only, not stratified) ======
frames_cds <- frame_psite(reads_psite_list, annotation_dt,
                           sample = input_samples,
                           multisamples = "average",
                           plot_style = "facet",
                           region = "cds",
                           colour = c("#333f50", "#39827c"))
ggsave(file.path(outdir, "frame_psite_cds.pdf"),
       plot = frames_cds[["plot"]],
       width = 10, height = 7)

# ====== 7. Metagene profiles  ======
meta_res <- metaprofile_psite(reads_psite_list, annotation_dt,
                              sample = input_samples,
                              multisamples = "average",
                              plot_style = "overlap", # overlap lines
                              utr5l = 20, cdsl = 40, utr3l = 20,
                              colour = c("#333f50", "#39827c", "gray70"))

ggsave(file.path(outdir, "metagene_profile_overlap.pdf"),
       plot = meta_res[["plot"]],
       width = 15, height = 8)



cds_coverage_example <- cds_coverage(reads_psite_list, annotation_dt, whole_transcriptome=FALSE)

ggsave(file.path(outdir, "metagene_profile_overlap.pdf"),
       plot = meta_res[["plot"]],
       width = 15, height = 8)


################################################
################ RIboseqQC ################
################################################

#!/usr/bin/env Rscript

library(RiboSeQC)

# ====== 1. Define paths ======
gtf_file   <- "/home/users/nus/ash.ps/scratch/refs/Homo_sapiens.GRCh38.114.gtf"
twobit_file <- "/home/users/nus/ash.ps/scratch/refs/hg38.2bit"
annot_dir  <- "/home/users/nus/ash.ps/scratch/JQQ/analysis-july/downstream/riboseqc"
bam_dir    <- "/home/users/nus/ash.ps/scratch/JQQ/analysis-july/mamba/alignment/star/sorted/ribo_only"
report_file <- file.path(annot_dir, "riboseqc_report.html")

# ====== 2. Create directories ======
dir.create(dirname(twobit_file), recursive = TRUE, showWarnings = FALSE)
dir.create(annot_dir, recursive = TRUE, showWarnings = FALSE)

# ====== 3. Download UCSC hg38.2bit if not present ======
if (!file.exists(twobit_file)) {
  message("[INFO] Downloading UCSC hg38.2bit...")
  download.file(
    url = "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit",
    destfile = twobit_file,
    mode = "wb"
  )
  message("[INFO] Download complete.")
} else {
  message("[INFO] UCSC hg38.2bit already exists. Skipping download.")
}

# ====== 4. Prepare annotation if not already present ======
annot_file <- file.path(annot_dir, "Homo_sapiens.GRCh38.114_Rannot")
if (!file.exists(paste0(annot_file, ".RData")) &&
    !file.exists(paste0(annot_file, ".Rannot"))) {
  
  message("[INFO] Preparing RiboSeQC annotation...")
  prepare_annotation_files(
    annotation_directory = annot_dir,
    gtf_file = gtf_file,
    twobit_file = twobit_file,
    scientific_name = "Homo.sapiens",
    annotation_name = "GRCh38_Ensembl114",
    forge_BSgenome = TRUE,
    create_TxDb = TRUE
  )
  message("[INFO] Annotation prepared.")
} else {
  message("[INFO] Annotation already prepared. Skipping.")
}

# ====== 5. List BAM files ======
bam_files <- list.files(path = bam_dir, pattern = "\\.bam$", full.names = TRUE)
if (length(bam_files) == 0) {
  stop("[ERROR] No BAM files found in: ", bam_dir)
}
sample_names <- sub("\\.bam$", "", basename(bam_files))

# ====== 6. Run RiboSeQC ======
message("[INFO] Running RiboSeQC_analysis...")
RiboseQC_analysis(
  annotation_file = annot_file,
  bam_files       = bam_files,
  sample_names    = sample_names,
  dest_names      = file.path(annot_dir, sample_names),
  report_file     = report_file,
  write_tmp_files = FALSE,  # set TRUE if you want intermediate RData results
  fast_mode       = TRUE    # set FALSE for full analysis (slower)
)

message("[INFO] RiboSeQC analysis complete.")
message("[INFO] Report saved at: ", report_file)

