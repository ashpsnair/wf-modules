library(data.table)
library(dplyr)
library(stringr)
library(purrr)
library(scarHRD)

dir.create("temp/scarhrd_inputs", recursive = TRUE, showWarnings = FALSE)
dir.create("results/hrd", recursive = TRUE, showWarnings = FALSE)

seg <- readRDS("results/ascat3_segments_matched.rds")

# ---------------------------------------
# helper to detect columns robustly
# ---------------------------------------
pick_col <- function(nms, pattern) {
  out <- nms[str_detect(tolower(nms), pattern)][1]
  ifelse(length(out) == 0 || is.na(out), NA, out)
}

nms <- colnames(seg)

sample_col <- pick_col(nms, "sample")
chr_col    <- pick_col(nms, "chrom")
start_col  <- pick_col(nms, "start")
end_col    <- pick_col(nms, "end")
major_col  <- pick_col(nms, "major|a_cn|major_copy")
minor_col  <- pick_col(nms, "minor|b_cn|minor_copy")
total_col  <- pick_col(nms, "total|total_cn|copy_number")
ploidy_col <- pick_col(nms, "ploidy")

cat("Detected columns:\n")
print(c(sample_col, chr_col, start_col, end_col, major_col, minor_col, total_col, ploidy_col))

if (is.na(sample_col) || is.na(chr_col) || is.na(start_col) || is.na(end_col)) {
  stop("Could not detect basic ASCAT columns.")
}
if (is.na(minor_col) && is.na(major_col)) {
  stop("Could not detect allele-specific CN columns.")
}

scar_df <- seg %>%
  transmute(
    SampleID = substr(.data[[sample_col]], 1, 12),
    Chromosome = as.character(.data[[chr_col]]),
    Start_position = as.integer(.data[[start_col]]),
    End_position   = as.integer(.data[[end_col]]),
    A_cn = if (!is.na(major_col)) as.numeric(.data[[major_col]]) else as.numeric(.data[[total_col]]) - as.numeric(.data[[minor_col]]),
    B_cn = if (!is.na(minor_col)) as.numeric(.data[[minor_col]]) else as.numeric(.data[[total_col]]) - as.numeric(.data[[major_col]]),
    total_cn = if (!is.na(total_col)) as.numeric(.data[[total_col]]) else A_cn + B_cn,
    ploidy = if (!is.na(ploidy_col)) as.numeric(.data[[ploidy_col]]) else NA_real_
  ) %>%
  filter(!is.na(Start_position), !is.na(End_position), !is.na(total_cn), !is.na(A_cn), !is.na(B_cn)) %>%
  mutate(
    Chromosome = ifelse(str_detect(Chromosome, "^chr"), Chromosome, paste0("chr", Chromosome))
  )

saveRDS(scar_df, "results/hrd/scarhrd_input_all.rds")

# ---------------------------------------
# run scarHRD sample by sample
# ---------------------------------------
samples <- unique(scar_df$SampleID)

run_one_sample <- function(sid) {
  x <- scar_df %>% filter(SampleID == sid)

  if (nrow(x) < 10) return(NULL)

  infile <- file.path("temp/scarhrd_inputs", paste0(sid, "_scar_input.tsv"))
  fwrite(x, infile, sep = "\t")

  res <- tryCatch({
    out <- scarHRD::scar_score(
      seg = infile,
      reference = "grch38",
      seqz = FALSE,
      chr.in.names = TRUE
    )

    out <- as.data.frame(out)
    out$patient_id <- sid
    out
  }, error = function(e) {
    message("scarHRD failed for ", sid, ": ", e$message)
    NULL
  })

  res
}

hrd_list <- lapply(samples, run_one_sample)
hrd_df <- bind_rows(hrd_list)

# standardize column names a bit
colnames(hrd_df) <- make.names(colnames(hrd_df))

# Expected columns usually include:
# HRD, Telomeric.AI, LST, HRD.sum
write.table(hrd_df, "results/hrd/hrd_scores.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

print(head(hrd_df))