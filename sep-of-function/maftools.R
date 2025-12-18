library(maftools)

input_dir <- "~/Downloads/pop-filter-multianno/"
output_dir <- "~/Downloads/pop-filter-multianno/mafs"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

create_maf <- function(files, output_name){
  if (length(files) > 0) {
    maf_df <- annovarToMaf(
      annovar = files,
      Center = NULL,
      refBuild = "hg38",
      tsbCol = NULL,
      ens2hugo = TRUE,
      basename = NULL,
      sep = "\t",
      MAFobj = FALSE,
      sampleAnno = NULL
    )
    out <- file.path(output_dir, paste0(output_name, ".maf"))
    write.table(maf_df, file = out, sep = "\t", quote = FALSE, row.names = FALSE)
    cat("Created", output_name, "MAF with", length(files), "files\n")
  } else {
    cat("No files for", output_name, "\n")
  }
}

all_files <- list.files(path = input_dir, pattern = "\\.hg38_multianno_pop_filt\\.txt$", full.names = TRUE)

# One combined MAF across all samples
create_maf(all_files, "all-genotypes")

# Extract prefix (before first underscore)
get_prefix <- function(file) {
  name <- basename(file)
  prefix <- sub("_.*", "", name)  # everything before first "_"
  return(prefix)
}

# Group files by prefix
prefix_groups <- split(all_files, sapply(all_files, get_prefix))

# Generate MAF per prefix
for (prefix in names(prefix_groups)) {
  create_maf(prefix_groups[[prefix]], prefix)
}


######## Visualization ##########


