if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

install.packages(c(
  "data.table", "dplyr", "tibble", "stringr", "ggplot2",
  "pheatmap", "survival", "survminer", "readr", "purrr"
))

BiocManager::install(c(
  "edgeR", "limma", "GSVA", "GSEABase", "maftools"
))

# Optional but useful for pathway gene sets
install.packages("msigdbr")

# scarHRD
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("sztup/scarHRD")