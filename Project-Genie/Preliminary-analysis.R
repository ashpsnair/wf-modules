# ------------------ 0.  Setup ------------------
setwd("~/NUS Dropbox/Aiswarya PS/MyWork/Project-Genie/02-Final-Analysis")
# install.packages(c("tidyverse", "arrow", "biomaRt"))  # once only
library(tidyverse)   # dplyr, readr, tibble, ggplot2 …
library(arrow)       # fast Parquet I/O
library(biomaRt)     # gene ↔︎ coordinate lookup (GRCh38 by default)

gistic_file  <- "~/NUS Dropbox/Aiswarya PS/MyWork/Project-Genie/00-Data/raw_data/cancer_filt_data_CNA.txt"   # 🔺 your GISTIC matrix
coord_cache  <- "chr13_coords.parquet"     # cached coordinates

# ------------------ 1.  Load GISTIC matrix ------------------
gistic <- read_tsv(gistic_file, col_types = cols()) %>%      # auto-detect col types
  column_to_rownames("Hugo_Symbol") %>%                      # rows = genes
  .[order(rownames(.)), ]                                    # keep row order sorted

# ------------------ 2.  Which samples have BRCA2 = -2? ------------------
brca2_row      <- gistic["BRCA2", ]
homdel_samples <- colnames(gistic)[brca2_row == -2]

cat(length(homdel_samples), "samples with BRCA2 -2\n")

# ------------------ 3.  Get chr 13 gene coordinates (cached) ------------
if (file.exists(coord_cache)) {
  coords <- read_parquet(coord_cache) %>% as_tibble()
} else {
  mart   <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
  coords <- getBM(
    attributes = c("hgnc_symbol", "chromosome_name",
                   "start_position", "end_position"),
    filters     = "hgnc_symbol",
    values      = rownames(gistic),
    mart        = mart
  ) %>%
    filter(chromosome_name == "13",
           !is.na(start_position)) %>%
    distinct(hgnc_symbol, .keep_all = TRUE) %>%              # drop dup symbols
    arrange(start_position) %>%
    rename(Hugo_Symbol = hgnc_symbol,
           start       = start_position,
           end         = end_position)
  
  write_parquet(coords, coord_cache)
}

# ------------------ 4.  Merge coordinates + chr13 CNV calls --------------
chr13 <- coords %>%
  column_to_rownames("Hugo_Symbol") %>%
  cbind(., gistic[rownames(.), , drop = FALSE])

# ------------------ 5.  Walk outward from BRCA2 for each sample ---------
results <- vector("list", length(homdel_samples))
brca2_idx <- which(rownames(chr13) == "BRCA2")

for (i in seq_along(homdel_samples)) {
  sample <- homdel_samples[i]
  calls  <- chr13[, sample]
  is_del <- calls == -2
  
  left_idx  <- brca2_idx
  while (left_idx  > 1              && is_del[left_idx  - 1]) left_idx  <- left_idx  - 1
  right_idx <- brca2_idx
  while (right_idx < nrow(chr13)    && is_del[right_idx + 1]) right_idx <- right_idx + 1
  
  left_gene  <- rownames(chr13)[left_idx]
  right_gene <- rownames(chr13)[right_idx]
  left_pos   <- chr13[left_idx , "start"]
  right_pos  <- chr13[right_idx, "end"]
  
  results[[i]] <- tibble(
    sample     = sample,
    left_gene  = left_gene,
    right_gene = right_gene,
    gene_count = right_idx - left_idx + 1,
    span_bp    = right_pos - left_pos
  )
}

summary_tbl <- bind_rows(results) %>%
  arrange(desc(span_bp))

print(head(summary_tbl))

# ------------------ 6.  Quick histogram ------------------
ggplot(summary_tbl, aes(span_bp / 1e6)) +
  geom_histogram(bins = 20) +
  labs(x = "BRCA2 homozygous-deletion span (Mb)",
       y = "# samples",
       title = "Extent of chr13 -2 blocks containing BRCA2") +
  theme_minimal(base_size = 12)
