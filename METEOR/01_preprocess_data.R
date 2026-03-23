############################################################
# PCAWG PREPROCESSING PIPELINE
# Shifu version 🔥 (robust + safe)
############################################################

library(tidyverse)
library(readxl)
library(stringr)
library(data.table)
library(ggplot2)

dir.create("output/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("output/plots", recursive = TRUE, showWarnings = FALSE)

############################################################
# 1. LOAD SIGNATURES
############################################################

sig <- read.delim("~/NUS Dropbox/Aiswarya PS/MyWork/METEOR/data/Assignment_Solution_Activities.txt")

# 🔥 ROBUST ID EXTRACTION
sig <- sig %>%
  mutate(
    sample_id = str_extract(Samples, "SP[0-9]+"),
    cancer_type = str_extract(Samples, "^[^:]+")
  ) %>%
  filter(!is.na(sample_id))

############################################################
# REMOVE DUPLICATES (VERY IMPORTANT)
############################################################

sig <- sig %>%
  distinct(sample_id, .keep_all = TRUE)

############################################################
# DEBUG
############################################################

cat("Example signature IDs:\n")
print(head(sig$sample_id))

############################################################
# 2. LOAD EXPRESSION
############################################################

expr <- fread("~/NUS Dropbox/Aiswarya PS/MyWork/METEOR/data/tophat_star_fpkm_uq.txt", data.table = FALSE)

# set gene names
rownames(expr) <- expr$feature
expr <- expr[, -1]

# transpose → samples x genes
expr_t <- as.data.frame(t(expr))

expr_t$sample_id <- rownames(expr_t)

############################################################
# CLEAN EXPRESSION IDS
############################################################

expr_t <- expr_t %>%
  mutate(sample_id = str_extract(sample_id, "SP[0-9]+")) %>%
  filter(!is.na(sample_id))

############################################################
# REMOVE DUPLICATES
############################################################

expr_t <- expr_t %>%
  distinct(sample_id, .keep_all = TRUE)

############################################################
# DEBUG: CHECK OVERLAP
############################################################

overlap <- intersect(sig$sample_id, expr_t$sample_id)

cat("Total signature samples:", nrow(sig), "\n")
cat("Total expression samples:", nrow(expr_t), "\n")
cat("Overlap:", length(overlap), "\n")

############################################################
# 🚨 HARD FILTER (CRITICAL)
############################################################

sig <- sig %>%
  filter(sample_id %in% overlap)

expr_t <- expr_t %>%
  filter(sample_id %in% overlap)

############################################################
# STRICT CHECK
############################################################

stopifnot(length(unique(sig$sample_id)) == length(unique(expr_t$sample_id)))
stopifnot(length(overlap) > 0)

############################################################
# SAVE CLEAN DATA
############################################################

saveRDS(sig, "output/processed/signatures_clean.rds")
saveRDS(expr_t, "output/processed/expression_clean.rds")

############################################################
# 3. PROCESS SIGNATURES → SBS3 PRESENCE
############################################################

sig <- sig %>%
  mutate(
    SBS3_status = ifelse(SBS3 > 0, "SBS3_positive", "SBS3_negative")
  )

############################################################
# QC PLOTS
############################################################

p1 <- ggplot(sig, aes(x = SBS3_status, fill = SBS3_status)) +
  geom_bar() +
  theme_minimal() +
  labs(
    title = "Sample Distribution by SBS3 Status",
    x = "",
    y = "Number of Samples"
  )
print(p1)
#ggsave("output/plots/SBS3_distribution.png", p1)

p2 <- ggplot(sig, aes(x = cancer_type, fill = SBS3_status)) +
  geom_bar(position = "fill") +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Cancer Type Composition",
    x = "Cancer Type",
    y = "Proportion"
  )

p2

#ggsave("output/plots/cancer_type_composition.png", p2)

############################################################
# 4. MERGE COHORT
############################################################

merged <- sig %>%
  inner_join(expr_t, by = "sample_id")

############################################################
# FINAL QC
############################################################

cat("Merged samples:", nrow(merged), "\n")

# 🔥 ensure no NA rows
merged <- merged %>%
  filter(!is.na(SBS3_status))

############################################################
# SAVE FINAL COHORT
############################################################

saveRDS(merged, "output/processed/merged_cohort.rds")

############################################################
# DONE
############################################################

cat("✅ Preprocessing complete (Shifu-grade) 🚀\n")

