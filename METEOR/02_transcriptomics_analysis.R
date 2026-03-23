############################################################
# PCAWG SBS3 → METABOLIC ANALYSIS PIPELINE
# Shifu 🔥 FIXED VERSION
############################################################

library(tidyverse)
library(limma)
library(ggplot2)
library(pheatmap)
library(ggpubr)
library(purrr)

dir.create("output/analysis", showWarnings = FALSE)

############################################################
# 1. LOAD DATA
############################################################

data <- readRDS("output/processed/merged_cohort.rds")

############################################################
# 2. DEFINE GROUPS
############################################################

data$SBS3_status <- factor(
  data$SBS3_status,
  levels = c("SBS3_negative","SBS3_positive")
)

table(data$SBS3_status)

############################################################
# 3. BUILD EXPRESSION MATRIX (ROBUST FINAL)
############################################################

# metadata columns
meta_cols <- c("Samples", "sample_id", "cancer_type", "SBS3_status")

# remove metadata
expr <- data[, !(colnames(data) %in% meta_cols)]

# 🔥 REMOVE MUTATIONAL SIGNATURE COLUMNS
expr <- expr[, !grepl("^SBS", colnames(expr))]

cat("Columns after removing SBS:", ncol(expr), "\n")

# convert to numeric safely
expr <- expr %>%
  mutate(across(everything(), as.numeric))

# remove columns that are all NA
expr <- expr[, colSums(!is.na(expr)) > 0]

cat("Final gene count:", ncol(expr), "\n")

############################################################
# BUILD MATRIX
############################################################

expr_mat <- t(as.matrix(expr))

# assign sample IDs
colnames(expr_mat) <- data$sample_id

cat("Expression matrix dim:", dim(expr_mat), "\n")

############################################################
# REMOVE LOW INFORMATION GENES
############################################################

# remove zero variance genes
keep <- apply(expr_mat, 1, var, na.rm = TRUE) > 0
expr_mat <- expr_mat[keep, ]

cat("Genes after variance filter:", nrow(expr_mat), "\n")

############################################################
# 4. DESIGN MATRIX (SAFE ALIGNMENT)
############################################################

design <- model.matrix(~ 0 + SBS3_status, data = data)
colnames(design) <- c("SBS3_neg", "SBS3_pos")

# assign rownames
rownames(design) <- data$sample_id

############################################################
# ALIGN EXPRESSION + DESIGN
############################################################

common_ids <- intersect(colnames(expr_mat), rownames(design))

cat("Common samples:", length(common_ids), "\n")

# 🚨 hard stop if broken
if (length(common_ids) == 0) {
  stop("❌ No matching sample IDs between expression and design")
}

# subset safely
expr_mat <- expr_mat[, common_ids, drop = FALSE]
design <- design[common_ids, , drop = FALSE]

############################################################
# FINAL SANITY CHECK
############################################################

stopifnot(ncol(expr_mat) == nrow(design))

cat("Final matrix:", dim(expr_mat), "\n")

############################################################
# 5. LIMMA DIFFERENTIAL EXPRESSION
############################################################

fit <- lmFit(expr_mat, design)

contrast <- makeContrasts(SBS3_pos - SBS3_neg, levels = design)

fit <- contrasts.fit(fit, contrast)
fit <- eBayes(fit)

res <- topTable(fit, number = Inf)

res$gene <- rownames(res)

############################################################
# 6. VOLCANO PLOT (CLEAN)
############################################################

res <- res %>%
  mutate(
    significance = case_when(
      adj.P.Val < 0.05 & logFC > 1  ~ "Up",
      adj.P.Val < 0.05 & logFC < -1 ~ "Down",
      TRUE ~ "NS"
    )
  )

volcano <- ggplot(res, aes(logFC, -log10(P.Value), color = significance)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Up"="red","Down"="blue","NS"="grey")) +
  theme_minimal() +
  ggtitle("SBS3+ vs SBS3- (PCAWG)") +
  labs(color = "")

volcano

#ggsave("output/plots/volcano_SBS3.png", volcano, width = 6, height = 5)


############################################################
# 7. GSVA PATHWAY ANALYSIS (HALLMARK / KEGG)
############################################################

library(GSVA)
library(msigdbr)

############################################################
# LOAD GENE SETS
############################################################

# 🔥 OPTION 1: Hallmark (recommended)
msig <- msigdbr(species = "Homo sapiens", category = "H")

# 🔥 OPTION 2: KEGG (uncomment if needed)
# msig <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")

gene_sets <- msig %>%
  split(x = .$gene_symbol, f = .$gs_name)

############################################################
# PREPARE EXPRESSION MATRIX
############################################################
############################################################
# FIX GENE IDS FOR GSVA (ENSEMBL → SYMBOL)
############################################################

library(biomaRt)

# remove version numbers
gene_ids <- gsub("\\..*", "", rownames(expr_mat))

############################################################
# CONNECT TO ENSEMBL
############################################################

mart <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl"
)

############################################################
# MAP IDS
############################################################

mapping <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = gene_ids,
  mart = mart
)

############################################################
# MERGE WITH EXPRESSION
############################################################

expr_df <- as.data.frame(expr_mat)
expr_df$ensembl_gene_id <- gene_ids

expr_df <- expr_df %>%
  left_join(mapping, by = "ensembl_gene_id") %>%
  filter(hgnc_symbol != "")

############################################################
# HANDLE DUPLICATES (IMPORTANT)
############################################################

expr_df <- expr_df %>%
  group_by(hgnc_symbol) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop")

############################################################
# REBUILD MATRIX
############################################################

expr_mat_symbol <- as.matrix(expr_df[, -1])
rownames(expr_mat_symbol) <- expr_df$hgnc_symbol

cat("Final GSVA matrix:", dim(expr_mat_symbol), "\n")

# GSVA expects: genes x samples
expr_gsva <- expr_mat

############################################################
# CREATE PARAM OBJECT (NEW WAY)
############################################################

gsva_param <- gsvaParam(
  exprData = expr_gsva,
  geneSets = gene_sets,
  kcdf = "Gaussian"
)

############################################################
# RUN GSVA
############################################################

gsva_scores <- gsva(gsva_param)

############################################################
# CONVERT TO DATAFRAME
############################################################

gsva_df <- as.data.frame(t(gsva_scores))
gsva_df$sample_id <- rownames(gsva_df)

############################################################
# MERGE WITH SBS3 STATUS
############################################################

gsva_df <- gsva_df %>%
  left_join(data %>% select(sample_id, SBS3_status), by = "sample_id")

############################################################
# 8. DIFFERENTIAL PATHWAY ANALYSIS (LIMMA AGAIN 🔥)
############################################################

gsva_mat <- t(as.matrix(gsva_df %>% select(-sample_id, -SBS3_status)))

design_gsva <- model.matrix(~ 0 + SBS3_status, data = gsva_df)
colnames(design_gsva) <- c("SBS3_neg", "SBS3_pos")

fit_gsva <- lmFit(gsva_mat, design_gsva)

contrast <- makeContrasts(SBS3_pos - SBS3_neg, levels = design_gsva)

fit_gsva <- contrasts.fit(fit_gsva, contrast)
fit_gsva <- eBayes(fit_gsva)

gsva_res <- topTable(fit_gsva, number = Inf)
gsva_res$pathway <- rownames(gsva_res)

write.csv(gsva_res, "output/gsva_pathway_results.csv", row.names = FALSE)

############################################################
# 9. PLOT TOP PATHWAYS
############################################################

top_pathways <- gsva_res %>%
  arrange(adj.P.Val) %>%
  head(12)

plot_df <- gsva_df %>%
  select(sample_id, SBS3_status, all_of(top_pathways$pathway)) %>%
  pivot_longer(
    -c(sample_id, SBS3_status),
    names_to = "pathway",
    values_to = "score"
  )

p_pathway <- ggplot(plot_df, aes(SBS3_status, score, fill = SBS3_status)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.3) +
  facet_wrap(~pathway, scales = "free") +
  theme_minimal() +
  ggtitle("GSVA Pathway Activity (Top pathways)")

ggsave("output/plots/GSVA_pathways.png", p_pathway, width = 10, height = 7)