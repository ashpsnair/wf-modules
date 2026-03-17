library(data.table)
library(dplyr)
library(tibble)
library(stringr)
library(GSVA)
library(msigdbr)
library(ggplot2)
library(pheatmap)
library(limma)

expr <- readRDS("results/expr_tpm_matched.rds")
labels <- fread("results/hrd/hrd_labels.tsv")

# ----------------------------------------
# match expression columns to labeled patients
# ----------------------------------------
expr_patients <- substr(colnames(expr), 1, 12)
keep <- expr_patients %in% labels$patient_id
expr <- expr[, keep, drop = FALSE]

# use one sample per patient if duplicates exist
ord <- !duplicated(substr(colnames(expr), 1, 12))
expr <- expr[, ord, drop = FALSE]

expr_patients <- substr(colnames(expr), 1, 12)
labels <- labels %>%
  filter(patient_id %in% expr_patients) %>%
  arrange(match(patient_id, expr_patients))

stopifnot(all(labels$patient_id == expr_patients))

# ----------------------------------------
# log transform TPM
# ----------------------------------------
expr_log2 <- log2(as.matrix(expr) + 1)

# ----------------------------------------
# build metabolic gene sets
# ----------------------------------------
msig <- msigdbr(species = "Homo sapiens", category = "C2")

met_sets_df <- msig %>%
  filter(
    str_detect(gs_name, "REACTOME_|KEGG_"),
    str_detect(gs_name, "GLYCOL|TCA|CITRATE|PENTOSE|FATTY_ACID|CHOLESTEROL|GLUTATHIONE|ONE_CARBON|FOLATE|PURINE|PYRIMIDINE|NUCLEOTIDE|TRYPTOPHAN|GLUTAMATE|AMINO_ACID|OXIDATIVE")
  ) %>%
  select(gs_name, gene_symbol) %>%
  distinct()

met_sets <- split(met_sets_df$gene_symbol, met_sets_df$gs_name)

# ----------------------------------------
# ssGSEA scores
# ----------------------------------------
ssgsea_scores <- gsva(
  expr = expr_log2,
  gset.idx.list = met_sets,
  method = "ssgsea",
  kcdf = "Gaussian",
  abs.ranking = TRUE
)

saveRDS(ssgsea_scores, "results/hrd/metabolic_ssgsea_scores.rds")

# ----------------------------------------
# correlate pathways with continuous HRD score
# ----------------------------------------
cor_df <- apply(ssgsea_scores, 1, function(x) {
  ct <- suppressWarnings(cor.test(x, labels$HRD_sum, method = "spearman"))
  c(rho = unname(ct$estimate), p = ct$p.value)
}) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("pathway") %>%
  mutate(FDR = p.adjust(p, method = "BH")) %>%
  arrange(FDR)

write.table(cor_df, "results/hrd/pathway_hrd_correlations.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

# ----------------------------------------
# compare HRD_high vs HRD_low
# ----------------------------------------
sel <- labels$HRD_group %in% c("HRD_high", "HRD_low")
mat_sel <- ssgsea_scores[, sel, drop = FALSE]
grp <- factor(labels$HRD_group[sel], levels = c("HRD_low", "HRD_high"))

design <- model.matrix(~ grp)
fit <- lmFit(mat_sel, design)
fit <- eBayes(fit)

res <- topTable(fit, coef = "grpHRD_high", number = Inf) %>%
  rownames_to_column("pathway")

write.table(res, "results/hrd/pathway_hrd_high_vs_low.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

# ----------------------------------------
# plots
# ----------------------------------------
top_cor <- cor_df %>%
  filter(FDR < 0.1) %>%
  slice_head(n = 20)

p1 <- ggplot(top_cor, aes(reorder(pathway, rho), rho)) +
  geom_col() +
  coord_flip() +
  theme_bw(base_size = 11) +
  labs(title = "Metabolic pathway correlation with HRD", x = "", y = "Spearman rho")

ggsave("results/hrd/top_pathway_correlations.png", p1, width = 8, height = 6, dpi = 300)

top_heat <- res %>%
  slice_head(n = 20) %>%
  pull(pathway)

ann <- data.frame(HRD_group = grp)
rownames(ann) <- colnames(mat_sel)

png("results/hrd/top20_pathway_heatmap.png", width = 1600, height = 1200, res = 180)
pheatmap(
  mat_sel[top_heat, , drop = FALSE],
  scale = "row",
  annotation_col = ann,
  show_colnames = FALSE,
  fontsize_row = 8
)
dev.off()