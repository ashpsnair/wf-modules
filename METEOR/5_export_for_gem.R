library(data.table)

expr <- readRDS("results/expr_tpm_matched.rds")
labels <- fread("results/hrd/hrd_labels.tsv")

expr_patients <- substr(colnames(expr), 1, 12)
keep <- expr_patients %in% labels$patient_id
expr <- expr[, keep, drop = FALSE]

# one sample per patient
ord <- !duplicated(substr(colnames(expr), 1, 12))
expr <- expr[, ord, drop = FALSE]

colnames(expr) <- substr(colnames(expr), 1, 12)

write.table(
  data.frame(gene = rownames(expr), expr, check.names = FALSE),
  file = "results/gem_input_tpm_matrix.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  labels,
  file = "results/gem_sample_metadata.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)