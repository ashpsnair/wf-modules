#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args) != 2) {
  stop("Two arguments must be supplied: base_dir and sample_name", call.=FALSE)
}

base_dir <- args[1]
sample_name <- args[2]

library(sigminer)
library(tibble)
library(pheatmap)

input_matrix <- file.path(base_dir, "output", "SBS", paste0(sample_name, ".SBS96.all"))
output <- file.path(base_dir, "sigminer")

dir.create(output, showWarnings = FALSE)

matrix <- read.delim(input_matrix, sep = "\t")
matrix <- as.matrix(setNames(data.frame(t(matrix[,-1])), matrix[,1]))

model <- sig_unify_extract(
  matrix,
  range = 1:10,
  approach = "bayes_nmf",
  nrun = 10
)

signatures <- as.data.frame(model$Signature.norm)
signatures <- cbind(MutationsType = rownames(signatures), signatures)
write.table(signatures, file.path(output, "denovo_sig.tsv"), sep = "\t", quote = FALSE, row.name = FALSE)

exposure <- tibble::rownames_to_column(as.data.frame(t(model$Exposure)), "Sample")
write.table(exposure, file.path(output, "denovo_exposure.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

write.table(model$Raw$summary_run, file.path(output, "run.log"), sep = "\t", quote = FALSE, row.name = FALSE)

sim <- get_sig_similarity(model, sig_db = "latest_SBS_GRCh38")

png(file.path(output, "similarity_heatmap.png"), width = 1200, height = 350)   

pheatmap(sim$similarity,
         color = colorRampPalette(c("#348ABD", "white", "#E24A33"))(70),
         border_color = FALSE)

dev.off()

print("R analysis completed.")
