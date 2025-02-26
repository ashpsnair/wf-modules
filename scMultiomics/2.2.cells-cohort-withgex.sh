
######## R code for gex classification; saved as gex-class.R ####

library(Seurat)

setwd("/home/users/nus/ash.ps/scratch/mulitomics/gex-class/")
output_folder <- "/home/users/nus/ash.ps/scratch/mulitomics/gex-class/"
output_file <- file.path(output_folder, "gex-classification.tsv")

# Load the data
data = Read10X_h5("/home/users/nus/ash.ps/scratch/mulitomics/analysis/10k_pbmc/outs/filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)

# Create a Seurat object using the gene expression data and correct name of the assay:
seurat_obj = CreateSeuratObject(data$`Gene Expression`, project = "10kpbmc", assay = "Gene Expression")

# Normalize the data
seurat_obj <- NormalizeData(seurat_obj, assay = "Gene Expression")

# Check if PRKAA1 and PRKAA2 are present in the dataset
ampk_genes <- c("PRKAA1", "PRKAA2")
present_genes <- ampk_genes[ampk_genes %in% rownames(seurat_obj)]

if (length(present_genes) > 0) {
  # Calculate mean expression for present AMPK genes
  mean_expression <- rowMeans(GetAssayData(seurat_obj, assay = "Gene Expression", layer = "counts")[present_genes,])
  
  # Create a new metadata column for AMPK classification
  seurat_obj$AMPK_status <- ifelse(
    GetAssayData(seurat_obj, assay = "Gene Expression", layer = "counts")["PRKAA1",] > mean_expression["PRKAA1"] |
    GetAssayData(seurat_obj, assay = "Gene Expression", layer = "counts")["PRKAA2",] > mean_expression["PRKAA2"],
    "AMPK_positive", "AMPK_negative"
  )
  
  # Create a dataframe with barcodes and AMPK classification
  ampk_classification <- data.frame(
    Index = rownames(seurat_obj@meta.data),
    Cell_type = seurat_obj$AMPK_status
  )
  
} else {
  print("AMPK genes (PRKAA1 and PRKAA2) not found in the dataset")
}


write.table(ampk_classification, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

