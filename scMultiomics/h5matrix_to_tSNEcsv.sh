# Load required libraries
library(Seurat)
library(Rtsne)

# Read the H5 file
data <- Read10X_h5("GSE232314_filtered_feature_bc_matrix.h5")

# Create a Seurat object
seurat_obj <- CreateSeuratObject(data, project = "brca-trial", assay = "RNA")

# Normalize the data
seurat_obj <- NormalizeData(seurat_obj)

# Find variable features
seurat_obj <- FindVariableFeatures(seurat_obj)

# Scale the data
seurat_obj <- ScaleData(seurat_obj)

# Run PCA
seurat_obj <- RunPCA(seurat_obj)

# Run t-SNE
seurat_obj <- RunTSNE(seurat_obj, dims = 1:30)

# Extract t-SNE coordinates
tsne_coords <- Embeddings(seurat_obj, reduction = "tsne")

# Create a data frame with cell barcodes and t-SNE coordinates
tsne_df <- data.frame(Barcode = rownames(tsne_coords),
                      tSNE_1 = tsne_coords[,1],
                      tSNE_2 = tsne_coords[,2])

# Save as CSV
write.csv(tsne_df, file = "t-SNE-Projection.csv", row.names = FALSE)
