#!/bin/bash

#PBS -l select=2:ncpus=64:mem=128g
#PBS -l walltime=10:00:00
#PBS -P 11003581
#PBS -j oe
#PBS -N run-auto-scMulti


# Change to the directory where the job was submitted
cd $PBS_O_WORKDIR

source /home/project/11003581/Tools/miniforge3/bin/activate
conda activate /home/project/11003581/conda-envs/scMultiomics

module load gcc
module load intel/2024.0
module load libfabric/1.11.0.4.125
module load hdf5/1.12.1-parallel-icc23-cmpi
module load gsl/2.7.1-gcc11
module load r/4.2.0 
module load samtools/1.15.1
module load bedtools/2.30.0
module load python/3.12.1-gcc11


# Set base directory and output directory
project="breast-cancer"
base_dir="/home/users/nus/ash.ps/scratch/mulitomics/data/breast-cancer"
output_dir="/home/users/nus/ash.ps/scratch/mulitomics/analysis/new"
SCOMATIC=/home/project/11003581/Tools/SComatic/

# Find and store file paths
h5_file=$(find "$base_dir" -name "*.h5" -type f)
bam_file=$(find "$base_dir" -name "*.bam" -type f)
bai_file=$(find "$base_dir" -name "*.bam.bai" -type f)

# Create R script
cat << EOF > "$output_dir/cell_type_annotation.R"
library(Seurat)
library(Azimuth)
library(SeuratData)

setwd("$output_dir")

options(timeout = 1000)

# Load dataset
data <- Read10X_h5("$h5_file", use.names = TRUE, unique.features = TRUE)

# Initialize Seurat object
seurat_obj <- CreateSeuratObject(data, project = "$project", assay = "RNA")

# Normalize the data
seurat_obj <- NormalizeData(seurat_obj, assay = "RNA")

# Define gene sets
ampk_up_targets <- c("PPARGC1A", "TRARG1", "SLC2A4")
ampk_down_targets <- c("ACC", "SREBF1")
ampk_genes <- c("PRKAA1", "PRKAA2")

# Calculate mean expression for each gene
gene_means <- rowMeans(seurat_obj@assays$RNA@data[c(ampk_up_targets, ampk_down_targets, ampk_genes), ])

# Function to check gene regulation
check_regulation <- function(gene_expr, gene_mean) {
  return(ifelse(gene_expr > gene_mean, "up", ifelse(gene_expr < gene_mean, "down", "no_change")))
}

# Apply regulation check
seurat_obj <- AddMetaData(seurat_obj, 
                          apply(seurat_obj@assays$RNA@data[c(ampk_up_targets, ampk_down_targets, ampk_genes), ], 2, 
                                function(x) sapply(seq_along(x), function(i) check_regulation(x[i], gene_means[i]))),
                          col.name = c(ampk_up_targets, ampk_down_targets, ampk_genes))

# Determine AMPK status
determine_ampk_status <- function(cell_data) {
  if (all(cell_data[ampk_up_targets] == "up") && all(cell_data[ampk_down_targets] == "down")) {
    return("ampk_activated")
  } else if (all(cell_data[c(ampk_up_targets, ampk_down_targets)] == "no_change")) {
    return("ampk_inactive")
  } else if (all(cell_data[ampk_genes] == "up")) {
    return("ampk_overexpressed")
  } else if (all(cell_data[ampk_genes] == "down")) {
    return("ampk_underexpressed")
  } else {
    return("other")
  }
}

# Apply AMPK status determination
seurat_obj$ampk_status <- apply(seurat_obj@meta.data[, c(ampk_up_targets, ampk_down_targets, ampk_genes)], 1, determine_ampk_status)

# Calculate mean expression for BRCA2
brca2_mean <- mean(seurat_obj@assays$RNA@data["BRCA2",])

# Determine BRCA2 status
determine_brca2_status <- function(brca2_expr) {
  if (brca2_expr == 0) {
    return("brca2_null")
  } else if (abs(brca2_expr - brca2_mean) / brca2_mean < 0.1) {
    return("brca2_het")
  } else {
    return("brca2_other")
  }
}

# Apply BRCA2 status determination
seurat_obj$brca2_status <- sapply(seurat_obj@assays$RNA@data["BRCA2",], determine_brca2_status)

# Combine AMPK and BRCA2 status
seurat_obj$combined_status <- paste(seurat_obj$ampk_status, seurat_obj$brca2_status, sep = "-")

# Create combined classification dataframe
combined_classification <- data.frame(
  Barcode = rownames(seurat_obj@meta.data),
  Combined_status = seurat_obj$combined_status
)

# Write combined classification
write.table(combined_classification, "combined_classification.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Merge with Azimuth results
merged_data <- merge(CellTypes_Azimuth, combined_classification, by = "Barcode", all = TRUE)

# Combine Cell_type values
merged_data$Cell_type <- paste(merged_data$CustomLevel, merged_data$Combined_status, sep = "-")

# Select required columns
result <- merged_data[, c("Barcode", "Cell_type")]

# Write final result
write.table(result, "merged_classification.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Generate stats
stats <- data.frame(
  Dataset = c(rep("CellTypes_Azimuth", length(table(CellTypes_Azimuth$CustomLevel))),
              rep("combined_classification", length(table(combined_classification$Combined_status))),
              rep("result", length(table(result$Cell_type)))),
  Cell_type = c(names(table(CellTypes_Azimuth$CustomLevel)),
                names(table(combined_classification$Combined_status)),
                names(table(result$Cell_type))),
  Count = c(as.vector(table(CellTypes_Azimuth$CustomLevel)),
            as.vector(table(combined_classification$Combined_status)),
            as.vector(table(result$Cell_type)))
)

# Write stats
write.table(stats, "stats.txt", sep = "\t", row.names = FALSE, quote = FALSE)

EOF

Rscript "$output_dir/cell_type_annotation.R" > log.txt 2>&1

echo "Cell type annotation completed successfully"
