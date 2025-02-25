'''
#Installation

source /home/project/11003581/Tools/miniforge3/bin/activate

conda create --prefix /home/project/11003581/conda-envs/scMultiomics python=3.7 r-base=4.2.1

conda activate /home/project/11003581/conda-envs/scMultiomics

## fixing R issue 
mkdir -p /home/project/11003581/R/library/4.2/
echo "R_LIBS_USER=/home/project/11003581/R/library/4.2/" >> ~/.Renviron
echo ".libPaths("/home/project/11003581/R/library/4.2/")" >> ~/.Rprofile

#or edit the Rprofile 
vim ~/.Rprofile
#add .libPaths("/home/project/11003581/R/library/4.2/")
#check
.libPaths()

BiocManager::install("TFBSTools", lib = "/home/project/11003581/R/library/4.2/")

install.packages("remotes")
install.packages("Seurat")
remotes::install_github("satijalab/azimuth", ref = "master")
install.packages("SeuratData")
remotes::install_github("10XGenomics/loupeR")


''' 

###### CHanging the data #######
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7777520
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE232314

############# R scipt ###########

library (remotes)
library(Seurat)
library(Azimuth)
library(SeuratData)
library(loupeR)
loupeR::setup()

setwd("/home/users/nus/ash.ps/scratch/mulitomics/celltype-annotation/")

options(timeout = 1000)
### Load datasets
data= Read10X_h5("/home/users/nus/ash.ps/scratch/mulitomics/data/brca/GSE232314_filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)

# initialize a Seurat object of the data:
data<- CreateSeuratObject(data, project = "brca-trial", assay = "RNA")

#Run Azimuth (takes the feature barcode matrix and Seurat object reference as inputs)
Prediction<-RunAzimuth(data, "pbmcMultiome")

#other ref available- pbmcMultiome

#extract level 3 cell type annotations generated by Azimuth:
CellTypes_Azimuth<-Prediction@meta.data
CellTypes_Azimuth$Barcode<-rownames(CellTypes_Azimuth)
CellTypes_Azimuth$CustomLevel<-CellTypes_Azimuth$predicted.celltype.l2


#Create output that can be added back into the Seurat object:
CellTypes_Azimuth<-CellTypes_Azimuth[,c("Barcode","CustomLevel")]
colnames(CellTypes_Azimuth)<-c("Barcode","CustomLevel")

#Add cell type annotations to the Seurat object as metadata:
rownames(CellTypes_Azimuth)<-CellTypes_Azimuth$Barcode
Prediction[["CustomLevel"]]<-CellTypes_Azimuth$Barcode[match(rownames(Prediction@meta.data), CellTypes_Azimuth$CustomLevel)]
Prediction<-AddMetaData(object = Prediction, metadata = CellTypes_Azimuth$CustomLevel, col.name = "CellTypes")


##### Creating tSNE file
seurat_obj <- NormalizeData(data)

# Find variable features
seurat_obj <- FindVariableFeatures(seurat_obj)

# Scale the data
seurat_obj <- ScaleData(seurat_obj)

# Run PCA
seurat_obj <- RunPCA(seurat_obj)

# Run t-SNE
#wihtout removing duplicates
seurat_obj <- RunTSNE(seurat_obj, dims = 1:30, check_duplicates = FALSE)
#by removing duplicates
#seurat_obj <- rmDuplicateGenes(seurat_obj)

# Extract t-SNE coordinates
tsne_coords <- Embeddings(seurat_obj, reduction = "tsne")

# Create a data frame with cell barcodes and t-SNE coordinates
tsne_df <- data.frame(Barcode = rownames(tsne_coords),
                      tSNE_1 = tsne_coords[,1],
                      tSNE_2 = tsne_coords[,2])

# Save as CSV
write.csv(tsne_df, file = "t-SNE-Projection.csv", row.names = FALSE)



tSNE_coordinates <- read.csv("t-SNE-Projection.csv", stringsAsFactors = F, header = T, row.names = 1)
tSNE_coordinates_mat <- as(tSNE_coordinates, "matrix")
Prediction[['tSNE']] <- CreateDimReducObject(embeddings = tSNE_coordinates_mat, key = "tSNE_", global = T, assay = "RNA")


#Run create_loupe_from_seurat command, which takes the Seurat object as input:
create_loupe_from_seurat(Prediction)

#write file
write.table(CellTypes_Azimuth, "CellTypes_Azimuth.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

########### Running the above R code as PBS job ##############

#!/bin/bash

#PBS -l select=2:ncpus=64:mem=128g
#PBS -l walltime=4:00:00
#PBS -P 11003581
#PBS -j oe
#PBS -N celltype_annotation


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

# Run the R script
Rscript celltype_annotation.R > log.txt 2>&1
