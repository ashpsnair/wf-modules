#!/bin/bash

#PBS -l select=2:ncpus=64:mem=128g
#PBS -l walltime=10:00:00
#PBS -P 11003581
#PBS -j oe
#PBS -N run-scMulti


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
project="bc"
base_dir="/home/users/nus/ash.ps/scratch/mulitomics/data/breast-cancer"
output_dir="/home/users/nus/ash.ps/scratch/mulitomics/analysis/trial1"
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
print("Seurat Object Created")

# Normalize the data
seurat_obj <- NormalizeData(seurat_obj, assay = "RNA")
print("Seurat Object Normalized")

############ START: classification based on gene expression ############

# Define gene sets
genes_of_interest <- c("BRCA2", "BRCA1", "HMOX1", "NRF1")

# Remove empty layers
for (assay in Assays(seurat_obj)) {
  for (layer in Layers(seurat_obj[[assay]])) {
    if (ncol(seurat_obj[[assay]][[layer]]) == 0 || ncol(seurat_obj[[assay]][[layer]]) == 1) {
      seurat_obj[[assay]][[layer]] <- NULL
      print(paste("Removed empty layer:", layer, "from assay:", assay))
    }
  }
}

# Check which genes are actually present in the dataset
available_genes <- intersect(genes_of_interest, rownames(seurat_obj[["RNA"]]))
print("Available genes:")
print(available_genes)

# Calculate expression for available genes
gene_expression <- seurat_obj[["RNA"]]$data[available_genes, , drop = FALSE]
gene_means <- rowMeans(gene_expression)

# Function to determine gene status
determine_gene_status <- function(expr, gene_name, mean_expr) {
  if (is.null(expr) || length(expr) == 0) return(paste0(gene_name, "_absent"))
  if (all(expr == 0)) return(paste0(gene_name, "_null"))
  
  deciles <- quantile(expr[expr > 0], probs = seq(0, 1, 0.1))
  
  if (mean(expr) >= deciles[8]) return(paste0(gene_name, "_high"))
  if (mean(expr) >= deciles[4] && mean(expr) <= deciles[6]) return(paste0(gene_name, "_low"))
  return(paste0(gene_name, "_other"))
}

# Apply gene status determination for available genes
for (gene in available_genes) {
  status_column <- paste0(gene, "_status")
  seurat_obj[[status_column]] <- apply(gene_expression[gene, , drop = FALSE], 2, determine_gene_status, gene_name = gene, mean_expr = gene_means[gene])
}

# Combine all statuses
seurat_obj$combined_status <- do.call(paste, c(lapply(available_genes, function(gene) seurat_obj[[paste0(gene, "_status")]]), sep = "-"))

# Create combined classification dataframe
combined_classification <- data.frame(
  Index = colnames(seurat_obj),
  Cell_type = seurat_obj$combined_status
)

# Write combined classification
write.table(combined_classification, "cell_classification.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


############ END: classification based on gene expression ############

# Print summary to console
print(stats_df)

# Calculate and print total number of cells
total_cells <- ncol(seurat_obj)
cat("\nTotal number of cells:", total_cells, "\n")

EOF

Rscript "$output_dir/cell_type_annotation.R" > log.txt 2>&1

echo "Cell type annotation completed successfully"

#############################################################
#Step 1: Splitting alignment file into cell-type-specific bams
#############################################################

output_dir1=$output_dir/Step1_BamCellTypes
mkdir -p $output_dir1

python $SCOMATIC/scripts/SplitBam/SplitBamCellTypes.py \
        --bam $bam_file \
        --meta ./cell_classification.tsv \
        --id $project \
        --n_trim 5 \
        --max_nM 5 \
        --max_NH 1 \
        --outdir $output_dir1

echo "BAM splitting completed successfully"

# Filtering and indexing BAM files
filtered_out=$output_dir1/filtered
mkdir -p "$filtered_out" 

# List of chromosomes to keep
chromosomes="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"

# Process each BAM file in the input directory
for input_bam in "$output_dir1"/*.bam; do
    base_name=$(basename "$input_bam")
    output_bam="${filtered_out}/${base_name%.bam}_filtered.bam"
    
    samtools view -h "$input_bam" $chromosomes | samtools sort -o "$output_bam"
    samtools index "$output_bam"
    
    echo "Processed and indexed: $base_name"
done

echo "All BAM files have been filtered and indexed."

########## STEP 2 ###########

REF=/home/project/11003581/Ref/Homo_sapiens/GATK/hg38/Homo_sapiens_assembly38.fasta
editing=$SCOMATIC/RNAediting/AllEditingSites.hg38.txt
PON=$SCOMATIC/PoNs/PoN.scRNAseq.hg38.tsv

output_dir1=$output_dir/Step1_BamCellTypes/filtered/
output_dir2=$output_dir/Step2_BaseCellCounts
output_dir4=$output_dir/Step4_VariantCalling
mkdir -p $output_dir2 $output_dir4

for bam in $(ls -d $output_dir1/*bam); do
  cell_type=$(basename $bam | awk -F'.' '{print $(NF-1)}')
  temp=$output_dir2/temp_${cell_type}
  mkdir -p $temp

  # Step 2: Collecting base count information
  python $SCOMATIC/scripts/BaseCellCounter/BaseCellCounter.py --bam $bam \
    --ref $REF \
    --chrom all \
    --out_folder $output_dir2 \
    --min_bq 30 \
    --tmp_dir $temp \
    --nprocs 128

########## STEP 3 ###########

  # Step 3.1: Detection of somatic mutations
  python $SCOMATIC/scripts/BaseCellCalling/BaseCellCalling.step1.py \
    --infile ${output_dir2}/${cell_type}.BaseCellCounts.tsv \
    --outfile ${output_dir4}/${cell_type} \
    --ref $REF

  # Step 3.2: Further processing
  python $SCOMATIC/scripts/BaseCellCalling/BaseCellCalling.step2.py \
    --infile ${output_dir4}/${cell_type}.calling.step1.tsv \
    --outfile ${output_dir4}/${cell_type} \
    --editing $editing \
    --pon $PON

  rm -rf $temp
done
