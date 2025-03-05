#!/bin/bash

#PBS -l select=2:ncpus=64:mem=128g
#PBS -l walltime=01:00:00
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
ampk_down_targets <- c("ACACA", "SREBF1")
ampk_genes <- c("PRKAA1", "PRKAA2")

# Calculate mean expression for AMPK-related genes
ampk_genes_all <- c(ampk_up_targets, ampk_down_targets, ampk_genes)
ampk_expression <- LayerData(seurat_obj, assay = "RNA", layer = "data")[ampk_genes_all, ]
gene_means <- rowMeans(ampk_expression)

# Function to determine AMPK status
determine_ampk_status <- function(cell_expr) {
  up_regulated <- all(cell_expr[ampk_up_targets] > gene_means[ampk_up_targets])
  down_regulated <- all(cell_expr[ampk_down_targets] < gene_means[ampk_down_targets])
  overexpressed <- all(cell_expr[ampk_genes] > gene_means[ampk_genes])
  underexpressed <- all(cell_expr[ampk_genes] < gene_means[ampk_genes])
  
  if (up_regulated && down_regulated) {
    return("ampk_activated")
  } else if (all(abs(cell_expr - gene_means) / gene_means < 0.1)) {
    return("ampk_inactive")
  } else if (overexpressed) {
    return("ampk_overexpressed")
  } else if (underexpressed) {
    return("ampk_underexpressed")
  } else {
    return("ampk_other")
  }
}

# Apply AMPK status determination
seurat_obj$ampk_status <- apply(ampk_expression, 2, determine_ampk_status)

# Calculate mean expression for BRCA2
brca2_expression <- LayerData(seurat_obj, assay = "RNA", layer = "data")["BRCA2", ]
brca2_mean <- mean(brca2_expression)

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
seurat_obj$brca2_status <- sapply(brca2_expression, determine_brca2_status)


# Combine AMPK and BRCA2 status
seurat_obj$combined_status <- paste(seurat_obj$ampk_status, seurat_obj$brca2_status, sep = "-")

# Create combined classification dataframe
combined_classification <- data.frame(
  Index = rownames(seurat_obj@meta.data),
  Cell_type = seurat_obj$combined_status
)

# Write combined classification
write.table(combined_classification, "cell_classification.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


# Generate statistics on the number of cells for each category
ampk_stats <- table(seurat_obj$ampk_status)
brca2_stats <- table(seurat_obj$brca2_status)
combined_stats <- table(seurat_obj$combined_status)

# Create a data frame with the statistics
stats_df <- data.frame(
  Category = c(names(ampk_stats), names(brca2_stats), names(combined_stats)),
  Count = c(ampk_stats, brca2_stats, combined_stats),
  Type = c(rep("AMPK", length(ampk_stats)), 
           rep("BRCA2", length(brca2_stats)), 
           rep("Combined", length(combined_stats)))
)

# Write statistics to a file
write.table(stats_df, "cell_category_stats.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

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



##################################################
########## Step 2 Collecting base count information
##################################################

REF=/home/project/11003581/Ref/Homo_sapiens/GATK/hg38/Homo_sapiens_assembly38.fasta

output_dir1=$output_dir/Step1_BamCellTypes/filtered/
output_dir2=$output_dir/Step2_BaseCellCounts
mkdir -p $output_dir2

for bam in $(ls -d $output_dir1/*bam);do
  
  # Cell type
  cell_type=$(basename $bam | awk -F'.' '{print $(NF-1)}')

  # Temp folder
  temp=$output_dir2/temp_${cell_type}
  mkdir -p $temp

  # Command line to submit to cluster
  python $SCOMATIC/scripts/BaseCellCounter/BaseCellCounter.py --bam $bam \
    --ref $REF \
    --chrom all \
    --out_folder $output_dir2 \
    --min_bq 30 \
    --tmp_dir $temp \
    --nprocs 16

  rm -rf $temp
done


##################################################
##########Step 4: Detection of somatic mutations
##################################################
# Step 4.1
output_dir4=$output_dir/Step4_VariantCalling
mkdir -p $output_dir4


python $SCOMATIC/scripts/BaseCellCalling/BaseCellCalling.step1.py \
          --infile ${output_dir3}/${project}.BaseCellCounts.AllCellTypes.tsv \
          --outfile ${output_dir4}/${project} \
          --ref $REF

# Step 4.2
editing=$SCOMATIC/RNAediting/AllEditingSites.hg38.txt
PON=$SCOMATIC/PoNs/PoN.scRNAseq.hg38.tsv

python $SCOMATIC/scripts/BaseCellCalling/BaseCellCalling.step2.py \
          --infile ${output_dir4}/${project}.calling.step1.tsv \
          --outfile ${output_dir4}/${project} \
          --editing $editing \
          --pon $PON
