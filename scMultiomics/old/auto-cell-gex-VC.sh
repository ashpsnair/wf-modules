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
output_dir="/home/users/nus/ash.ps/scratch/mulitomics/analysis/breast-cancer"
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

### Load datasets
data <- Read10X_h5("$h5_file", use.names = TRUE, unique.features = TRUE)

# initialize a Seurat object of the data:
seurat_obj <- CreateSeuratObject(data, project = "$project", assay = "RNA")

# Run Azimuth for cell type annotation
Prediction <- RunAzimuth(seurat_obj, reference = "pbmcref")


#extract level 3 cell type annotations generated by Azimuth:
CellTypes_Azimuth<-Prediction@meta.data
CellTypes_Azimuth\$Barcode<-rownames(CellTypes_Azimuth)
CellTypes_Azimuth\$CustomLevel<-CellTypes_Azimuth\$predicted.celltype.l2


#Create output that can be added back into the Seurat object:
CellTypes_Azimuth<-CellTypes_Azimuth[,c("Barcode","CustomLevel")]
colnames(CellTypes_Azimuth)<-c("Barcode","CustomLevel")


# Write file
write.table(CellTypes_Azimuth, "CellTypes_Azimuth.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Normalize the data
seurat_obj <- NormalizeData(seurat_obj, assay = "RNA")

# Check if PRKAA1 and PRKAA2 are present in the dataset
ampk_genes <- c("PRKAA1", "PRKAA2")
present_genes <- ampk_genes[ampk_genes %in% rownames(seurat_obj)]

if (length(present_genes) > 0) {
  # Calculate mean expression for present AMPK genes
  mean_expression <- rowMeans(GetAssayData(seurat_obj, assay = "RNA", slot = "counts")[present_genes,])
  
  # Create a new metadata column for AMPK classification
  seurat_obj\$AMPK_status <- ifelse(
    GetAssayData(seurat_obj, assay = "RNA", slot = "counts")["PRKAA1",] > mean_expression["PRKAA1"] |
    GetAssayData(seurat_obj, assay = "RNA", slot = "counts")["PRKAA2",] > mean_expression["PRKAA2"],
    "AMPK_positive", "AMPK_negative"
  )

# Create a dataframe with barcodes and AMPK classification
  ampk_classification <- data.frame(
    Index = rownames(seurat_obj@meta.data),
    Cell_type = seurat_obj\$AMPK_status
  )
  
  write.table(ampk_classification, "gex-classification.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
} else {
  print("AMPK genes (PRKAA1 and PRKAA2) not found in the dataset")
}

### Merging barcodes
merged_data <- merge(CellTypes_Azimuth, ampk_classification, by = "Index", all = TRUE, suffixes = c("_azimuth", "_ampk"))

# Combine Cell_type values
merged_data\$Cell_type <- ifelse(!is.na(merged_data\$Cell_type_azimuth) & !is.na(merged_data\$Cell_type_ampk),
                                paste(merged_data\$Cell_type_azimuth, merged_data\$Cell_type_ampk, sep = "-"),
                                ifelse(!is.na(merged_data\$Cell_type_azimuth), merged_data\$Cell_type_azimuth, merged_data\$Cell_type_ampk))

# Select only the required columns
result <- merged_data[, c("Index", "Cell_type")]

# Write the result to a file
write.table(result, "merged_classification.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Generate stats.txt
stats <- data.frame(
  Dataset = c(rep("CellTypes_Azimuth", length(table(CellTypes_Azimuth\$Cell_type))),
              rep("ampk_classification", length(table(ampk_classification\$Cell_type))),
              rep("result", length(table(result\$Cell_type)))),
  Cell_type = c(names(table(CellTypes_Azimuth\$Cell_type)),
                names(table(ampk_classification\$Cell_type)),
                names(table(result\$Cell_type))),
  Count = c(as.vector(table(CellTypes_Azimuth\$Cell_type)),
            as.vector(table(ampk_classification\$Cell_type)),
            as.vector(table(result\$Cell_type)))
)

EOF

Rscript "$output_dir/cell_type_annotation.R" > log.txt 2>&1

echo "Cell type annotation completed successfully"

##################################################
### Step 1: Splitting alignment file into cell-type-specific bams
##################################################

output_dir1=$output_dir/Step1_BamCellTypes
mkdir -p $output_dir1

python $SCOMATIC/scripts/SplitBam/SplitBamCellTypes.py \
        --bam $bam_file \
        --meta ./gex-classification.tsv \
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

### Step 5 intersect

#!/bin/bash

#PBS -l select=2:ncpus=64:mem=128g
#PBS -l walltime=6:00:00
#PBS -P 11003581
#PBS -j oe
#PBS -N run-scomatic-intersect


# Change to the directory where the job was submitted
cd $PBS_O_WORKDIR

source /home/project/11003581/Tools/miniforge3/bin/activate
conda activate /home/project/11003581/conda-envs/SComatic

output_dir="/home/users/nus/ash.ps/scratch/mulitomics/analysis/breast-cancer"
SCOMATIC=/home/project/11003581/Tools/SComatic/

output_dir4=$output_dir/Step4_VariantCalling

module load bedtools/2.30.0
bedtools intersect -header -a ${output_dir4}/ampk-neg.step2.tsv -b $SCOMATIC/bed_files_of_interest/UCSC.k100_umap.without.repeatmasker.bed | awk '$1 ~ /^#/ || $6 == "PASS"' > ${output_dir4}/ampk-neg.step2.pass.tsv