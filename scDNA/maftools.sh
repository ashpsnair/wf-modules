#!/bin/bash
#PBS -l select=3:ncpus=64:mem=256g
#PBS -l walltime=04:00:00
#PBS -P 11003581
#PBS -N maftools-3401
#PBS -j oe
#PBS -M ash.ps@nus.edu.sg
#PBS -m abe

cd $PBS_O_WORKDIR

module load r/4.2.0

Rscript -e '
library(maftools)
library(stringr)
library(NMF)
library(pheatmap)
library(BSgenome.Hsapiens.UCSC.hg38)

sample_name <- "3401"
base_dir <- "/home/users/nus/ash.ps/scratch/scDNA/analysis/annotation/mafs"
output_graphs <- file.path(base_dir, "3401-maftools-out")

cat("Creating output dir:", output_graphs, "\n")
dir.create(output_graphs, showWarnings = TRUE, recursive = TRUE)

if (!dir.exists(output_graphs)) stop("Output directory does not exist: ", output_graphs)

maf_file <- file.path(base_dir, paste0(sample_name, ".maf"))
cat("Reading MAF from:", maf_file, "\n")

if (!file.exists(maf_file)) stop("MAF file not found at path: ", maf_file)

laml <- read.maf(maf = maf_file)

png(file.path(output_graphs, "maf_summary_plot.png"), width = 1200, height = 800, res = 150)
plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = "median", dashboard = TRUE, titvRaw = FALSE)
dev.off()

png(file.path(output_graphs, "onco_plot.png"), width = 1200, height = 800, res = 150)
oncoplot(maf = laml, top = 25)
dev.off()

png(file.path(output_graphs, "titv.png"), width = 1200, height = 800, res = 150)
laml.titv <- titv(maf = laml, plot = TRUE, useSyn = TRUE)
dev.off()

png(file.path(output_graphs, "tmb-tcga.png"), width = 1200, height = 800, res = 150)
laml.mutload <- tcgaCompare(maf = laml, cohortName = sample_name , logscale = TRUE, capture_size = 50)
dev.off()

png(file.path(output_graphs, "somatic-interactions.png"), width = 1200, height = 800, res = 150)
somaticInteractions(maf = laml, top = 25, pvalue = c(0.05, 0.1))
dev.off()

png(file.path(output_graphs, "pfam.png"), width = 1200, height = 800, res = 150)
laml.pfam <- pfamDomains(maf = laml, AACol = "aaChange", top = 10)
dev.off()

png(file.path(output_graphs, "drug-gene.png"), width = 1200, height = 800, res = 150)
dgi <- drugInteractions(maf = laml, fontSize = 0.75)
dev.off()

png(file.path(output_graphs, "onco-pathways.png"), width = 1200, height = 800, res = 150)
OncogenicPathways(laml)
dev.off()

laml.tnm <- trinucleotideMatrix(maf = laml, add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg38")
png(file.path(output_graphs, "apobec-enrichment.png"), width = 1200, height = 800, res = 150)
plotApobecDiff(tnm = laml.tnm, maf = laml, pVal = 0.05)
dev.off()

laml.sig <- oncodrive(maf = laml, AACol = "aaChange", minMut = 5, pvalMethod = "zscore")
png(file.path(output_graphs, "oncodrive.png"), width = 1200, height = 800, res = 150)
plotOncodrive(res = laml.sig, fdrCutOff = 0.1, useFraction = TRUE, labelSize = 0.5)
dev.off()

png(file.path(output_graphs, "brca2-lollipop.png"), width = 1200, height = 800, res = 150)
lollipopPlot(maf = laml, gene = "BRCA2", AACol = "aaChange", showMutationRate = TRUE)
dev.off()
'
