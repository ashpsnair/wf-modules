'''
## Bioconductor packages required by mutSigExtractor
install.packages('BiocManager')
BiocManager::install('BSgenome') ## Install genome parser
BiocManager::install('BSgenome.Hsapiens.UCSC.hg19') ## Install the default genome

## randomForest is required by CHORD
install.packages('randomForest')

## Install CHORD and mutSigExtractor directly from github using devtools
install.packages("devtools")
devtools::install_github('https://github.com/UMCUGenetics/mutSigExtractor/')
devtools::install_github('https://github.com/UMCUGenetics/CHORD/')

## Make sure to install and load the desired ref genome first
install.packages('BiocManager')
BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')

'''

library(BSgenome.Hsapiens.UCSC.hg38)
library(CHORD)

#CHORD can be supplied with vcf files containing somatic (no germline) SNVs, indels, and SVs per sample.
contexts <- extractSigsChord(
  vcf.snv = '/path/to/vcf/with/snvs/',
  vcf.sv = '/path/to/vcf/with/svs/',
  sv.caller = 'gridss', ## Can be 'gridss' (default) or 'manta'

## This is optional (for when your vcfs are not filtered for PASS variants)
  vcf.filters=list(snv="PASS",indel="PASS",sv="PASS"),
  ref.genome=BSgenome.Hsapiens.UCSC.hg38
)

#Predicting HRD
chordPredict(contexts)


######### Running for multiple samples
## Get file paths from the vcf/ subdirectory
vcf_files <- data.frame(
  snv_indel=list.files('vcf', pattern='*snv_indel.vcf.gz', full.names=T),
  sv=list.files('vcf', pattern='*sv.vcf.gz', full.names=T)
)

## Assign sample names to files
vcf_files$sample <- sapply(strsplit(basename(vcf_files$snv_indel),'_'),`[`,1)

vcf_files
