#!/bin/bash

#PBS -l select=1:ncpus=32:mem=128g
#PBS -l walltime=1:00:00
#PBS -P 11003581
#PBS -N maftools-p3l
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load r/4.2.0

Rscript -e "

library(maftools)
counts = maftools::gtMarkers(t_bam = '/home/project/11003581/Data/Ash/P3L-lab-analysis/sorted_bams/p3292l_hac_sorted.bam',
                             n_bam = '/home/project/11003581/Data/Ash/P3L-lab-analysis/sorted_bams/WT_MCF10A_hac_sorted.bam',
                             build = 'hg38',
                             prefix = 'chr')

library(ASCAT)

ascat.bc = maftools::prepAscat(t_counts = '/home/project/11003581/Data/Ash/P3L-lab-analysis/R-maftools/p3292l_hac_sorted_nucleotide_counts.tsv',
                               n_counts = '/home/project/11003581/Data/Ash/P3L-lab-analysis/R-maftools/WT_MCF10A_hac_sorted_nucleotide_counts.tsv',
                               sample_name = 'tumor')

ascat.bc = ASCAT::ascat.loadData(
  Tumor_LogR_file = "tumor.tumour.logR.txt",
  Tumor_BAF_file = "tumor.tumour.BAF.txt",
  Germline_LogR_file = "tumor.normal.logR.txt",
  Germline_BAF_file = "tumor.normal.BAF.txt",
  chrs = c(1:22, "X", "Y"),
  sexchromosomes = c("X", "Y")
)

ASCAT::ascat.plotRawData(ASCATobj = ascat.bc, img.prefix = "tumor")
ascat.bc = ASCAT::ascat.aspcf(ascat.bc)
ASCAT::ascat.plotSegmentedData(ascat.bc)
ascat.output = ASCAT::ascat.runAscat(ascat.bc) 

maftools::segmentLogR(tumor_logR = "tumor.tumour.logR.txt", sample_name = "tumor")

"


'''
#!/bin/bash

#PBS -l select=1:ncpus=32:mem=128g
#PBS -l walltime=1:00:00
#PBS -P 11003581
#PBS -N plotmosdepth
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load r/4.2.0

Rscript -e "
library(maftools)
plotMosdepth(
  t_bed = "/home/project/11003581/Data/Ash/P3L-lab-analysis/mosdepth/p3l/p3l.regions.bed.gz",
  n_bed = "/home/project/11003581/Data/Ash/P3L-lab-analysis/mosdepth/wt/wt.regions.bed.gz",
  segment = TRUE,
  sample_name = "tumor"
)


"

'''