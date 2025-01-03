'''
#Installation source- https://bioconductor.org/packages/devel/bioc/vignettes/maftools/inst/doc/cnv_analysis.html

'''

#!/bin/bash

#PBS -l select=2:ncpus=64:mem=256gb:ngpus=1
#PBS -l walltime=10:00:00
#PBS -P 11003581
#PBS -N ascat-p3l
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

# Load the R module
module load r

# Run the R script
Rscript -e "

library(ASCAT)

ascat.prepareHTS(
  tumourseqfile = '/home/project/11003581/Data/Ash/P3L-lab-analysis/sorted_bams/p3292l_hac_sorted.bam',
  normalseqfile = '/home/project/11003581/Data/Ash/P3L-lab-analysis/sorted_bams/WT_MCF10A_hac_sorted.bam',
  tumourname = 'p3l',
  normalname = 'wt-mcf10a',
  allelecounter_exe = allelecounter,
  skip_allele_counting_normal = FALSE,
  skip_allele_counting_tumour = FALSE,
  alleles.prefix = 'chr',
  loci.prefix = 'chr',
  gender = 'XY',
  genomeVersion = 'hg38',
  nthreads = 64,
  tumourLogR_file = 'tumor_logR.txt',
  tumourBAF_file = 'tumor_BAF.txt',
  normalLogR_file = 'normal_logR.txt',
  normalBAF_file = 'normal_BAF.txt',
  loci_binsize = 500,
  min_base_qual= 10,
  additional_allelecounter_flags='-f 0')

"