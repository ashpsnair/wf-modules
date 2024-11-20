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

"