#!/bin/bash

#PBS -l select=1:ncpus=32:mem=128g
#PBS -l walltime=1:00:00
#PBS -P 11003581
#PBS -N p3l-coverage
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load samtools

#samtools view -b /home/project/11003581/Data/Ash/P3L-lab-analysis/sorted_bams/p3292l_hac_sorted.bam chr13 chr12 > /home/project/11003581/Data/Ash/P3L-lab-analysis/coverage/chr12_13_reads.bam

#samtools depth /home/project/11003581/Data/Ash/P3L-lab-analysis/coverage/chr12_13_reads.bam > /home/project/11003581/Data/Ash/P3L-lab-analysis/coverage/chr12_13_coverage.txt

awk 'BEGIN { bin_size=10000; } 
     { 
       bin=int($2/bin_size); 
       sum[bin]+=$3; 
       count[bin]++; 
     } 
     END { 
       for (b in sum) 
         printf "chr12\t%d\t%.2f\n", b*bin_size, sum[b]/count[b]; 
     }' chr12_13_coverage.txt > binned_coverage.txt


#python script
import matplotlib.pyplot as plt

# Load coverage data
coverage_data = []
with open('chr12_13_coverage.txt', 'r') as f:
    for line in f:
        coverage_data.append(int(line.split('\t')[2]))  # Get depth

# Create histogram
plt.hist(coverage_data, bins=50)
plt.title('Coverage Histogram for Chromosome 13')
plt.xlabel('Coverage Depth')
plt.ylabel('Number of Bases')
plt.show()