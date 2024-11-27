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

#samtools depth /home/project/11003581/Data/Ash/P3L-lab-analysis/sorted_bams/p3292l_hac_sorted.bam > /home/project/11003581/Data/Ash/P3L-lab-analysis/coverage/p3l_coverage.txt

samtools coverage /home/project/11003581/Data/Ash/P3L-lab-analysis/coverage/chr12_13_reads.bam -o /home/project/11003581/Data/Ash/P3L-lab-analysis/coverage/cov-out



#!/bin/bash

#PBS -l select=1:ncpus=64:mem=256g:ngpus=2
#PBS -l walltime=1:00:00
#PBS -P 11003581
#PBS -N plot-depth
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

# Load necessary modules
module load gcc
module load python/3.12.1-gcc11

# Create a Python script for plotting
cat << 'EOF' > plot_depth.py

import pandas as pd
import matplotlib.pyplot as plt

# Load the depth data
depth_data = pd.read_csv('/home/project/11003581/Data/Ash/P3L-lab-analysis/coverage/chr13_depth.txt', 
                          sep='\t', header=None, 
                          names=['Chromosome', 'Position', 'Depth'])

# Filter for chromosome 13
chr13_data = depth_data[depth_data['Chromosome'] == 'chr13']

# Create a new figure
plt.figure(figsize=(12, 6))

# Plotting depth against genomic position for chromosome 13
plt.plot(chr13_data['Position'], chr13_data['Depth'], color='blue')

# Adding labels and title
plt.xlabel('Genomic Position (bp)')
plt.ylabel('Depth')
plt.title('Depth of Coverage for Chromosome 13')
plt.grid()

# Show the plot
plt.show()

EOF

# Run the Python script
python plot_depth.py
