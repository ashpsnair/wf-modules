#!/bin/bash

#PBS -l select=1:ngpus=1:ncpus=64:mem=128g
#PBS -l walltime=15:00:00
#PBS -P 11003581
#PBS -N process-p3l-shallow

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

# Load necessary modules
module load gcc
module load python/3.12.1-gcc11
module load samtools

## Tools Directory
dorado_bin_path="/home/project/11003581/Tools/dorado-0.7.3-linux-x64/bin"
minimap2_path="/home/project/11003581/Tools/minimap2-2.28_x64-linux/"

fast5_dir="/home/project/11003581/Data/HC/ONT/P3292L/TJPROJ6/TGS/haiwai/haiwai/HW_ONT_qc/X401SC23084120-Z01-F001/data_release/X401SC23084120-Z01-F001/raw_data/P3292L/20231210_1504_6G_PAS94472_b08ff6e9/fast5_pass/"
output_dir="/home/project/11003581/Data/P3L-shallow"
samplename="p3l_shallow"
ref_fasta="/home/project/11003581/Ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna" 

# Convert each fast5 to its relative converted output
mkdir -p "$output_dir/pod5/"
pod5 convert fast5 "$fast5_dir"/*.fast5 --output "$output_dir/pod5/" --one-to-one "$fast5_dir/"

####### Running dorado
mkdir -p "$output_dir/dorado_output"
$dorado_bin_path/dorado basecaller hac "$output_dir/pod5/" \
--trim all \
--reference "$ref_fasta" \
> "$output_dir/dorado_output/${samplename}.bam"

##### Sorting & Indexing BAM files ######
mkdir -p "$output_dir/sorted_bams/"
samtools sort "$output_dir/dorado_output/${samplename}.bam" -o "$output_dir/sorted_bams/${samplename}_sorted.bam"
samtools index "$output_dir/sorted_bams/${samplename}_sorted.bam"

########## SV calling using Sniffles #####
mkdir -p "$output_dir/sniffles_out"
mkdir -p "$output_dir/sniffles_out/plots/"

sniffles -i "$output_dir/sorted_bams/${samplename}_sorted.bam" \
        -v "$output_dir/sniffles_out/${samplename}_sniffles.vcf" \
        --reference "$ref_fasta"

####### Sniffles Plot
python3 -m sniffles2_plot -i "$output_dir/sniffles_out/${samplename}_sniffles.vcf" -o "$output_dir/sniffles_out/plots/"

######### Running Mosdepth
module load miniforge3
conda activate /home/project/11003581/conda-envs/spectre

mkdir -p "$output_dir/mosdepth"
mosdepth -t 8 -x -b 1000 -Q 20 "$output_dir/mosdepth/${samplename}" "$output_dir/sorted_bams/${samplename}_sorted.bam"

########## R Script for Plotting Mosdepth ##########
Rscript -e "
setwd('$output_dir/mosdepth/')
library(maftools)
plotMosdepth(
  t_bed = '${samplename}.regions.bed.gz',
  n_bed = '/home/project/11003581/Data/Ash/P3L-lab-analysis/mosdepth/wt/wt.regions.bed.gz',
  segment = TRUE,
  sample_name = '$samplename'
)
"

#Running Spectre
mkdir -p "$output_dir/spectre"

# Ensure the path to spectre is correct
/home/project/11003581/conda-envs/spectre/bin/spectre CNVCaller \
  --coverage "$output_dir/mosdepth/${samplename}.regions.bed.gz" \
  --sample-id p3l \
  --only-chr chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY \
  --output-dir "$output_dir/spectre" \
  --reference "$ref_fasta"