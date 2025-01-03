#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  2 15:48:03 2025

@author: ash
"""

import os
import pandas as pd
import shutil
import pysam

# Define the directory containing the VCF files
vcf_dir = "/Users/ash/Downloads/MASH-vcfs/filter-vaf20-MASH/"

df= pd.read_excel("/Users/ash/Downloads/MASH-vcfs/patient-wise.csv.xlsx")

output_dir= "/Users/ash/Downloads/MASH-vcfs/pooled-patientwise-vaf20/"

# Loop through each row in the DataFrame
for index, row in df.iterrows():
    patient_id = row['patient_id']
    samplename = row['samplename']
    
    # Create a directory for the patient_id in the output directory if it doesn't exist
    patient_dir = os.path.join(output_dir, patient_id)
    os.makedirs(patient_dir, exist_ok=True)
    
    # Define the source file path
    source_file = os.path.join(vcf_dir, f"{samplename}-mutect2_filtered.vcf")
    
    # Check if the source file exists and copy it to the patient directory
    if os.path.isfile(source_file):
        shutil.copy(source_file, patient_dir)
        print(f"Copied {source_file} to {patient_dir}")
    else:
        print(f"File {source_file} does not exist.")

print("All files processed.")

##################################################################################
####### python script for segregating the samples
##################################################################################

import os
import shutil


# Create cohort directories
cohorts = ["NFL-F0", "FL-F0", "FL-F1-2", "FL-F3-4"]
for cohort in cohorts:
    os.makedirs(os.path.join(output_dir, cohort), exist_ok=True)

# Define patient ID folders and their corresponding cohort directories
patient_folders = {
    "NFL-F0": ["HEP152", "HEP229", "A009", "HEP549", "C010", "D009"],
    "FL-F0": ["HEP269", "A008", "B020", "PB_56025"],
    "FL-F1-2": ["HEP319", "C019", "F009"],
    "FL-F3-4": ["PB_15_28610", "HEP241", "HEP276", "HEP321", "B007", "B023", "D008", "D013"]
}


# Move patient ID folders to their respective cohort directories
for cohort, patients in patient_folders.items():
    for patient in patients:
        src_path = os.path.join(output_dir, patient)
        dest_path = os.path.join(output_dir, cohort)
        if os.path.exists(src_path):  # Check if the patient folder exists
            shutil.move(src_path, dest_path)
            print(f'Moved: {src_path} to {dest_path}')
        else:
            print(f'Folder not found: {src_path}')



##################################################################################
####### shell script for bcftools isec
##################################################################################

base_dir="/Users/ash/Downloads/MASH-vcfs/pooled-patientwise-vaf20/"

# Loop through each subfolder
for dir in "$base_dir"/*/*/; do
    cd "$dir" || continue  # Change to the subfolder, skip if it fails

    echo "Processing directory: $dir"
    
    # Get a list of VCF files and compress them
    vcf_files=(*.vcf)
    
    # Compress all VCF files
    bgzip "${vcf_files[@]}"

    # Index the compressed VCF files
    for vcf in *.vcf.gz; do
        bcftools index "$vcf"
    done
    
    echo "compressed and indexed"
    
    # Get a list of VCF files (only .vcf.gz)
    vcf_files=(*.vcf.gz)
    
    # Get the number of compressed VCF files
    num_files=${#vcf_files[@]}
    
    # Perform bcftools isec using the compressed VCF files
    output_prefix="uniq-"
    
    # Create output directory if it doesn't exist
    mkdir -p output

    # Run bcftools isec and save outputs in the 'output' directory
    echo "Running bcftools isec with ${num_files} files..."
    
   # Use bcftools isec on the .vcf.gz files directly
    bcftools isec -n-"$((num_files - 1))" -p output "${vcf_files[@]}"
    
    bcftools isec -n="$num_files" -p common "${vcf_files[@]}"
    
    # Rename output files with prefix uniq- based on original VCF filenames
    for i in "${!vcf_files[@]}"; do
        original_name="${vcf_files[i]}"
        base_name="${original_name%.vcf.gz}"  # Remove .vcf.gz extension for base name
        
        # Ensure that the expected output file exists before renaming
        if [ -e "output/000${i}.vcf" ]; then
            mv "output/000${i}.vcf" "output/${output_prefix}${base_name}.vcf"  # Rename with prefix
            echo "Renamed output/000${i}.vcf to output/${output_prefix}${base_name}.vcf"
        else
            echo "Expected output file output/000${i}.vcf does not exist."
        fi        
    done

    cd .. # Go back to the parent directory

done

##################################################################################
#################### shell script to run from terminal to Count variants each files 
##################################################################################

### Counting the filtered vcf files
for file in */*/*.vcf.gz; do

    # Count the number of variants in the current file
    variant_count=$(bcftools view -H "$file" | wc -l)

    # Output the result to both the terminal and the output file
    echo "$file:$variant_count" | tee -a "raw-counts.txt"
done

### Counting the common vcf files
for file in */*/common/*.vcf; do

    # Count the number of variants in the current file
    variant_count=$(bcftools view -H "$file" | wc -l)

    # Output the result to both the terminal and the output file
    echo "$file:$variant_count" | tee -a "common-counts.txt"
done

### Counting the variants after substracting commons
for file in */*/output/*.vcf; do

    # Count the number of variants in the current file
    variant_count=$(bcftools view -H "$file" | wc -l)

    # Output the result to both the terminal and the output file
    echo "$file:$variant_count" | tee -a "filtered-counts.txt"
done

##################################################################################
####### python script for copying sigprofile inputs
##################################################################################

import os
import shutil
import glob

# Define the base directory containing patient-wise VCF files
base_dir = "/Users/ash/Downloads/MASH-vcfs/patient-wise-vaf30/"

# Loop through each subfolder in the base directory
for dir_name in os.listdir(base_dir):
    dir_path = os.path.join(base_dir, dir_name)
    
    # Check if it's a directory
    if os.path.isdir(dir_path):
        # Define the path for the new sigprofile directory
        sigprofile_dir = os.path.join(dir_path, 'sigprofile')
        
        # Check if the directory already exists
        if not os.path.exists(sigprofile_dir):
            os.makedirs(sigprofile_dir)  # Create the sigprofile directory
            print(f"Created directory: {sigprofile_dir}")

        else:
            print(f"Directory already exists: {sigprofile_dir}")
        
        # Define the pattern for .vcf files in the output subdirectory
        vcf_pattern = os.path.join(dir_path, '*', 'output', '*.vcf')
        vcf_files = glob.glob(vcf_pattern)  # Find all matching .vcf files

        # Copy each .vcf file to the sigprofile directory
        for vcf_file in vcf_files:
            shutil.copy(vcf_file, sigprofile_dir)  # Copy file
            print(f"Copied {vcf_file} to {sigprofile_dir}")



##################################################################################
####### python script for sigprofile generator
##################################################################################

from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen


matrices = matGen.SigProfilerMatrixGeneratorFunc("NFL-F0", "GRCh38", "/Users/ash/Downloads/MASH-vcfs/patient-wise-vaf30/NFL-F0/sigprofile/", plot=True, exome=False, bed_file=None, chrom_based=False, tsb_stat=False, seqInfo=False, cushion=100)


matrices = matGen.SigProfilerMatrixGeneratorFunc("FL-F0", "GRCh38", "/Users/ash/Downloads/MASH-vcfs/patient-wise-vaf30/FL-F0/sigprofile/", plot=True, exome=False, bed_file=None, chrom_based=False, tsb_stat=False, seqInfo=False, cushion=100)


matrices = matGen.SigProfilerMatrixGeneratorFunc("FL-F1-2", "GRCh38", "/Users/ash/Downloads/MASH-vcfs/patient-wise-vaf30/FL-F1-2/sigprofile/", plot=True, exome=False, bed_file=None, chrom_based=False, tsb_stat=False, seqInfo=False, cushion=100)


matrices = matGen.SigProfilerMatrixGeneratorFunc("FL-F3-4", "GRCh38", "/Users/ash/Downloads/MASH-vcfs/patient-wise-vaf30/FL-F3-4/sigprofile/", plot=True, exome=False, bed_file=None, chrom_based=False, tsb_stat=False, seqInfo=False, cushion=100)



##################################################################################
####### R script for sigminer
##################################################################################

library(sigminer)
library(tibble)

matrix <- read.delim("~/Downloads/MASH-vcfs/patient-wise-vaf30/FL-F0/sigprofile/output/SBS/FL-F0.SBS96.all", sep = "\t")
output <- "~/Downloads/MASH-vcfs/patient-wise-vaf30/FL-F0/sigprofile/sigminer/"

dir.create(output)

# Formats the matrix for sigminer if necessary
matrix <- as.matrix(setNames(data.frame(t(matrix[,-1])), matrix[,1]))

# Signature extraction using Bayesian NMF from K=1 to K=10
model <- sig_unify_extract(
  matrix,
  range = 1:10,
  approach = "bayes_nmf",
  nrun = 10
)


# Save H and W matrices
signatures <- as.data.frame(model$Signature.norm)
signatures <- cbind(MutationsType = rownames(signatures), signatures)
write.table(signatures, paste0(output, "denovo_sig.tsv"), sep = "\t", quote = FALSE, row.name = FALSE)
exposure <- tibble::rownames_to_column(as.data.frame(t(model$Exposure)), "Sample")
write.table(exposure, paste0(output, "denovo_exposure.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Save run logs
write.table(model$Raw$summary_run, paste0(output, "run.log"), sep = "\t", quote = FALSE, row.name = FALSE)

# [optional] Visualize de novo signatures similarity to COSMIC
library(pheatmap)
library(repr)

sim <- get_sig_similarity(model, sig_db = "latest_SBS_GRCh38")

# Save heatmap as SVG without defining aspect ratio
png(paste0(output, "similarity_heatmap.png"), width = 1200, height = 350)   

pheatmap(sim$similarity,
         color = colorRampPalette(c("#348ABD", "white", "#E24A33"))(70),
         border_color = FALSE)

dev.off()  # Close the png device



##################################################################################
####### R script for sigprofile assigner
##################################################################################

from SigProfilerAssignment import Analyzer as Analyze

signatures = '~/Downloads/MASH-vcfs/patient-wise-vaf30/FL-F0/sigprofile/sigminer/denovo_sig.tsv'
activities = '~/Downloads/MASH-vcfs/patient-wise-vaf30/FL-F0/sigprofile/sigminer/denovo_exposure.tsv'
samples = '~/Downloads/MASH-vcfs/patient-wise-vaf30/FL-F0/sigprofile/output/SBS/FL-F0.SBS96.all'
output = '/Users/ash/Downloads/MASH-vcfs/patient-wise-vaf30/FL-F0/sigprofile/sigminer/'

Analyze.decompose_fit(samples=samples,output=output,signatures=signatures,genome_build="GRCh38",signature_database='/Users/ash/opt/anaconda3/lib/python3.12/site-packages/SigProfilerAssignment/data/Reference_Signatures/GRCh38/COSMIC_v3.3_SBS_GRCh38.txt')



from SigProfilerAssignment import Analyzer as Analyze

signatures = '~/Downloads/MASH-vcfs/patient-wise-vaf30/FL-F1-2/sigprofile/sigminer/denovo_sig.tsv'
activities = '~/Downloads/MASH-vcfs/patient-wise-vaf30/FL-F1-2/sigprofile/sigminer/denovo_exposure.tsv'
samples = '~/Downloads/MASH-vcfs/patient-wise-vaf30/FL-F1-2/sigprofile/output/SBS/FL-F1-2.SBS96.all'
output = '/Users/ash/Downloads/MASH-vcfs/patient-wise-vaf30/FL-F1-2/sigprofile/sigminer/'

Analyze.decompose_fit(samples=samples,output=output,signatures=signatures,genome_build="GRCh38",signature_database='/Users/ash/opt/anaconda3/lib/python3.12/site-packages/SigProfilerAssignment/data/Reference_Signatures/GRCh38/COSMIC_v3.3_SBS_GRCh38.txt')


from SigProfilerAssignment import Analyzer as Analyze

signatures = '~/Downloads/MASH-vcfs/patient-wise-vaf30/FL-F3-4/sigprofile/sigminer/denovo_sig.tsv'
activities = '~/Downloads/MASH-vcfs/patient-wise-vaf30/FL-F3-4/sigprofile/sigminer/denovo_exposure.tsv'
samples = '~/Downloads/MASH-vcfs/patient-wise-vaf30/FL-F3-4/sigprofile/output/SBS/FL-F3-4.SBS96.all'
output = '/Users/ash/Downloads/MASH-vcfs/patient-wise-vaf30/FL-F3-4/sigprofile/sigminer/'

Analyze.decompose_fit(samples=samples,output=output,signatures=signatures,genome_build="GRCh38",signature_database='/Users/ash/opt/anaconda3/lib/python3.12/site-packages/SigProfilerAssignment/data/Reference_Signatures/GRCh38/COSMIC_v3.3_SBS_GRCh38.txt')



from SigProfilerAssignment import Analyzer as Analyze

signatures = '~/Downloads/MASH-vcfs/patient-wise-vaf30/NFL-F0/sigprofile/sigminer/denovo_sig.tsv'
activities = '~/Downloads/MASH-vcfs/patient-wise-vaf30/NFL-F0/sigprofile/sigminer/denovo_exposure.tsv'
samples = '~/Downloads/MASH-vcfs/patient-wise-vaf30/NFL-F0/sigprofile/output/SBS/NFL-F0.SBS96.all'
output = '/Users/ash/Downloads/MASH-vcfs/patient-wise-vaf30/NFL-F0/sigprofile/sigminer/'

Analyze.decompose_fit(samples=samples,output=output,signatures=signatures,genome_build="GRCh38",signature_database='/Users/ash/opt/anaconda3/lib/python3.12/site-packages/SigProfilerAssignment/data/Reference_Signatures/GRCh38/COSMIC_v3.3_SBS_GRCh38.txt')















