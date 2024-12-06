####################################################
### Script for concatenating fq files into a folder
####################################################

#!/bin/bash

#PBS -l select=1
#PBS -l walltime=2:00:00
#PBS -P 11003581
#PBS -N concatenating-fq
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

#!/bin/bash

#PBS -l select=1
#PBS -l walltime=2:00:00
#PBS -P 11003581
#PBS -N concatenating-fq
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

# Define the paths
rawdata="/home/project/11003581/Data/YS1/01.RawData"
analysis_folder="/home/project/11003581/Data/YS-analysis/"

# Create the analysis folder if it doesn't exist
mkdir -p "$analysis_folder/raw_data"

# Log file for output
log_file="$PBS_O_WORKDIR/cat-fq.log"
size_info_file="$PBS_O_WORKDIR/size-info.txt"

# Clear previous log files
> "$log_file"
> "$size_info_file"

# Loop through each subdirectory in rawdata
for dir in "$rawdata"/*; do
    if [ -d "$dir" ]; then
        # Extract the folder name
        folder_name=$(basename "$dir")
        
        # List all _1.fq.gz files in the current directory
        fq1_files=("$dir"/*_1.fq.gz)
        
        # Check if there are any _1.fq.gz files
        if [ ${#fq1_files[@]} -gt 0 ]; then
            # Calculate size of all _1.fq.gz files before concatenation
            total_size_before=$(du -sh "${fq1_files[@]}" | awk '{sum += $1} END {print sum}')
            echo "Folder: $folder_name" >> "$log_file"
            echo "Size of all *_1.fq.gz files: $total_size_before" >> "$size_info_file"
            
            # Concatenate all _1.fq.gz files into one
            cat "${fq1_files[@]}" > "$analysis_folder/raw_data/${folder_name}_1.fq.gz"
            
            # Get size of concatenated file
            size_after=$(du -sh "$analysis_folder/raw_data/${folder_name}_1.fq.gz" | awk '{print $1}')
            echo "Size of ${folder_name}_1.fq.gz after concatenation: $size_after" >> "$size_info_file"
        fi
        
        # List all _2.fq.gz files in the current directory (if needed)
        fq2_files=("$dir"/*_2.fq.gz)
        
        # Check if there are any _2.fq.gz files
        if [ ${#fq2_files[@]} -gt 0 ]; then
            # Calculate size of all _2.fq.gz files before concatenation
            total_size_before_2=$(du -sh "${fq2_files[@]}" | awk '{sum += $1} END {print sum}')
            echo "Size of all *_2.fq.gz files: $total_size_before_2" >> "$size_info_file"
            
            # Concatenate all _2.fq.gz files into one
            cat "${fq2_files[@]}" > "$analysis_folder/raw_data/${folder_name}_2.fq.gz"
            
            # Get size of concatenated file
            size_after_2=$(du -sh "$analysis_folder/raw_data/${folder_name}_2.fq.gz" | awk '{print $1}')
            echo "Size of ${folder_name}_2.fq.gz after concatenation: $size_after_2" >> "$size_info_file"
        fi
    fi
done

# Indicate completion in log file
echo "Concatenation completed." >> "$log_file"

####################################################
# Creating samplesheet.csv in analysis folder
####################################################

import os
import csv

# Define the directory containing .fq.gz files
data_directory = '/home/project/11003581/Data/YS-analysis/raw_data/'

# Define the sex mapping dictionary
sex_mapping = {
    'T01': 'XY', 'N01': 'XY',
    'T02': 'XY', 'N02': 'XY',
    'T03': 'XX', 'N03': 'XX',
    'T04': 'XY', 'N04': 'XY',
    'T05': 'XY', 'N05': 'XY',
    'T06': 'XX', 'N06': 'XX',
    'T07': 'XY', 'N07': 'XY',
    'T08': 'XY', 'N08': 'XY',
    'T09': 'XY', 'N09': 'XY',
    'T10': 'XY', 'N10': 'XY',
    'T11': 'XY', 'N11': 'XY',
    'T12': 'XY', 'N12': 'XY',
    'T13': 'XY', 'N13': 'XY',
    'T14': 'XY', 'N14': 'XY'
}

# Prepare to collect data for the CSV
data_rows = []
unique_entries = set()  # To track unique entries

# Loop through all files in the data directory
for filename in os.listdir(data_directory):
    if filename.endswith('_1.fq.gz') or filename.endswith('_2.fq.gz'):
        # Extract patient identifier (e.g., N01 or T01)
        patient_id = filename.split('_')[0]
        sample_type = patient_id[0]  # N or T
        patient_number = patient_id[1:]  # 01, 02, etc.
        
        # Determine sex and status
        sex = sex_mapping.get(patient_id, '')
        status = 1 if sample_type == "T" else 0
        
        # Lane is fixed as lane_1 for all entries
        lane = "lane_1"
        
        # Determine fastq file names based on the suffix
        fastq_1 = os.path.join(data_directory, f"{patient_id}_1.fq.gz")
        fastq_2 = os.path.join(data_directory, f"{patient_id}_2.fq.gz")
        
        # Create a unique entry key based on patient ID and sample type
        entry_key = (patient_number, fastq_1, fastq_2)  # Use patient_number instead of patient_id
        
        # Check if this entry is unique
        if entry_key not in unique_entries:
            unique_entries.add(entry_key)
            patient_formatted = f"patient{patient_number}"  # Format patient as patient01, patient02, etc.
            data_rows.append([patient_formatted, sex, status, patient_id, lane, fastq_1, fastq_2])

# Sort data rows by patient column (first column)
data_rows.sort(key=lambda x: x[0])  # Sort by formatted patient name

# Define the output CSV file path
output_csv_path = os.path.join(data_directory, "samplesheet.csv")

# Write to CSV file
with open(output_csv_path, mode='w', newline='') as csv_file:
    writer = csv.writer(csv_file)
    
    # Write header
    writer.writerow(['patient', 'sex', 'status', 'sample', 'lane', 'fastq_1', 'fastq_2'])
    
    # Write data rows
    writer.writerows(data_rows)

# Print number of rows in the samplesheet.csv
print(f"CSV file '{output_csv_path}' created successfully with {len(data_rows)} rows.")

####################################################
# Running sarek
####################################################

#!/bin/bash

#PBS -l select=1:ncpus=64:mem=256g
#PBS -l walltime=20:00:00
#PBS -P 11003581
#PBS -N YS-8-14-sarek
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR


module load java/17.0.6-jdk
module load singularity

/home/project/11003581/Tools/nextflow run nf-core/sarek \
   -profile singularity \
   --input /home/project/11003581/Data/YS-analysis/batch2_8_14/samplesheet.csv \
   --outdir /home/project/11003581/Data/YS-analysis/batch2_8_14/ \
   --tools mutect2,ascat,manta,snpeff,vep,merge,haplotypecaller,msisensorpro \
   -resume

####################################################
# Running tumor only for all normals
####################################################

#!/bin/bash

#PBS -l select=3:ncpus=64:mem=128g
#PBS -l walltime=24:00:00
#PBS -P 11003581
#PBS -N YS-TO
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR


module load java/17.0.6-jdk
module load singularity

/home/project/11003581/Tools/nextflow run nf-core/sarek -r 3.4.4 \
   -profile singularity \
   --input /home/users/nus/ash.ps/scratch/YS-tumor-only/YS8-14/samplesheet.csv \
   --outdir /home/users/nus/ash.ps/scratch/YS-tumor-only/YS8-14/ \
   --tools mutect2,ascat,manta,snpeff,vep,msisensorpro \
   --pon /home/project/11003581/Ref/pons/somatic-hg38_1000g_pon.hg38.vcf.gz \
   --pon_tbi /home/project/11003581/Ref/pons/somatic-hg38_1000g_pon.hg38.vcf.gz.tbi


####################################################
# Annotating the files
####################################################
#!/bin/bash

#PBS -l select=2:ncpus=64:mem=256g
#PBS -l walltime=18:00:00
#PBS -P 11003581
#PBS -N YS-annotation-except 13
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR


module load java/17.0.6-jdk
module load singularity

/home/project/11003581/Tools/nextflow run nf-core/sarek -r 3.4.4 \
   -profile singularity \
   --step annotate \
   --input /home/users/nus/ash.ps/scratch/YS-rerun-nf/YS5-8/samplesheet.csv \
   --outdir /home/users/nus/ash.ps/scratch/YS-rerun-nf/YS5-8/ \
   --tools snpeff,vep,merge


patient,sample,variantcaller,vcf
patient05,T05_vs_N05,mutect2,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS5-8//variant_calling/mutect2/T05_vs_N05/T05_vs_N05.mutect2.filtered.vcf.gz
patient06,T06_vs_N06,mutect2,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS5-8//variant_calling/mutect2/T06_vs_N06/T06_vs_N06.mutect2.filtered.vcf.gz
patient07,T07_vs_N07,mutect2,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS5-8//variant_calling/mutect2/T07_vs_N07/T07_vs_N07.mutect2.filtered.vcf.gz
patient08,T08_vs_N08,mutect2,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS5-8//variant_calling/mutect2/T08_vs_N08/T08_vs_N08.mutect2.filtered.vcf.gz
patient09,T09_vs_N09,mutect2,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS9-11//variant_calling/mutect2/T09_vs_N09/T09_vs_N09.mutect2.filtered.vcf.gz
patient10,T10_vs_N10,mutect2,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS9-11//variant_calling/mutect2/T10_vs_N10/T10_vs_N10.mutect2.filtered.vcf.gz
patient11,T11_vs_N11,mutect2,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS9-11//variant_calling/mutect2/T11_vs_N11/T11_vs_N11.mutect2.filtered.vcf.gz
patient12,T12_vs_N12,mutect2,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS12-14//variant_calling/mutect2/T12_vs_N12/T12_vs_N12.mutect2.filtered.vcf.gz
patient14,T14_vs_N14,mutect2,/home/users/nus/ash.ps/scratch/YS-rerun-nf/YS12-14//variant_calling/mutect2/T14_vs_N14/T14_vs_N14.mutect2.filtered.vcf.gz

















