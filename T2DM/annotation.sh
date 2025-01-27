#!/bin/bash

#PBS -l select=1:ncpus=128:mem=128g
#PBS -l walltime=06:00:00
#PBS -P 11003581
#PBS -N t2dm-annotation
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

# Define the input and output directories
input_dir="/home/users/nus/ash.ps/scratch/T2DM/PD-analysis/filtered-vcfs/"
output_dir="/home/users/nus/ash.ps/scratch/T2DM/PD-analysis/annotated/"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

for vcf_file in "$input_dir"/*_filtered.vcf; do
  # Extract the filename without the '.mutect2.vcf' extension
  samplename=$(basename "$vcf_file" _filtered.vcf)

  # Print the sample name
  echo "Processing sample: $samplename"

  # Create a corresponding output directory
  mkdir -p "$output_dir/$samplename"

  # Run the table_annovar.pl script with the appropriate parameters
  perl /home/project/11003581/Tools/annovar/table_annovar.pl "$vcf_file" \
      /home/project/11003581/Tools/annovar/humandb/ \
      -buildver hg38 \
      -out "$output_dir/$samplename/$samplename" \
      -protocol refGene,cosmic70,exac03,gnomad_genome,ALL.sites.2015_08,SAS.sites.2015_08,EAS.sites.2015_08,esp6500siv2_all,avsnp151,icgc28,dbnsfp30a,intervar_20180118,clinvar_20240611 \
      -operation g,f,f,f,f,f,f,f,f,f,f,f,f \
      -remove -vcfinput -polish -nastring .

done