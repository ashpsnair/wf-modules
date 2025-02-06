#!/bin/bash

#PBS -l select=1:ncpus=32
#PBS -l walltime=6:00:00
#PBS -P 11003581
#PBS -N annovar
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR


perl /home/project/11003581/Tools/annovar/table_annovar.pl /home/users/nus/ash.ps/scratch/YS-analysis/mutect2-vcfs/T01_vs_N01.mutect2.filtered.vcf \
    /home/project/11003581/Tools/annovar/humandb/ \
    -buildver hg38 \
    -out /home/users/nus/ash.ps/scratch/YS-analysis/vcf-to-maf/P01 \
    -protocol refGeneWithVer \
    -operation g \
    -remove -polish -vcfinput -nastring .



''' Additional databases
ljb26_all : whole-exome SIFT, PolyPhen2 HDIV, PolyPhen2 HVAR, LRT, MutationTaster, MutationAssessor, FATHMM, MetaSVM, MetaLR, VEST, CADD, GERP++, PhyloP and SiPhy scores from dbNSFP version 2.6
cosmic70 : cosmic databse
exac03 : ExAC 65000 exome allele frequency data for ALL, AFR (African), AMR (Admixed American), EAS (East Asian), FIN (Finnish), NFE (Non-finnish European), OTH (other), SAS (South Asian)). version 0.3. Left normalization done.
gnomad_genome : 
1000g2015aug


/home/project/11003581/Tools/annovar/annotate_variation.pl -buildver hg38 -downdb -webfrom annovar gnomad_genome /home/project/11003581/Tools/annovar/humandb/ 
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar avsnp151 humandb/ 
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp47a humandb/

annotate_variation.pl -buildver hg38 -downdb -webfrom annovar icgc21 /home/project/11003581/Tools/annovar/humandb/

'''

### Automated annovar script

#!/bin/bash

#PBS -l select=1:ncpus=16
#PBS -l walltime=00:30:00
#PBS -P 11003581
#PBS -N auto-annovar-hg38
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

# Define the input and output directories
input_dir="/home/users/nus/ash.ps/scratch/MASH-analysis/vcf-inputs/"
output_dir="/home/users/nus/ash.ps/scratch/MASH-analysis/annovar"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

for vcf_file in "$input_dir"/*mutect2.vcf; do
  # Extract the filename without the '.mutect2.vcf' extension
  samplename=$(basename "$vcf_file" mutect2.vcf)

  # Print the sample name
  echo "Processing sample: $samplename"

  # Create a corresponding output directory
  mkdir -p "$output_dir/$samplename"

  # Run the table_annovar.pl script with the appropriate parameters
  perl /home/project/11003581/Tools/annovar/table_annovar.pl "$vcf_file" \
      /home/project/11003581/Tools/annovar/humandb/ \
      -buildver hg38 \
      -out "$output_dir/$samplename/$samplename" \
      -protocol refGene,cosmic70,exac03,gnomad_genome,1000g2015aug \
      -operation gx,g,g,f,f,f \
      -remove -vcfinput -csvout -polish -nastring .

done



########### annovar to maf
#!/bin/bash

#PBS -l select=1:ncpus=16:mem=128g
#PBS -l walltime=03:30:00
#PBS -P 11003581
#PBS -N create-maf
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load r

# Run the R command directly within the bash script using Rscript
Rscript -e "
library(maftools)

# Define the input and output directories
input_dir='/home/users/nus/ash.ps/scratch/MASH-analysis/annovar/'
output_dir='/home/users/nus/ash.ps/scratch/MASH-analysis/'

# Get the list of .hg38_multianno.txt files in the specified location (including subfolders)
annovar_outputs <- list.files(path = input_dir, pattern = '\\\\.hg38_multianno\\\\.txt$', recursive = TRUE, full.names = TRUE)

# Run the annovarToMaf function
MASH_maf <- annovarToMaf(
  annovar_outputs,
  Center = NULL,
  refBuild = 'hg38',
  tsbCol = NULL,
  table = 'refGene',
  ens2hugo = TRUE,
  basename = NULL,
  sep = '\\t',
  MAFobj = FALSE,
  sampleAnno = NULL
)

# Save the output (modify as needed)
write.table(MASH_maf, file=paste0(output_dir, '/MASH_maf_output.txt'), sep='\\t', row.names=FALSE, quote=FALSE)
"

echo "Processing complete. Output saved to $output_dir/MASH_maf_output.txt"


############## Split and left align

bcftools norm "$input_dir"/*mutect2.vcf.gz --fasta-ref /home/project/11003581/Ref/goldenPaths/hg38.fa.gz -m -O z --threads 128


## Automated annovar script- rerun 2

#!/bin/bash

#PBS -l select=1:ncpus=128:mem=128g
#PBS -l walltime=10:00:00
#PBS -P 11003581
#PBS -N auto-annovar-YS
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

# Define the input and output directories
input_dir="/home/users/nus/ash.ps/scratch/YS-analysis/filtered-vcfs/"
output_dir="/home/users/nus/ash.ps/scratch/YS-analysis/annotated/"

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


/home/project/11003581/Tools/annovar/humandb/hg38_1000g2015aug.txt
