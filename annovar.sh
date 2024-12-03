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


/home/users/nus/ash.ps/scratch/YS-analysis/vcf-to-maf/annovar/P01.hg19_multianno.vcf


setwd('/home/users/nus/ash.ps/scratch/YS-analysis/vcf-to-maf/')
library(maftools)

annovarToMaf('/home/users/nus/ash.ps/scratch/YS-analysis/vcf-to-maf/annovar/P01.hg19_multianno.txt', Center = NULL, refBuild = 'hg19', tsbCol = NULL,
table = 'refGene',basename = '/home/users/nus/ash.ps/scratch/YS-analysis/vcf-to-maf/P01.output_maf')

''' Additional databases
ljb26_all : whole-exome SIFT, PolyPhen2 HDIV, PolyPhen2 HVAR, LRT, MutationTaster, MutationAssessor, FATHMM, MetaSVM, MetaLR, VEST, CADD, GERP++, PhyloP and SiPhy scores from dbNSFP version 2.6
GTEx_v8_eQTL :	Expression Quantitative Trait Loci (eQTL) across tissues based on GTEx v8
cosmic70 : cosmic databse
exac03 : ExAC 65000 exome allele frequency data for ALL, AFR (African), AMR (Admixed American), EAS (East Asian), FIN (Finnish), NFE (Non-finnish European), OTH (other), SAS (South Asian)). version 0.3. Left normalization done.
gnomad_genome : 
1000g2015aug


/home/project/11003581/Tools/annovar/annotate_variation.pl -buildver hg38 -downdb -webfrom annovar gnomad_genome /home/project/11003581/Tools/annovar/humandb/ 
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar avsnp151 humandb/ 
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp47a humandb/

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
input_dir="/home/users/nus/ash.ps/scratch/YS-analysis/mutect2-vcfs"
output_dir="/home/users/nus/ash.ps/scratch/YS-analysis/annovar"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

for vcf_file in "$input_dir"/*mutect2.vcf; do
  # Extract the filename without the '.mutect2.vcf' extension
  samplename=$(basename "$vcf_file" .mutect2.vcf)

  # Print the sample name
  echo "Processing sample: $samplename"

  # Create a corresponding output directory
  mkdir -p "$output_dir/$samplename"

  # Run the table_annovar.pl script with the appropriate parameters
  perl /home/project/11003581/Tools/annovar/table_annovar.pl "$vcf_file" \
      /home/project/11003581/Tools/annovar/humandb/ \
      -buildver hg38 \
      -out "$output_dir/$samplename/$samplename" \
      -protocol refGene,GTEx_v8_eQTL,cosmic70,exac03,gnomad_genome,1000g2015aug \
      -operation gx,g,g,f,f,f \
      -remove -csvout -polish -nastring .

done