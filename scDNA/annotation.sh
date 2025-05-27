#!/bin/bash

#PBS -l select=3:ncpus=64:mem=256g
#PBS -l walltime=06:00:00
#PBS -P 11003581
#PBS -N scdna-annotate
#PBS -j oe
#PBS -M ash.ps@nus.edu.sg

cd $PBS_O_WORKDIR
module load bcftools/1.15.1

base_dir="/home/users/nus/ash.ps/scratch/scDNA/analysis"
annot_dir="${base_dir}/annotation"

mkdir -p "${annot_dir}"/{annotated}

###############################################
### annotation
###############################################

for cat in 2204 3401; do
  input_dir="${base_dir}/filtered-vcfs/$cat"
  output_dir="${annot_dir}/annotated/$cat"
  mkdir -p "$output_dir"

  for vcf_file in "$input_dir"/*.vcf.gz; do
    samplename=$(basename "$vcf_file" .vcf.gz)
    echo "Annotating: $samplename"
    mkdir -p "$output_dir/$samplename"

    perl /home/project/11003581/Tools/annovar/table_annovar.pl "$vcf_file" \
      /home/project/11003581/Tools/annovar/humandb/ \
      -buildver hg38 \
      -out "$output_dir/$samplename/$samplename" \
      -protocol refGene,cosmic70,exac03,gnomad_genome,ALL.sites.2015_08,SAS.sites.2015_08,EAS.sites.2015_08,esp6500siv2_all,avsnp151,icgc28,dbnsfp30a,intervar_20180118,clinvar_20240611 \
      -operation g,f,f,f,f,f,f,f,f,f,f,f,f \
      -remove -vcfinput -polish -nastring .
  done
done

