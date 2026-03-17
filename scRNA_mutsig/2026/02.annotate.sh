#!/bin/bash

#PBS -l select=2:ncpus=64:mem=128g
#PBS -l walltime=10:00:00
#PBS -P 11003581
#PBS -j oe
#PBS -N scomatic_annotation

set -euo pipefail
cd $PBS_O_WORKDIR

module load bcftools/1.15.1
module load r/4.2.0

###############################################
# PATHS
###############################################

BASE=/home/users/nus/ash.ps/scratch/mulitomics/example/analysis

INPUT=${BASE}/Step4_VariantCalling
ANNOT=${BASE}/annotated
POPVCF=${BASE}/pop-filter-vcfs
POPANN=${BASE}/pop-filter-multianno
MAF=${BASE}/mafs

mkdir -p $ANNOT $POPVCF $POPANN $MAF

ANNOVAR=/home/project/11003581/Tools/annovar
DB=${ANNOVAR}/humandb


###############################################
# TSV → VCF + ANNOVAR
###############################################

for tsv in ${INPUT}/*.calling.step1.tsv
do

sample=$(basename $tsv .calling.step1.tsv)

mkdir -p ${ANNOT}/${sample}

VCF=${ANNOT}/${sample}/${sample}.vcf

echo "Converting $sample"

{
echo "##fileformat=VCFv4.2"
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"

awk 'BEGIN{FS=OFS="\t"} !/^#/ {print $1,$2,".",$4,$5,".","PASS","."}' $tsv

} > $VCF


perl ${ANNOVAR}/table_annovar.pl \
$VCF \
$DB \
-buildver hg38 \
-out ${ANNOT}/${sample}/${sample} \
-remove \
-protocol refGene,cosmic70,gnomad_genome,avsnp151,dbnsfp30a,clinvar_20240611 \
-operation g,f,f,f,f,f \
-vcfinput -polish


done


###############################################
# POPULATION FILTERING
###############################################

find $ANNOT -name "*.hg38_multianno.vcf" | while read vcf
do

out=${POPVCF}/$(basename ${vcf%.vcf}_pop_filt.vcf)

bcftools view $vcf | awk '
BEGIN{FS=OFS="\t"}
/^#/ {print; next}

{
split($8,a,";")

gnomad="."

for(i in a){

if(a[i]~/gnomAD_genome_ALL/){
split(a[i],b,"=")
gnomad=b[2]
}

}

if(gnomad=="." || gnomad<=0.01)
print
}' > $out

done


###############################################
# CREATE MAF
###############################################

Rscript - <<EOF

library(maftools)

input_dir="$POPANN"
output_dir="$MAF"

dir.create(output_dir,showWarnings=FALSE)

files=list.files(input_dir,
pattern="multianno_pop_filt.txt",
full.names=TRUE)

maf=annovarToMaf(files,
refBuild="hg38",
ens2hugo=TRUE)

write.table(maf,
file=paste0(output_dir,"/combined.maf"),
sep="\t",
row.names=FALSE,
quote=FALSE)

EOF