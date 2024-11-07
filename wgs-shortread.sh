
ref_fasta=/home/project/11003581/Ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
output_dir=/home/project/11003581/Data/Ash/wgs-tcga-pipeline/align-out
gnomad_vcf=af-only-gnomad.vcf.gz

#Matched tumor normal

gatk Mutect2 \
    -R $ref_fasta \
    -I tumor.bam \
    -I normal.bam \
    -normal normal_sample_name \
    --germline-resource $gnomad_vcf \
    -O $output_dir/somatic.vcf.gz