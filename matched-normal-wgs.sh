#!/bin/bash

#PBS -l select=1:ngpus=2:ncpus=64:mpiprocs=2 
#PBS -l walltime=10:00:00
#PBS -P 11003581
#PBS -N process-p3l-sup
#PBS -j oe

####### mandatory fields ########
wt_samplename=WT_MCF10A
case_samplename=p3292l
dorado_model=sup

#filepaths
output_dir="/home/project/11003581/Data/Ash/P3292L/"
ref_fasta="/home/project/11003581/Ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
wt_pod5_dir="/home/users/nus/ash.ps/scratch/MCF10A-lab/pod5/WT/"
case_pod5_dir="/home/users/nus/ash.ps/scratch/MCF10A-lab/pod5/P3L/"

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load python
module load bcftools
module load samtools

##Tools Directory
dorado_bin_path="/home/project/11003581/Tools/dorado-0.7.3-linux-x64/bin"
minimap2_path="/home/project/11003581/Tools/minimap2-2.28_x64-linux/"

####### Running dorado
mkdir $output_dir/dorado_output

#running for WT
$dorado_bin_path/dorado basecaller $dorado_model $wt_pod5_dir \
--trim all \
--reference $ref_fasta \
> $output_dir/dorado_output/${wt_samplename}_${dorado_model}.bam &

#running for case
$dorado_bin_path/dorado basecaller $dorado_model $case_pod5_dir \
--trim all \
--reference $ref_fasta \
> $output_dir/dorado_output/${case_samplename}_${dorado_model}.bam &

wait

$dorado_bin_path/dorado summary $output_dir/dorado_output/${wt_samplename}_${dorado_model}.bam > $output_dir/dorado_output/{$wt_samplename}_dorado_summary.txt &
$dorado_bin_path/dorado summary $output_dir/dorado_output/${case_samplename}_${dorado_model}.bam > $output_dir/dorado_output/{$case_samplename}_dorado_summary.txt &

wait

echo "dorado done"

##### Sorting & Indexing BAM files ######
mkdir -p $output_dir/sorted_bams

samtools sort $output_dir/dorado_output/${wt_samplename}_${dorado_model}.bam -o $output_dir/sorted_bams/${wt_samplename}_${dorado_model}_sorted.bam &
samtools sort $output_dir/dorado_output/${case_samplename}_${dorado_model}.bam -o $output_dir/sorted_bams/${case_samplename}_${dorado_model}_sorted.bam &

wait 

samtools index $output_dir/sorted_bams/${wt_samplename}_${dorado_model}_sorted.bam &
samtools index $output_dir/sorted_bams${case_samplename}_${dorado_model}_sorted.bam &

echo "sorting indexing done"

########## SV calling using Sniffles #####
mkdir -p $output_dir/sniffles_out

sniffles -i $output_dir/sorted_bams/${wt_samplename}_${dorado_model}_sorted.bam \
        -v $output_dir/sniffles_out/${wt_samplename}_${dorado_model}_sniffles.vcf \
        --reference $ref_fasta &

sniffles -i $output_dir/sorted_bams${case_samplename}_${dorado_model}_sorted.bam \
        -v $output_dir/sniffles_out/${case_samplename}_${dorado_model}_sniffles.vcf \
        --reference $ref_fasta &

wait 

echo "sniffles done"

####### Sniffles Plot
/home/users/nus/ash.ps/anaconda3/bin/python3.8 -m sniffles2_plot -i $output_dir/sniffles_out/${wt_samplename}_${dorado_model}_sniffles.vcf -o $output_dir/sniffles_out/plot &
/home/users/nus/ash.ps/anaconda3/bin/python3.8 -m sniffles2_plot -i $output_dir/sniffles_out/${case_samplename}_${dorado_model}_sniffles.vcf -o $output_dir/sniffles_out/plot &

wait

###### intersecting sniffles
mkdir -p $output_dir/sniffles_out/sniffles_intersect/

bcftools view -Oz -o $output_dir/sniffles_out/sniffles_intersect/${wt_samplename}_${dorado_model}_sniffles.vcf.gz $output_dir/sniffles_out/${wt_samplename}_${dorado_model}_sniffles.vcf &
bcftools view -Oz -o $output_dir/sniffles_out/sniffles_intersect/${case_samplename}_${dorado_model}_sniffles.vcf.gz $output_dir/sniffles_out/${case_samplename}_${dorado_model}_sniffles.vcf &

wait

bcftools index $output_dir/sniffles_out/sniffles_intersect/${wt_samplename}_${dorado_model}_sniffles.vcf.gz &
bcftools index $output_dir/sniffles_out/sniffles_intersect/${case_samplename}_${dorado_model}_sniffles.vcf.gz &

bcftools isec -p $output_dir/sniffles_out/sniffles_intersect/intersect_output -n=1 \
$output_dir/sniffles_out/sniffles_intersect/${wt_samplename}_${dorado_model}_sniffles.vcf.gz \
$output_dir/sniffles_out/sniffles_intersect/${case_samplename}_${dorado_model}_sniffles.vcf.gz

echo "sniffles intersect done"

####### Running CuteSV
mkdir -p $output_dir/cutesv_output  # Use -p to avoid error if the directory already exists

# Run the first cuteSV command in the background
cuteSV $output_dir/sorted_bams/${wt_samplename}_${dorado_model}_sorted.bam $ref_fasta $output_dir/cutesv_output/${wt_samplename}_${dorado_model}_cutesv.vcf $output_dir/cutesv_output \
    --max_cluster_bias_INS 100 \
    --diff_ratio_merging_INS 0.3 \
    --max_cluster_bias_DEL 100 \
    --diff_ratio_merging_DEL 0.3 &

# Run the second cuteSV command in the background
cuteSV $output_dir/sorted_bams/${case_samplename}_${dorado_model}_sorted.bam $ref_fasta $output_dir/cutesv_output/${case_samplename}_${dorado_model}_cutesv.vcf $output_dir/cutesv_output \
    --max_cluster_bias_INS 100 \
    --diff_ratio_merging_INS 0.3 \
    --max_cluster_bias_DEL 100 \
    --diff_ratio_merging_DEL 0.3 &

# Wait for both processes to finish
wait
