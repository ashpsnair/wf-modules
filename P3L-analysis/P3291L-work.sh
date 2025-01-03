#!/bin/bash

#PBS -l select=1:ngpus=1
#PBS -l walltime=10:00:00
#PBS -P 11003581
#PBS -N process-p3l-lab
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load python

##Tools Directory
dorado_bin_path="/home/project/11003581/Tools/dorado-0.7.3-linux-x64/bin"
minimap2_path="/home/project/11003581/Tools/minimap2-2.28_x64-linux/"

#sample details
output_dir="/home/users/nus/ash.ps/scratch/P3L-lab-analysis/"
ref_fasta="/home/project/11003581/Ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

'''
#Codes for converting fast5 to pod5
wt_pod5_dir="/home/users/nus/ash.ps/scratch/MCF10A-lab/pod5/WT/"
p3l_pod5_dir="/home/users/nus/ash.ps/scratch/MCF10A-lab/pod5/P3L/"
wt_fast5_dir="/home/users/nus/ash.ps/scratch/MCF10A-lab/WT/MCF10A_WT_fast5_pass"
p3l_fast5_dir="/home/users/nus/ash.ps/scratch/MCF10A-lab/P3L/MCF10A_P3L_2_fast5_pass"

#creating output dirs
mkdir $wt_pod5_dir
mkdir $p3l_pod5_dir

#Convert each fast5 to its relative converted output. 
pod5 convert fast5 $wt_fast5_dir/*.fast5 --output $wt_pod5_dir/ --one-to-one $wt_fast5_dir/
pod5 convert fast5 $p3l_fast5_dir/*.fast5 --output $p3l_pod5_dir/ --one-to-one $p3l_fast5_dir/
'''


wt_pod5_dir="/home/users/nus/ash.ps/scratch/MCF10A-lab/pod5/WT/"
p3l_pod5_dir="/home/users/nus/ash.ps/scratch/MCF10A-lab/pod5/P3L/"

####### Running dorado
mkdir $output_dir/dorado_output

$dorado_bin_path/dorado basecaller hac $wt_pod5_dir \
--trim all \
--reference $ref_fasta \
> $output_dir/dorado_output/WT_MCF10A_hac.bam

$dorado_bin_path/dorado basecaller hac $p3l_pod5_dir \
--trim all \
--reference $ref_fasta \
> $output_dir/dorado_output/p3292l_hac.bam


$dorado_bin_path/dorado summary $output_dir/dorado_output/WT_p53ko_hac.bam > $output_dir/dorado_output/WT_MCF10A_dorado_summary.txt
$dorado_bin_path/dorado summary $output_dir/dorado_output/p3292l_hac.bam > $output_dir/dorado_output/p3292l_dorado_summary.txt


##### Sorting & Indexing BAM files ######
mkdir $output_dir/sorted_bams

samtools sort $output_dir/dorado_output/WT_MCF10A_hac.bam -o $output_dir/sorted_bams/WT_MCF10A_hac_sorted.bam
samtools index $output_dir/sorted_bams/WT_MCF10A_hac_sorted.bam

samtools sort $output_dir/dorado_output/p3292l_hac.bam -o $output_dir/sorted_bams/p3292l_hac_sorted.bam
samtools index $output_dir/sorted_bams/p3292l_hac_sorted.bam

########## SV calling using Sniffles #####
mkdir $output_dir/sniffles_out

module load python/3.12.1-gcc11

sniffles -i $output_dir/sorted_bams/WT_MCF10A_hac_sorted.bam \
        -v $output_dir/sniffles_out/WT_MCF10A_sniffles.vcf \
        --reference $ref_fasta

sniffles -i $output_dir/sorted_bams/p3292l_hac_sorted.bam \
        -v $output_dir/sniffles_out/p3292l_sniffles.vcf \
        --reference $ref_fasta

'''
#SV calling using mosaic argument
sniffles -i $output_dir/sorted_bams/p3292l_hac_sorted.bam \
        -v $output_dir/sniffles_out/p3292l_mosaic_sniffles.vcf \
        --mosaic \
        --reference $ref_fasta

'''
####### Sniffles Plot
/home/users/nus/ash.ps/anaconda3/bin/python3.8 -m sniffles2_plot -i p3292l_sniffles.vcf -o ./output

####### Running CuteSV

mkdir -p $output_dir/cutesv_output  # Use -p to avoid error if the directory already exists

# Run the first cuteSV command in the background
cuteSV $output_dir/sorted_bams/p3292l_hac_sorted.bam $ref_fasta $output_dir/cutesv_output/p3292l_cutesv.vcf $output_dir/cutesv_output \
    --max_cluster_bias_INS 100 \
    --diff_ratio_merging_INS 0.3 \
    --max_cluster_bias_DEL 100 \
    --diff_ratio_merging_DEL 0.3 &

# Run the second cuteSV command in the background
cuteSV $output_dir/sorted_bams/WT_MCF10A_hac_sorted.bam $ref_fasta $output_dir/cutesv_output/WT_MCF10A_cutesv.vcf $output_dir/cutesv_output \
    --max_cluster_bias_INS 100 \
    --diff_ratio_merging_INS 0.3 \
    --max_cluster_bias_DEL 100 \
    --diff_ratio_merging_DEL 0.3 &

# Wait for both processes to finish
wait

echo "Both cuteSV commands have completed."