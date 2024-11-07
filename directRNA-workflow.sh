#!/bin/bash

#PBS -l select=1:ngpus=1
#PBS -l walltime=120:00:00
#PBS -P personal-ash.ps
#PBS -N direct-rna-wf-check
#PBS -j oe
#PBS -q ai

module load samtools/1.15.1

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

##Tools Directory
dorado_bin_path="/home/project/11003581/Tools/dorado-0.7.3-linux-x64/bin"
nanocount_bin_path="/home/users/nus/ash.ps/anaconda3/bin"
minimap2_path="/home/project/11003581/Tools/minimap2-2.28_x64-linux/"

#sample
sample_csv="/home/users/nus/ash.ps/scratch/directRNA-workflow/sample.csv"
output_dir="/home/users/nus/ash.ps/scratch/directRNA-workflow"
ref_fasta="/home/project/11003581/Ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
transcript_fasta="/home/project/11003581/Ref/gencode.v46.transcripts.fa"

#Running dorado for first sample
# Read the first data line from the CSV file (skip the header)
samplename=$(awk -F ',' 'NR==2 {print $1}' "$sample_csv")  # Assuming samplename is in the first column
pod5_path=$(awk -F ',' 'NR==2 {print $2}' "$sample_csv")  # Assuming fastq_path is in the 2nd column
fastq_path=$(awk -F ',' 'NR==2 {print $3}' "$sample_csv")  # Assuming fastq_path is in the 3rd column

#Running dorado command
echo "##########################################################"
echo "Running the dorado command"
echo "##########################################################"

# Print the dorado command
mkdir $output_dir/dorado_output

$dorado_bin_path/dorado basecaller sup,m6A,pseU $pod5_path \
--estimate-poly-a -v \
--reference $ref_fasta \
--mm2-preset splice -k 14 \
> $output_dir/dorado_output/${samplename}_sup.bam

$dorado_bin_path/dorado summary $output_dir/dorado_output/${samplename}_hac.bam > $output_dir/r_aligned_bams/${samplename}_dorado_summary.txt

##### Sorting & Indexing BAM files ######
mkdir $output_dir/sorted_bams

samtools sort -@ 4 $output_dir/dorado_output/${samplename}_hac.bam \
-o $output_dir/sorted_bams/${samplename}_sorted.bam

########## Trasnscript alignment #####
mkdir $output_dir/t_aligned_bams

#Running minimap for transcript align
echo "##########################################################"
echo "Running minimap for transcript align"
echo "##########################################################"

#concatenate fastq files
zcat "$fastq_path"/*.fastq.gz > "$output_dir/${samplename}_merged.fastq"
$minimap2_path/minimap2 -ax splice \
-t 4 -ax map-ont -N 10 $transcript_fasta $output_dir/${samplename}_merged.fastq | samtools view -bh > $output_dir/t_aligned_bams/${samplename}_t_align.bam

########## Running Nanocount #######
mkdir $output_dir/nanocount_out

#Running Nanocount
echo "##########################################################"
echo "Nanocount"
echo "##########################################################"
#nanocount command
$nanocount_bin_path/NanoCount \
-i $output_dir/t_aligned_bams/${samplename}_t_align.bam \
-o $output_dir/nanocount_out/${samplename}_counts.txt
