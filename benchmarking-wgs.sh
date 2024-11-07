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
minimap2_path="/home/project/11003581/Tools/minimap2-2.28_x64-linux/"

#sample
sample_csv="/home/users/nus/ash.ps/scratch/directRNA-workflow/sample.csv"
output_dir="/home/users/nus/ash.ps/scratch/directRNA-workflow"
ref_fasta="/home/project/11003581/Ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

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

$dorado_bin_path/dorado basecaller hac,4mC_5mC,5mCG_5hmCG,5mC_5hmC,6mA $pod5_path \
--trim all \
--reference $ref_fasta \
> $output_dir/dorado_output/${samplename}_hac.bam

$dorado_bin_path/dorado summary $output_dir/dorado_output/${samplename}_hac.bam > $output_dir/dorado_output/${samplename}_dorado_summary.txt

##### Sorting & Indexing BAM files ######
mkdir $output_dir/sorted_bams

samtools sort -@ 4 $output_dir/dorado_output/${samplename}_hac.bam \
-o $output_dir/sorted_bams/${samplename}_sorted.bam

########## SV calling using Sniffles #####
mkdir $output_dir/sniffles_out

module load python/3.12.1-gcc11

sniffles -i $output_dir/sorted_bams/${samplename}_sorted.bam \
        -v $output_dir/sniffles_out/${samplename}_sniffles.vcf \
        --reference $ref_fasta


######### SNV calling using Deepvariant ############
mkdir $output_dir/deepvariant_out
module load singularity

# Pull the image.
BIN_VERSION="1.6.1"
singularity pull docker://google/deepvariant:"${BIN_VERSION}-gpu"

# Run DeepVariant.
singularity run --nv -B /usr/lib/locale/:/usr/lib/locale/ \
  docker://google/deepvariant:"${BIN_VERSION}-gpu" \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=WGS \ 
  --ref=$ref_fasta \
  --reads=$output_dir/sorted_bams/${samplename}_sorted.bam \
  --output_vcf=$output_dir/deepvariant_out/${samplename}_deepvar.vcf.gz \
  --output_gvcf=$output_dir/deepvariant_out/${samplename}_deepvar.g.vcf.gz \
  --num_shards=32 