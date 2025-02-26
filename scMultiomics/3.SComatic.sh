'''
### Installation

source /home/project/11003581/Tools/miniforge3/bin/activate
conda create --prefix /home/project/11003581/conda-envs/SComatic python=3.7 samtools datamash bedtools
conda activate /home/project/11003581/conda-envs/SComatic

#git clone the repo to a location

pip install -r /home/project/11003581/Tools/SComatic/requirements.txt

about-time==3.1.1
numpy==1.21.6
numpy-groupies==0.9.14
pandas==1.3.5
pybedtools==0.8.1
pysam==0.16.0.1
rpy2==2.9.4
scipy==1.7.3

Rscript /home/project/11003581/Tools/SComatic/r_requirements_install.R

### gunzip the ref files

gunzip PoNs/PoN.scRNAseq.hg38.tsv.gz
gunzip PoNs/PoN.scATACseq.hg38.tsv.gz 
gunzip RNAediting/AllEditingSites.hg38.txt.gz


## scripts
https://github.com/cortes-ciriano-lab/SComatic/tree/main/scripts

'''
### Requirement 
# 1 annotaion file as below in tsv format (do not include doublets)
Index Cell_type
AAACCTGCATGCTAGT  Epithelial
AAACCTGGTAGCCTAT  Epithelial
AAACCTGGTTGTCGCG  Epithelial
AAACCTGTCATGTGGT  Epithelial
AAACCTGTCCTTGGTC  Epithelial
AAACCTGTCGGATGTT  T_cell
AAACCTGTCGTACGGC  T_cell
AAACCTGTCTTGCAAG  T_cell
AAACGGGAGACGCACA  T_cell


########### PBS SCRIPT ###########
#!/bin/bash
  
#PBS -l select=2:ncpus=64:mem=128g
#PBS -l walltime=6:00:00
#PBS -P 11003581
#PBS -j oe
#PBS -N run-scomatic


# Change to the directory where the job was submitted
cd $PBS_O_WORKDIR

source /home/project/11003581/Tools/miniforge3/bin/activate
conda activate /home/project/11003581/conda-envs/SComatic

module load r

#Step 1: Splitting alignment file into cell-type-specific bams

SCOMATIC=/home/project/11003581/Tools/SComatic/
sample=10kpbmc
output_dir=/home/users/nus/ash.ps/scratch/mulitomics/SComatic/
output_dir1=$output_dir/Step1_BamCellTypes
mkdir -p $output_dir1

python $SCOMATIC/scripts/SplitBam/SplitBamCellTypes.py \
        --bam /home/users/nus/ash.ps/scratch/mulitomics/analysis/10k_pbmc/outs/gex_possorted_bam.bam \
        --meta /home/users/nus/ash.ps/scratch/mulitomics/gex-class/gex-classification.tsv \
        --id ${sample} \
        --n_trim 5 \
        --max_nM 5 \
        --max_NH 1 \
        --outdir $output_dir1

#### Pre step 2 - filtering out chr contigs


module load samtools/1.15.1

filtered_out=$output_dir1/filtered
mkdir -p "$filtered_out"

# List of chromosomes to keep
chromosomes="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"

# Process each BAM file in the input directory
for input_bam in "$input_dir"/*.bam; do
    # Get the base name of the input file
    base_name=$(basename "$input_bam")
    
    # Define the output BAM file path
    output_bam="${output_dir}/${base_name%.bam}_filtered.bam"
    
    # Filter BAM file
    samtools view -h "$input_bam" $chromosomes | samtools sort -o "$output_bam"
    
    # Index the filtered BAM file
    samtools index "$output_bam"
    
    echo "Processed and indexed: $base_name"
done

echo "All BAM files have been filtered and indexed."




#Step 2: Collecting base count information

REF=/home/project/11003581/Ref/Homo_sapiens/GATK/hg38/Homo_sapiens_assembly38.fasta

output_dir1=$output_dir/Step1_BamCellTypes/filtered/
output_dir2=$output_dir/Step2_BaseCellCounts
mkdir -p $output_dir2

for bam in $(ls -d $output_dir1/*bam);do
  
  # Cell type
  cell_type=$(basename $bam | awk -F'.' '{print $(NF-1)}')

  # Temp folder
  temp=$output_dir2/temp_${cell_type}
  mkdir -p $temp

  # Command line to submit to cluster
  python $SCOMATIC/scripts/BaseCellCounter/BaseCellCounter.py --bam $bam \
    --ref $REF \
    --chrom all \
    --out_folder $output_dir2 \
    --min_bq 30 \
    --tmp_dir $temp \
    --nprocs 16

  rm -rf $temp
done

#Step 3: Merging base count matrices

sample=Example
output_dir3=$output_dir/Step3_BaseCellCountsMerged
mkdir -p $output_dir3

python $SCOMATIC/scripts/MergeCounts/MergeBaseCellCounts.py --tsv_folder ${output_dir2} \
  --outfile ${output_dir3}/${sample}.BaseCellCounts.AllCellTypes.tsv


#Step 4: Detection of somatic mutations

# Step 4.1
output_dir4=$output_dir/Step4_VariantCalling
mkdir -p $output_dir4

sample=Example
REF=$SCOMATIC/example_data/chr10.fa

python $SCOMATIC/scripts/BaseCellCalling/BaseCellCalling.step1.py \
          --infile ${output_dir3}/${sample}.BaseCellCounts.AllCellTypes.tsv \
          --outfile ${output_dir4}/${sample} \
          --ref $REF

# Step 4.2
editing=$SCOMATIC/RNAediting/AllEditingSites.hg38.txt
PON=$SCOMATIC/PoNs/PoN.scRNAseq.hg38.tsv

python $SCOMATIC/scripts/BaseCellCalling/BaseCellCalling.step2.py \
          --infile ${output_dir4}/${sample}.calling.step1.tsv \
          --outfile ${output_dir4}/${sample} \
          --editing $editing \
          --pon $PON


#Optional step - TNM
sample=Example
output_dir8=$output_dir/TrinucleotideContext
output_dir4=$output_dir/Step4_VariantCalling # Already defined in previous steps
mkdir -p $output_dir8

echo ${output_dir4}/${sample}.calling.step1.tsv > ${output_dir8}/step1_files.txt

python $SCOMATIC/scripts/TrinucleotideBackground/TrinucleotideContextBackground.py \
        --in_tsv ${output_dir8}/step1_files.txt \
        --out_file ${output_dir8}/TrinucleotideBackground.txt