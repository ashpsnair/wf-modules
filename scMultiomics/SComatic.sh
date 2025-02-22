'''
### Installation

source /home/project/11003581/Tools/miniforge3/bin/activate
conda create --prefix /home/project/11003581/conda-envs/SComatic python=3.7 r-base=3.6.1 samtools datamash bedtools
conda activate /home/project/11003581/conda-envs/SComatic

pip install -r requirements.txt
[requirements.txt = https://github.com/cortes-ciriano-lab/SComatic/blob/main/requirements.txt ]

Rscript r_requirements_install.R
[ r_requirements_install.R = https://github.com/cortes-ciriano-lab/SComatic/blob/main/r_requirements_install.R ]

### Get the reference files 
N
https://github.com/cortes-ciriano-lab/SComatic/tree/main/RNAediting

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

#Step 1: Splitting alignment file into cell-type-specific bams

SCOMATIC=/home/project/11003581/Tools/SComatic/
sample=10kpbmc
output_dir1=$output_dir/Step1_BamCellTypes
mkdir -p $output_dir1

python $SCOMATIC/scripts/SplitBam/SplitBamCellTypes.py --bam $SCOMATIC/example_data/Example.scrnaseq.bam \
        --meta $SCOMATIC/example_data/Example.cell_barcode_annotations.tsv \
        --id ${sample} \
        --n_trim 5 \
        --max_nM 5 \
        --max_NH 1 \
        --outdir $output_dir1

#Step 2: Collecting base count information

REF=$SCOMATIC/example_data/chr10.fa

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
    --nprocs 1

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