#!/bin/bash
#PBS -P 11003581
#PBS -l select=1:ncpus=64:mem=64g
#PBS -l walltime=04:00:00
#PBS -N fc_per_bam_merge
#PBS -j oe

set -e
cd $PBS_O_WORKDIR

LOGFILE="featurecounts_per_bam_merge.log"
exec > >(tee -a "$LOGFILE") 2>&1

echo "[$(date)] Job started on $(hostname)"
echo "[$(date)] Working directory: $PBS_O_WORKDIR"

# Load Conda env
module load miniforge3
source activate /home/project/11003581/conda-envs/riboseq_qc

# Paths
STAR_SORTED="/home/users/nus/ash.ps/scratch/JQQ/analysis-july/mamba/alignment/star/sorted"
GTF="/home/users/nus/ash.ps/scratch/refs/Homo_sapiens.GRCh38.114.gtf"
OUTDIR="/home/users/nus/ash.ps/scratch/JQQ/analysis-july/downstream/fc_per_sample"
MERGED="$OUTDIR/cds_counts_all.txt"

mkdir -p "$OUTDIR"

echo "[$(date)] Step 1: Running featureCounts for each BAM..."
for bam in "$STAR_SORTED"/*.genome.sorted.bam; do
    sample=$(basename "$bam" .genome.sorted.bam)
    echo "Processing: $sample"
    
    featureCounts \
        -T 64 \
        -t CDS \
        -g gene_id \
        -a "$GTF" \
        -o "$OUTDIR/${sample}_cds_counts.txt" \
        "$bam"

    # Simplify output: keep Geneid and count column only
    awk 'NR>2 {print $1 "\t" $NF}' "$OUTDIR/${sample}_cds_counts.txt" \
        > "$OUTDIR/${sample}_counts_simple.txt"
done

echo "[$(date)] Step 2: Merging simplified counts into one matrix..."
Rscript - <<EOF
files <- list.files(path="$OUTDIR", pattern="_counts_simple.txt$", full.names=TRUE)
data_list <- lapply(files, function(f) {
  sample_name <- sub("_counts_simple.txt$", "", basename(f))
  df <- read.table(f, header=FALSE, col.names=c("Geneid", sample_name))
  return(df)
})
merged <- Reduce(function(x, y) merge(x, y, by="Geneid", all=TRUE), data_list)
write.table(merged, "$MERGED", sep="\t", row.names=FALSE, quote=FALSE)
EOF

echo "[$(date)] Step 3: Merged counts saved to $MERGED"
echo "[$(date)] Job finished."
