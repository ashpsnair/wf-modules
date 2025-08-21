#PBS -l select=1:ncpus=128:mem=256g
#PBS -l walltime=04:00:00
#PBS -P 11003581
#PBS -N slice-bam-rad51
#PBS -j oe

cd $PBS_O_WORKDIR
module load samtools

REFERENCE="/home/project/11003581/Ref/sarek-refs/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"
BASE_INPUT_DIR="/home/users/nus/ash.ps/scratch/sep-fun-analysis/RAD51_S181P/preprocessing/recalibrated"
OUTPUT_DIR="/home/users/nus/ash.ps/scratch/sep-fun-analysis/RAD51_S181P/preprocessing/chr15-bams"
mkdir -p "$OUTPUT_DIR"

find "$BASE_INPUT_DIR" -type f -name "*.cram" | while read -r cramfile; do
  filename=$(basename "$cramfile" .cram)
  outbam="$OUTPUT_DIR/${filename}_chr15.bam"

  # Detect whether contigs are '15' or 'chr15'
  if samtools idxstats "$cramfile" | cut -f1 | grep -qx "chr15"; then
    CHR="chr15"
  elif samtools idxstats "$cramfile" | cut -f1 | grep -qx "15"; then
    CHR="15"
  else
    echo "Could not find chromosome 15 name in $cramfile"; continue
  fi

  echo "Processing $cramfile -> $outbam ($CHR)"
  # Extract region, keep header, sort just to be safe, then index
  samtools view -@ 16 -T "$REFERENCE" -b -h "$cramfile" "$CHR" \
    | samtools sort -@ 16 -o "$outbam" -
  samtools index -@ 16 "$outbam"
done
