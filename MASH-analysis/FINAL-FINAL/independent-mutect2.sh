#!/bin/bash
#PBS -l select=2:ncpus=64:mem=128gb
#PBS -l walltime=12:00:00
#PBS -N mutect2_joint
#PBS -j oe

cd $PBS_O_WORKDIR

# Load modules
module load samtools/1.15.1
module load python/3.10.9
module load parallel

export SHEET="/home/users/nus/enambis/scratch/NCCS-MASH/gatk-mutect2/samplesheet.csv"
export OUTDIR="/home/users/nus/enambis/scratch/NCCS-MASH/gatk-mutect2/crams"
mkdir -p "$OUTDIR"

tail -n +2 "$SHEET" | cut -d',' -f5 | parallel -j 16 --jobs 16 '
  cram={}
  base=$(basename "$cram" .recal.cram)
  sorted_cram="$OUTDIR/${base}.recal.sorted.cram"
  echo "➡ Sorting $cram to $sorted_cram"
  samtools sort -T "$OUTDIR/tmp_${base}" -@ 8 -O cram -o "$sorted_cram" "$cram"
  samtools index "$sorted_cram"
'



###-------------------------------------

###-------------------------------------

#!/bin/bash
#PBS -l select=2:ncpus=64:mem=128gb
#PBS -l walltime=12:00:00
#PBS -N mutect2_joint
#PBS -j oe

cd $PBS_O_WORKDIR

# Load modules
module load samtools/1.15.1
module load python/3.10.9

# Tool and reference paths
GATK="/home/project/11003581/Tools/gatk-4.6.1.0/gatk"
SHEET="/home/users/nus/enambis/scratch/NCCS-MASH/gatk-mutect2/samplesheet.csv"
REFERENCE="/home/project/11003581/Ref/sarek-refs/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"
GERMLINE_RESOURCE="/home/project/11003581/Ref/sarek-refs/Homo_sapiens/GATK/GRCh38/Annotation/intervals/af-only-gnomad.hg38.vcf.gz"
PON="/home/project/11003581/Ref/pons/1000g_pon.hg38.vcf.gz"
INTERVALS="/home/project/11003581/Ref/sarek-refs/Homo_sapiens/GATK/GRCh38/Annotation/intervals/wgs_calling_regions.hg38.bed"

# CRAM directory (all sorted & indexed)
CRAM_DIR="/home/users/nus/enambis/scratch/NCCS-MASH/gatk-mutect2/crams"

# Output root directory
JOINT_OUTDIR="/home/users/nus/enambis/scratch/NCCS-MASH/gatk-mutect2/joint-calls"
mkdir -p "$JOINT_OUTDIR"
GLOBAL_LOG="$JOINT_OUTDIR/mutect2_pipeline.log"
echo "===== Starting Mutect2 Pipeline: $(date) =====" > "$GLOBAL_LOG"

# Process each patient
patients=$(tail -n +2 "$SHEET" | cut -d',' -f1 | sort | uniq)

for patient in $patients; do
  echo "🔍 Processing patient: $patient" | tee -a "$GLOBAL_LOG"

  normal_sample=$(awk -F',' -v p="$patient" '$1 == p && $3 == 0 {print $4; exit}' "$SHEET")
  normal_cram="${CRAM_DIR}/${normal_sample}.recal.sorted.cram"

  OUTDIR="${JOINT_OUTDIR}/${patient}"
  mkdir -p "$OUTDIR"

  # Loop through tumor samples
  awk -F',' -v p="$patient" '$1 == p && $3 == 1' "$SHEET" | while IFS=',' read -r _ _ _ _ tumor_sample _ _; do
    tumor_cram="${CRAM_DIR}/${tumor_sample}.recal.sorted.cram"
    prefix="${patient}_${tumor_sample}_vs_${normal_sample}"

    echo "▶ $prefix" | tee -a "$GLOBAL_LOG"

    echo "[1/5] Mutect2 ..." | tee -a "$GLOBAL_LOG"
    $GATK Mutect2 \
      -R "$REFERENCE" \
      -I "$tumor_cram" -tumor "$tumor_sample" \
      -I "$normal_cram" -normal "$normal_sample" \
      --germline-resource "$GERMLINE_RESOURCE" \
      --panel-of-normals "$PON" \
      -L "$INTERVALS" \
      --f1r2-tar-gz "$OUTDIR/${prefix}.f1r2.tar.gz" \
      -O "$OUTDIR/${prefix}.unfiltered.vcf.gz" 2>> "$GLOBAL_LOG"

    echo "[2/5] LearnReadOrientationModel ..." | tee -a "$GLOBAL_LOG"
    $GATK LearnReadOrientationModel \
      -I "$OUTDIR/${prefix}.f1r2.tar.gz" \
      -O "$OUTDIR/${prefix}.read-orientation-model.tar.gz" 2>> "$GLOBAL_LOG"

    echo "[3/5] GetPileupSummaries ..." | tee -a "$GLOBAL_LOG"
    $GATK GetPileupSummaries \
      -I "$tumor_cram" \
      -V "$GERMLINE_RESOURCE" \
      -L "$INTERVALS" \
      -R "$REFERENCE" \
      -O "$OUTDIR/${prefix}.pileups.table" 2>> "$GLOBAL_LOG"

    echo "[4/5] CalculateContamination ..." | tee -a "$GLOBAL_LOG"
    $GATK CalculateContamination \
      -I "$OUTDIR/${prefix}.pileups.table" \
      -O "$OUTDIR/${prefix}.contamination.table" \
      --tumor-segmentation "$OUTDIR/${prefix}.segments.table" 2>> "$GLOBAL_LOG"

    echo "[5/5] FilterMutectCalls ..." | tee -a "$GLOBAL_LOG"
    $GATK FilterMutectCalls \
      -V "$OUTDIR/${prefix}.unfiltered.vcf.gz" \
      -R "$REFERENCE" \
      --contamination-table "$OUTDIR/${prefix}.contamination.table" \
      --tumor-segmentation "$OUTDIR/${prefix}.segments.table" \
      --ob-priors "$OUTDIR/${prefix}.read-orientation-model.tar.gz" \
      -O "$OUTDIR/${prefix}.filtered.vcf.gz" 2>> "$GLOBAL_LOG"

    echo "✅ Done: $prefix at $(date)" | tee -a "$GLOBAL_LOG"
    echo "----------------------------------------------------------" | tee -a "$GLOBAL_LOG"
  done
done

echo "🎉 All patients done. Combined log: $GLOBAL_LOG"

