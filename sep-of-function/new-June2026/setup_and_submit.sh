#!/bin/bash
# Run this ON THE CLUSTER from the directory containing master_samplesheet.csv
# (the full fastq sheet with all 10 WT_3_XX samples)

BASE=/home/users/nus/ash.ps/scratch/sep-fun/nf-sarek/variant_calling/WT-old
MASTER=$BASE/master_samplesheet.csv   # place your full fastq sheet here first
NORMAL_R1=/home/project/11003581/AV-lab/MCF10A_grandparent/fastq/SRR25739470_1.fastq.gz
NORMAL_R2=/home/project/11003581/AV-lab/MCF10A_grandparent/fastq/SRR25739470_2.fastq.gz

SAMPLES=(WT_P30_1 WT_P30_2 WT_P30_3 WT_P30_4 WT_P30_5 WT_P30_6 WT_P30_7 WT_P30_8 WT_P30_9 WT_P30_10)

for SAMPLE in "${SAMPLES[@]}"; do
  SDIR=$BASE/$SAMPLE
  mkdir -p $SDIR

  # header
  echo "patient,sex,status,sample,lane,fastq_1,fastq_2" > $SDIR/samplesheet.csv
  # tumor lines for this sample, patient forced to $SAMPLE (already matches)
  awk -F',' -v s="$SAMPLE" 'NR>1 && $4==s {print}' $MASTER >> $SDIR/samplesheet.csv
  # normal row, patient matches this sample -> no cross-sample reuse, SM tag will be consistent
  echo "${SAMPLE},XX,0,gMCF10A,1,${NORMAL_R1},${NORMAL_R2}" >> $SDIR/samplesheet.csv

  cat > $SDIR/run.pbs << EOF
#!/bin/bash
#PBS -l select=1:ncpus=64:mem=128g
#PBS -l walltime=24:00:00
#PBS -P 11003581
#PBS -N wtold-${SAMPLE}
#PBS -M ash.ps@nus.edu.sg
#PBS -j oe

cd \$PBS_O_WORKDIR
module purge
module load java/17.0.6-jdk
module load nextflow/25.10.4
module load singularity

export NXF_HOME=/home/project/11003581/Tools/nf-home
export NXF_SINGULARITY_CACHEDIR=/home/users/nus/ash.ps/scratch/sing/cache/
export SINGULARITY_CACHEDIR=/home/users/nus/ash.ps/scratch/sing/cache/

WORKDIR=${SDIR}

nextflow run nf-core/sarek -r 3.9.0 \\
   -profile singularity \\
   --input \$WORKDIR/samplesheet.csv \\
   --outdir \$WORKDIR \\
   --tools mutect2,strelka,ascat,manta \\
   -params-file /home/project/11003581/Tools/nf-sarek-params.yaml
EOF

  qsub $SDIR/run.pbs
done