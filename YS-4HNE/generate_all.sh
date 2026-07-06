'''
This creates analysis/{4HNE_day0,4HNE_day35,MGO_day0,MGO_day17,MGO_day35}/setup_and_submit.sh, each pre-filled with its own TREAT/DAY. Then, per folder:
cd analysis/4HNE_day0 && bash setup_and_submit.sh
This makes one subfolder per Sample, writes samplesheet.csv (normal + all lanes, flowcell-disambiguated), writes run_sarek.pbs, and submits via qsub — for every sample in that treatment/day group.

'''

#!/bin/bash
set -euo pipefail

ANALYSIS_ROOT=/home/users/nus/ash.ps/scratch/YS-4HNE/analysis
mkdir -p "$ANALYSIS_ROOT"

declare -A COMBOS=(
  [4HNE_day0]="4HNE 0"
  [4HNE_day35]="4HNE 35"
  [MGO_day0]="MGO 0"
  [MGO_day17]="MGO 17"
  [MGO_day35]="MGO 35"
)

for FOLDER in "${!COMBOS[@]}"; do
  read -r TREAT DAY <<< "${COMBOS[$FOLDER]}"
  DIR="$ANALYSIS_ROOT/$FOLDER"
  mkdir -p "$DIR"

  cat > "$DIR/setup_and_submit.sh" <<SCRIPT
#!/bin/bash
set -euo pipefail

TREAT="$TREAT"
DAY="$DAY"
BASE=/home/users/nus/ash.ps/scratch/YS-4HNE
METADATA=\$BASE/metadata.tsv
FASTQ_LIST=\$BASE/fastq-list.txt
NORMAL_FQ1=/home/project/11003581/AV-lab/MCF10A_grandparent/fastq/SRR25739470_1.fastq.gz
NORMAL_FQ2=/home/project/11003581/AV-lab/MCF10A_grandparent/fastq/SRR25739470_2.fastq.gz
PARAMS=/home/project/11003581/Tools/nf-sarek-params.yaml
GROUPDIR=\$(cd "\$(dirname "\$0")" && pwd)

# short sample names (strip _FDSW... suffix), e.g. A_2_10_FDSW202588819 -> A_2_10
SAMPLES=\$(awk -F'\\t' -v t="\$TREAT" -v d="\$DAY" 'NR>1 && \$3==t && \$4==d {s=\$6; sub(/_FDSW.*/,"",s); print s}' "\$METADATA" | sort -u)

for SAMPLE in \$SAMPLES; do
  # map short name back to the full raw Sample value used in metadata column 6
  RAWSAMPLE=\$(awk -F'\\t' -v t="\$TREAT" -v d="\$DAY" -v short="\$SAMPLE" \
    'NR>1 && \$3==t && \$4==d {s=\$6; sub(/_FDSW.*/,"",s); if (s==short) {print \$6; exit}}' "\$METADATA")

  SDIR="\$GROUPDIR/\$SAMPLE"
  mkdir -p "\$SDIR"
  SHEET="\$SDIR/samplesheet.csv"
  echo "patient,sex,status,sample,lane,fastq_1,fastq_2" > "\$SHEET"
  echo "\$SAMPLE,XX,0,gMCF10A_\${SAMPLE},L1,\$NORMAL_FQ1,\$NORMAL_FQ2" >> "\$SHEET"

  while read -r STUB; do
    [[ -z "\$STUB" ]] && continue
    FQ1=\$(grep -F -- "\${STUB}_1.fq.gz" "\$FASTQ_LIST" || true)
    FQ2=\$(grep -F -- "\${STUB}_2.fq.gz" "\$FASTQ_LIST" || true)
    if [[ -z "\$FQ1" || -z "\$FQ2" ]]; then
      echo "WARNING: missing fastq for \$STUB" >&2
      continue
    fi
    if [[ \$STUB =~ _([A-Z0-9]+)_(L[0-9]+)\$ ]]; then
      LANE_ID="\${BASH_REMATCH[2]}_\${BASH_REMATCH[1]}"
    else
      LANE_ID="\$STUB"
    fi
    echo "\$SAMPLE,XX,1,\$SAMPLE,\$LANE_ID,\$FQ1,\$FQ2" >> "\$SHEET"
  done < <(awk -F'\\t' -v t="\$TREAT" -v d="\$DAY" -v s="\$RAWSAMPLE" '\$3==t && \$4==d && \$6==s {print \$1}' "\$METADATA")

  PBS="\$SDIR/run_sarek.pbs"
  cat > "\$PBS" <<EOF
#!/bin/bash
#PBS -l select=1:ncpus=64:mem=128g
#PBS -l walltime=24:00:00
#PBS -P 11003581
#PBS -N ${FOLDER}-\${SAMPLE}
#PBS -M ash.ps@nus.edu.sg
#PBS -j oe

cd \\\$PBS_O_WORKDIR
module purge
module load java/17.0.6-jdk
module load nextflow/25.10.4
module load singularity

export NXF_HOME=/home/project/11003581/Tools/nf-home
export NXF_SINGULARITY_CACHEDIR=/home/users/nus/ash.ps/scratch/sing/cache/
export SINGULARITY_CACHEDIR=/home/users/nus/ash.ps/scratch/sing/cache/

WORKDIR=\$SDIR

nextflow run nf-core/sarek -r 3.9.0 \\\\
   -profile singularity \\\\
   --input \\\$WORKDIR/samplesheet.csv \\\\
   --outdir \\\$WORKDIR \\\\
   --tools mutect2,strelka,ascat,manta \\\\
   -params-file \$PARAMS
EOF

  echo "Submitting ${FOLDER}/\$SAMPLE"
  cd "\$SDIR" && qsub "\$PBS"
done
SCRIPT

  chmod +x "$DIR/setup_and_submit.sh"
  echo "Created $DIR/setup_and_submit.sh"
done