'''
################# to be run inside the directory to start the sigprofile run. #################

qsub \
    -o /home/users/nus/ash.ps/scratch/YS-4HNE/downstream/mutsig/logs \
    -j oe \
    -v PROJECT_NAME=4HNE,BASE_DIR=/home/users/nus/ash.ps/scratch/YS-4HNE/downstream/mutsig/ \
    /home/project/11003581/Tools/scripts/sigprofiler.pbs

'''


#!/bin/bash
#PBS -l select=1:ncpus=128:mem=256g
#PBS -l walltime=06:00:00
#PBS -P 11003581
#PBS -N SigProf_Matrix
#PBS -j oe
#PBS -M ash.ps@nus.edu.sg

set -euo pipefail

module purge
module load gcc/15.2.0-nscc
module load python/3.12.1-gcc11

source /home/project/11003581/conda-envs/SigProfiler/bin/activate

# These are passed in at submission time
: "${PROJECT_NAME:?PROJECT_NAME is required}"
: "${BASE_DIR:?BASE_DIR is required}"

VCF_DIR="${BASE_DIR}/vcfs"
OUTPUT_DIR="${BASE_DIR}"

MATRIX_SCRIPT="/home/project/11003581/Tools/scripts/matrix_generator.py"
EXTRACT_SCRIPT="/home/project/11003581/Tools/scripts/run_extractor.py"
PBS_SCRIPT="/home/project/11003581/Tools/scripts/sigprofiler.pbs"

LOG_DIR="${OUTPUT_DIR}/logs"
mkdir -p "${LOG_DIR}"

if [[ -z "${MUTATION_TYPE:-}" ]]; then
    LOG_NAME="MATRIX"
else
    LOG_NAME="${MUTATION_TYPE}"
fi

LOG_FILE="${LOG_DIR}/${PROJECT_NAME}_${LOG_NAME}_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "=================================================="
echo "Started : $(date)"
echo "Job ID  : ${PBS_JOBID}"
echo "Node    : $(hostname)"
echo "Project : ${PROJECT_NAME}"

if [[ -z "${MUTATION_TYPE:-}" ]]; then
    echo "Mode    : MATRIX"
else
    echo "Mode    : ${MUTATION_TYPE}"
fi

echo "BASE_DIR: ${BASE_DIR}"
echo "VCF_DIR  : ${VCF_DIR}"
echo "OUTPUT   : ${OUTPUT_DIR}"
echo "=================================================="

JOB_TMP="${OUTPUT_DIR}/.tmp_${PBS_JOBID}"
mkdir -p "${JOB_TMP}"

export TMPDIR="${JOB_TMP}"
export MPLCONFIGDIR="${JOB_TMP}/matplotlib"
export XDG_CACHE_HOME="${JOB_TMP}/.cache"

mkdir -p "${MPLCONFIGDIR}" "${XDG_CACHE_HOME}"

cleanup() {
    echo
    echo "Cleaning temporary files..."
    rm -rf "${JOB_TMP}"
    echo "Temporary files removed."
}
trap cleanup EXIT

export PROJECT_NAME
export BASE_DIR
export VCF_DIR
export OUTPUT_DIR

if [[ -z "${MUTATION_TYPE:-}" ]]; then
    echo
    echo "Generating mutation matrices..."
    echo

    python "${MATRIX_SCRIPT}"

    echo
    echo "Submitting SBS job..."
    qsub -N "SigProf_${PROJECT_NAME}_SBS" \
        -v PROJECT_NAME="${PROJECT_NAME}",BASE_DIR="${BASE_DIR}",MUTATION_TYPE=SBS \
        "${PBS_SCRIPT}"

    echo "Submitting DBS job..."
    qsub -N "SigProf_${PROJECT_NAME}_DBS" \
        -v PROJECT_NAME="${PROJECT_NAME}",BASE_DIR="${BASE_DIR}",MUTATION_TYPE=DBS \
        "${PBS_SCRIPT}"

    echo "Submitting ID job..."
    qsub -N "SigProf_${PROJECT_NAME}_ID" \
        -v PROJECT_NAME="${PROJECT_NAME}",BASE_DIR="${BASE_DIR}",MUTATION_TYPE=ID \
        "${PBS_SCRIPT}"

    echo
    echo "All extraction jobs submitted."

else
    echo
    echo "Running ${MUTATION_TYPE} extraction..."
    echo

    python "${EXTRACT_SCRIPT}"
fi

echo
echo "Finished : $(date)"
echo "Log file : ${LOG_FILE}"

deactivate



