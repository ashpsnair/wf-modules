#!/bin/bash
#PBS -l select=1:ncpus=128:mem=256g
#PBS -l walltime=06:00:00
#PBS -P 11003581
#PBS -N SigProfiler
#PBS -j oe
#PBS -M ash.ps@nus.edu.sg

set -euo pipefail

########################################
# Load modules
########################################

module purge
module load gcc/15.2.0-nscc
module load python/3.12.1-gcc11

########################################
# Activate SigProfiler environment
########################################

source /home/project/11003581/conda-envs/SigProfiler/bin/activate

########################################
# USER SETTINGS (EDIT THESE)
########################################

PROJECT_NAME="4HNE"

VCF_DIR="/home/users/nus/ash.ps/scratch/trial-mutsig/vcfs"

OUTPUT_DIR="/home/users/nus/ash.ps/scratch/trial-mutsig"

PYTHON_SCRIPT="/home/project/11003581/Tools/scripts/run_sigprofiler.py"

########################################
# Logging
########################################

LOG_DIR="${OUTPUT_DIR}/logs"
mkdir -p "${LOG_DIR}"

LOG_FILE="${LOG_DIR}/${PROJECT_NAME}_$(date +%Y%m%d_%H%M%S).log"

# Write everything to both terminal and log
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "=================================================="
echo "SigProfiler PBS Job"
echo "Started : $(date)"
echo "Project : ${PROJECT_NAME}"
echo "Job ID  : ${PBS_JOBID}"
echo "Node    : $(hostname)"
echo "=================================================="

########################################
# Temporary directory
########################################

JOB_TMP="${OUTPUT_DIR}/.tmp_${PBS_JOBID}"

mkdir -p "${JOB_TMP}"

export TMPDIR="${JOB_TMP}"

export MPLCONFIGDIR="${JOB_TMP}/matplotlib"
mkdir -p "${MPLCONFIGDIR}"

export XDG_CACHE_HOME="${JOB_TMP}/.cache"
mkdir -p "${XDG_CACHE_HOME}"

cleanup() {

    echo
    echo "Cleaning temporary files..."

    rm -rf "${JOB_TMP}"

    echo "Temporary files removed."

}

trap cleanup EXIT

########################################
# Export variables for Python
########################################

export PROJECT_NAME
export VCF_DIR
export OUTPUT_DIR

########################################
# Run SigProfiler
########################################

echo
echo "Running SigProfiler..."
echo

python "${PYTHON_SCRIPT}"

echo
echo "=================================================="
echo "Job completed successfully"
echo "Finished : $(date)"
echo "Log file : ${LOG_FILE}"
echo "=================================================="

########################################
# Deactivate environment
########################################

deactivate