#!/bin/bash
#PBS -l select=1:ncpus=128:mem=128gb
#PBS -l walltime=01:00:00
#PBS -P 11003581
#PBS -N run-t2dm-mutsig
#PBS -j oe

# Set variables
base_dir="/home/users/nus/ash.ps/scratch/T2DM/PD-analysis/pop-filter-vcfs"
sample_name="PD"

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

# Define log file
LOG_FILE="sigprofiler_run_$(date +%Y%m%d_%H%M%S).log"

# Start logging
{
    echo "Job started at $(date)"
    echo "Working directory: $PBS_O_WORKDIR"

    module load miniforge3
    conda activate /home/project/11003581/conda-envs/sigprofile

    echo "Running matrix generation..."
    python matrix_generation.py "$base_dir" "$sample_name"

    echo "Matrix generation completed. Starting parallel processes..."

    # Run signature extraction and R script in parallel
    python sigprofile_extraction.py "$base_dir" "$sample_name" > signature_extraction.log 2>&1 &
    Rscript r_analysis.R "$base_dir" "$sample_name" > r_analysis.log 2>&1 &

    # Wait for both processes to complete
    wait

    echo "Signature extraction and R analysis completed. Starting signature assignment..."

    # Run signature assignment
    python signature_assignment.py "$base_dir" "$sample_name" > signature_assignment.log 2>&1

    echo "All processes completed."
    echo "Job finished at $(date)"

} 2>&1 | tee -a "$LOG_FILE"

echo "Log file created: $LOG_FILE"
