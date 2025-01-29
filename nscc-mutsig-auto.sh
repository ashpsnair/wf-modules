#!/bin/bash
#PBS -l select=1:ncpus=128:mem=128gb
#PBS -l walltime=01:00:00
#PBS -P 11003581
#PBS -N run-t2dm-mutsig
#PBS -j oe

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

    # Create a temporary Python script
    cat << EOF > run_sigprofiler.py

import os
import sys
import traceback
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
from SigProfilerExtractor import sigpro as sig

def main():
    sample_name = "PD"
    base_dir = "/home/users/nus/ash.ps/scratch/T2DM/PD-analysis/pop-filter-vcfs"
    input_matrix = os.path.join(base_dir, "output", "SBS", f"{sample_name}.SBS96.all")

    print("Starting matrix generation...")
    matrices = matGen.SigProfilerMatrixGeneratorFunc(
        sample_name, "GRCh38", base_dir, plot=True, exome=False,
        bed_file=None, chrom_based=False, tsb_stat=False, seqInfo=False,
        cushion=100
    )
    print("Matrix generation completed.")

    print("Starting signature extraction...")
    sig.sigProfilerExtractor(
        "matrix", os.path.join(base_dir, "sigpro_extract"), input_matrix,
        reference_genome="GRCh38", opportunity_genome="GRCh38",
        minimum_signatures=1, maximum_signatures=10, nmf_replicates=100,
        cpu=1
    )
    print("Signature extraction completed.")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print("An error occurred:")
        print(traceback.format_exc())
        sys.exit(1)

EOF

    echo "Running Python script..."
    python run_sigprofiler.py

    echo "Cleaning up temporary Python script..."
    rm run_sigprofiler.py

    echo "Job finished at $(date)"

} 2>&1 | tee -a "$LOG_FILE"

echo "Log file created: $LOG_FILE"




####################### trying it on YS-normal data with multiprocessing ##################
#!/bin/bash
#PBS -l select=1:ncpus=128:mem=128gb
#PBS -l walltime=01:00:00
#PBS -P 11003581
#PBS -N run-ys-normal-mutsig
#PBS -j oe

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

    # Create a temporary Python script
    cat << EOF > run_sigprofiler.py

import os
import sys
import traceback
import multiprocessing
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
from SigProfilerExtractor import sigpro as sig

def main():
    # Set the multiprocessing start method to 'spawn'
    multiprocessing.set_start_method('spawn', force=True)

    sample_name = "YS-normal"
    base_dir = "/home/users/nus/ash.ps/scratch/YS-tumor-only/normal-analysis/pop-filter-vcfs"
    input_matrix = os.path.join(base_dir, "output", "SBS", f"{sample_name}.SBS96.all")

    print("Starting matrix generation...")
    matrices = matGen.SigProfilerMatrixGeneratorFunc(
        sample_name, "GRCh38", base_dir, plot=True, exome=False,
        bed_file=None, chrom_based=False, tsb_stat=False, seqInfo=False,
        cushion=100
    )
    print("Matrix generation completed.")

    print("Starting signature extraction...")
    sig.sigProfilerExtractor(
        "matrix", os.path.join(base_dir, "sigpro_extract"), input_matrix,
        reference_genome="GRCh38", opportunity_genome="GRCh38",
        minimum_signatures=1, maximum_signatures=10, nmf_replicates=100,
        cpu=128  # Use all available CPUs for multithreading
    )
    print("Signature extraction completed.")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print("An error occurred:")
        print(traceback.format_exc())
        sys.exit(1)


EOF

    echo "Running Python script..."
    python run_sigprofiler.py

    echo "Cleaning up temporary Python script..."
    rm run_sigprofiler.py

    echo "Job finished at $(date)"

} 2>&1 | tee -a "$LOG_FILE"

echo "Log file created: $LOG_FILE"