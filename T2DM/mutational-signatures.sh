#!/bin/bash
#PBS -l select=1:ncpus=128:mem=128gb
#PBS -l walltime=01:00:00
#PBS -P 11003581
#PBS -N run-t2dm-mutsig
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load miniforge3
conda activate /home/project/11003581/conda-envs/sigprofile

# Create a temporary Python script
cat << EOF > run_sigprofiler.py

from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
from SigProfilerExtractor import sigpro as sig
from SigProfilerAssignment import Analyzer as Analyze

sample_name= "PD"
base_dir="/home/users/nus/ash.ps/scratch/T2DM/PD-analysis/pop-filter-vcfs"

input_matrix= base_dir+ "/output/SBS/"+ sample_name+ ".SBS96.all"

matrices = matGen.SigProfilerMatrixGeneratorFunc(sample_name, "GRCh38", base_dir, plot=True, exome=False, bed_file=None, chrom_based=False, tsb_stat=False, seqInfo=False, cushion=100)

sig.sigProfilerExtractor("matrix", base_dir+"/sigpro_extract", input_matrix, reference_genome="GRCh38", opportunity_genome="GRCh38", minimum_signatures=1, maximum_signatures=10, nmf_replicates=100)

EOF

# Run the Python script
python run_sigprofiler.py

# Clean up the temporary Python script
rm run_sigprofiler.py