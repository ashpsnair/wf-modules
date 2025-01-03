'''
module load miniforge3
conda create --prefix /home/project/11003581/conda-envs/sigprofile python=3.9.12

conda activate /home/project/11003581/conda-envs/sigprofile
pip install SigProfilerExtractor

conda install r-base r-devtools r-reticulate -c conda-forge -y

python
from SigProfilerMatrixGenerator import install as genInstall
genInstall.install('GRCh38', rsync=True, bash=True)

'''

##################################

#!/bin/bash

#PBS -l select=1
#PBS -l walltime=5:00:00
#PBS -P 11003581
#PBS -N sigprofile-extractor
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load miniforge3
conda activate /home/project/11003581/conda-envs/sigprofile

# Run the Python script
python -c "

from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen

matrices = matGen.SigProfilerMatrixGeneratorFunc('WH', 'GRCh38', '/home/users/nus/ash.ps/scratch/MASH-vcfs/WH', plot=True, exome=False, bed_file=None, chrom_based=False, tsb_stat=False, seqInfo=False, cushion=100)


"

