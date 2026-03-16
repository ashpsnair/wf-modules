source: https://github.com/pipekorsi/isomut2py


#Install, setup
cd /home/project/11003581/Tools
git clone https://github.com/pipekorsi/isomut2py.git
cd isomut2/src

cd /home/project/11003581/Tools
python3 -m venv isomut2py_env



module load gcc/11.2.0
module load python/3.7.13
module load samtools/1.22.1
source isomut2py_env/bin/activate

isomut2_path="/home/project/11003581/Tools/isomut2"



################### Test script

### python script
from isomut2py import IsoMut2py

iso = IsoMut2py(
    bam_files=[
        "/home/project/11003581/bams/sample1.sorted.bam",
        "/home/project/11003581/bams/sample2.sorted.bam"
    ],
    reference_genome="/home/project/11003581/reference/hg38.fa",
    output_dir="/home/project/11003581/analysis/isomut2py/results",

    # either this OR put isomut2 in PATH
    isomut2_path="/home/project/11003581/Tools/isomut2/isomut2",

    min_var_reads=5,
    min_freq=0.1,
    max_freq=0.9,
    base_quality_limit=20,
    mapping_quality_limit=30,
    n_threads=8
)

iso.run()


### PBS script

#!/bin/bash
#PBS -N isomut2py
#PBS -l select=1:ncpus=8:mem=32gb
#PBS -l walltime=24:00:00
#PBS -j oe

cd $PBS_O_WORKDIR

module load gcc/11.2.0
module load python/3.7.13
module load samtools/1.22.1

source /home/project/11003581/Tools/isomut2py_env/bin/activate

python run_isomut2py.py
