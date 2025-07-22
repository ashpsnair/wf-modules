module load singularity/4.3.1
module load nextflow/25.04.6
module unload java
module load java/17.0.6-jdk

#Set memory caps in your shell (recommended by nf-core):
export NXF_OPTS='-Xms1g -Xmx4g'

export SINGULARITY_BINDPATH=/home/users/nus/ash.ps/scratch,/home/project/11003581/Tools/singularity-cache

export NXF_SINGULARITY_CACHEDIR=/home/project/11003581/Tools/singularity-cache
export SINGULARITY_CACHEDIR=/home/project/11003581/Tools/singularity-cache

