### Building singulairty image in HPC
module load singularity

singularity build --sandbox myContainerSandbox docker://ubuntu


#Launch an interactive shell in your sandbox
singularity shell myContainerSandbox

#Inside the container, install your desired bioinformatics tools using package managers or by compiling from source. 
apt-get update && apt-get install blast

#When done, exit the shell (Ctrl+D) and convert your sandbox into a read-only Singularity image:
sudo singularity build myContainer.simg myContainerSandbox/