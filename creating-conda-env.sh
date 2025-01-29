'''
module load miniforge3
conda create --prefix /home/project/11003581/conda-envs/mut-profile python=3.12 -y

conda activate /home/project/11003581/conda-envs/mut-profile

conda install r-base r-devtools r-reticulate -c conda-forge -y
conda install bioconda::svviz -y

wget https://github.com/Illumina/Nirvana/releases/download/v3.18.1/Nirvana-3.18.1-net6.0.zip
'''


/home/project/11003581/Tools/Nirvana-v3.18.1/Nirvana

dotnet bin/Release/net6.0/Nirvana.dll \
     -c Data/Cache/GRCh37 \
     --sd Data/SupplementaryAnnotation/GRCh37 \
     -r Data/References/Homo_sapiens.GRCh37.Nirvana.dat \
     -i HiSeq.10000.vcf.gz \
     -o HiSeq.10000