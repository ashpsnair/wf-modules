module load gcc
module load python/3.12.1-gcc11
pip install --upgrade awscli

##downloading data to nscc

aws s3 sync s3://ashok-lab-nccs-data/ /home/users/nus/ash.ps/scratch/LR-NCCS/ --exact-timestamps


##getting genomes

aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg38/ ./GATK/hg38/



