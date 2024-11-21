module load gcc
module load python/3.12.1-gcc11
pip install --upgrade awscli

##downloading data to nscc

aws s3 cp s3://ashok-lab-nccs-data/ /home/users/nus/ash.ps/scratch/LR-NCCS/ --recursive


s3/ashok-lab-nccs-data