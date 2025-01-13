module load gcc
module load python/3.12.1-gcc11
#pip install --upgrade awscli

##downloading data to nscc

aws s3 sync s3://ashok-lab-nccs-data/ /home/users/nus/ash.ps/scratch/LR-NCCS/ --exact-timestamps


##getting genomes

aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg38/ ./GATK/hg38/

#get the size of bucket
aws s3 ls s3://ashok-lab-nccs-data/ --human-readable --summarize

#### Filtering files that are in standard storage only and saving them in text file
aws s3api list-objects --bucket ashok-lab-nccs-data --query "Contents[?StorageClass=='STANDARD'].{Key: Key}" --output text > standard_files.txt

#### copying the files that are in the the text files
aws s3 sync s3://ashok-lab-nccs-data <destination> --exclude "*" --include "$(cat standard_files.txt | tr '\n' ' ')"

