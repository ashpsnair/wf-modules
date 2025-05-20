module load gcc
module load python/3.12.1-gcc11
#pip install --upgrade awscli

##downloading data to nscc
aws s3 sync s3://ashok-lab-nccs-data/ /home/users/nus/ash.ps/scratch/LR-NCCS/ --exact-timestamps

### uploading data from NSCC to AWS
aws s3 sync /home/users/nus/ash.ps/scratch/LR-NCCS/ s3://ashok-lab-nccs-data/ --exact-timestamps

##getting genomes
aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg38/ ./GATK/hg38/

#get the size of bucket
aws s3 ls s3://ashok-lab-nccs-data/ --human-readable --summarize

#### Filtering files that are in standard storage only and saving them in text file
aws s3api list-objects --bucket ashok-lab-nccs-data --query "Contents[?StorageClass=='STANDARD'].{Key: Key}" --output text > standard_files.txt

#### copying the files that are in the the text files
aws s3 sync s3://ashok-lab-nccs-data <destination> --exclude "*" --include "$(cat standard_files.txt | tr '\n' ' ')"

########### Getting list of files/ grep (patting pattersnß)

#!/bin/bash
module load gcc
module load python/3.12.1-gcc11

ids=(
    WHL383 WHT190 WHT215 WHT221 WHT227 WHT273 WHT391
    WHT412 WHT462 WHT478 WHT533 WHT750 WHT798 WHT807
    WHT688 WHT694 WHT855 WHT849 WHT858 WHT917 WHT928
)

echo "ID,File Count,File Names" > s3_fastq_data.csv

for id in "${ids[@]}"; do
    files=$(aws s3 ls s3://ashok-lab-nccs-data | grep "$id" | grep "\.gz$" | awk '{print $4}')
    
    file_count=$(echo "$files" | grep -v '^$' | wc -l)
    file_names=$(echo "$files" | tr '\n' ',' | sed 's/,$//')
    
    echo "$id,$file_count,\"$file_names\"" >> s3_fastq_data.csv
done

echo "CSV file 's3_fastq_data.csv' has been created with the requested information."

######### copying from online server to local machine

rsync -a ash.ps@aspire2a.nus.edu.sg:/home/users/nus/ash.ps/scratch/NCCS-MASH/FINAL/intersect/true-normal .

aws s3 sync s3://mirxes-scome/scWGS_20241218_191/ /home/users/nus/ash.ps/scratch/scDNA/ --exact-timestamps
 