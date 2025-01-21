https://pmc.ncbi.nlm.nih.gov/articles/PMC11187843/
https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP513973&o=library_name_s%3Aa&s=SRR29413843,SRR29413844,SRR29413845,SRR29413846,SRR29413847,SRR29413848,SRR29413849,SRR29413850,SRR29413851,SRR29413852,SRR29413853,SRR29413854,SRR29413855,SRR29413856,SRR29413857,SRR29413858,SRR29413859,SRR29728988,SRR29728989,SRR29728990,SRR29728991,SRR29728992,SRR29728993,SRR29728994,SRR29728995,SRR29728996,SRR29728997,SRR29728998,SRR29728999,SRR29729000,SRR29729001,SRR29729002,SRR29729003,SRR29729004,SRR29729005,SRR29729006

######################################################
####################### DOWNLOAD SRA files 

######################################################
#Download PD files

#!/bin/bash
#PBS -l select=1
#PBS -l walltime=24:00:00
#PBS -P 11003581
#PBS -N run-dwnld-pd
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load gcc
module load python/3.12.1-gcc11

# Create a temporary Python script
cat << EOF > download-t2dm.py

import subprocess

# samples correspond to the project
sra_numbers = [
    "SRR29413853", "SRR29413854", "SRR29728994", "SRR29728993", "SRR29728992",
    "SRR29728991", "SRR29728990", "SRR29728989", "SRR29413858", "SRR29728988",
    "SRR29729004", "SRR29729003", "SRR29413848", "SRR29413856", "SRR29413855"
]

# this will download the .sra files 
for sra_id in sra_numbers:
    print ("Currently downloading: " + sra_id)
    prefetch = "/home/project/11003581/Tools/sratoolkit.3.1.1-ubuntu64/bin/prefetch --max-size 70G " + sra_id
    print ("The command used was: " + prefetch)
    subprocess.call(prefetch, shell=True)
    
    # this will extract the .sra files from above into a folder named 'fastq'
    print ("Generating fastq for: " + sra_id)
    fastq_dump = "/home/project/11003581/Tools/sratoolkit.3.1.1-ubuntu64/bin/fastq-dump --outdir /home/users/nus/ash.ps/scratch/T2DM/PD/fastq/ --gzip --skip-technical --max-size 80g --readids --read-filter pass --dumpbase --split-3 --clip /home/users/nus/ash.ps/scratch/T2DM/PD/" + sra_id + "/" + sra_id + ".sra"
    print ("The command used was: " + fastq_dump)
    subprocess.call(fastq_dump, shell=True)




EOF

# Run the Python script
python download-t2dm.py

# Clean up the temporary Python script
rm download-t2dm.py



######################################################
#Downlod PD-DM filesß





######################################################
########## Extract PD-DM samples

#!/bin/bash
#PBS -l select=4:ncpus=64:mem=128g
#PBS -l walltime=15:00:00
#PBS -P 11003581
#PBS -N run-extract
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load gcc
module load python/3.12.1-gcc11

# Create a temporary Python script
cat << EOF > extract.py

import subprocess

# samples correspond to the project
sra_numbers = [
    "SRR29413850", "SRR29729001", "SRR29728999", "SRR29728998", 
    "SRR29728997", "SRR29728996","SRR29413847", "SRR29413852", "SRR29413851", 
    "SRR29729002"
]


# this will download the .sra files to ~/ncbi/public/sra/ (will create directory if not present)
for sra_id in sra_numbers:
    # this will extract the .sra files from above into a folder named 'fastq'
    print ("Generating fastq for: " + sra_id)
    fastq_dump = "/home/project/11003581/Tools/sratoolkit.3.1.1-ubuntu64/bin/fastq-dump --gzip --outdir /home/users/nus/ash.ps/scratch/T2DM/PD-DM/fastq/ --skip-technical --read-filter pass --dumpbase --split-3 --clip /home/users/nus/ash.ps/scratch/T2DM/PD-DM/" + sra_id + "/" + sra_id + ".sra"
    print ("The command used was: " + fastq_dump)
    subprocess.call(fastq_dump, shell=True)

EOF

# Run the Python script
python extract.py

# Clean up the temporary Python script
rm extract.py


######################################################
######## PD samples
 https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP513973&o=library_name_s%3Ad%3Bacc_s%3Aa&s=SRR29413853,SRR29413854,SRR29728994,SRR29728993,SRR29728992,SRR29728991,SRR29728990,SRR29728989,SRR29413858,SRR29728988,SRR29729004,SRR29729003,SRR29413848,SRR29413856,SRR29413855



 sra_numbers = [
    "SRR29413853", "SRR29413854", "SRR29728994", "SRR29728993", "SRR29728992",
    "SRR29728991", "SRR29728990", "SRR29728989", "SRR29413858", "SRR29728988",
    "SRR29729004", "SRR29729003"
]

######################################################
########## extracting PD samples
#!/bin/bash
#PBS -l select=1:ncpus=64:mem=128g
#PBS -l walltime=15:00:00
#PBS -P 11003581
#PBS -N run-extract-pd
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load gcc
module load python/3.12.1-gcc11

# Create a temporary Python script
cat << EOF > extract.py

import subprocess

# samples correspond to the project
sra_numbers = [
 "SRR29413853", "SRR29413854", "SRR29728994", "SRR29728993", "SRR29728992",
    "SRR29728991", "SRR29728990", "SRR29728989", "SRR29413858", "SRR29728988",
    "SRR29729004", "SRR29729003"
]


# this will download the .sra files to ~/ncbi/public/sra/ (will create directory if not present)
for sra_id in sra_numbers:
    # this will extract the .sra files from above into a folder named 'fastq'
    print ("Generating fastq for: " + sra_id)
    fastq_dump = "/home/project/11003581/Tools/sratoolkit.3.1.1-ubuntu64/bin/fastq-dump --gzip --outdir /home/users/nus/ash.ps/scratch/T2DM/PD/fastq/ --skip-technical --read-filter pass --dumpbase --split-3 --clip /home/users/nus/ash.ps/scratch/T2DM/PD/" + sra_id + "/" + sra_id + ".sra"
    print ("The command used was: " + fastq_dump)
    subprocess.call(fastq_dump, shell=True)

EOF

# Run the Python script
python extract.py

# Clean up the temporary Python script
rm extract.py