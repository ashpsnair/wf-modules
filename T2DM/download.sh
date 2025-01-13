https://pmc.ncbi.nlm.nih.gov/articles/PMC11187843/
https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP513973&o=library_name_s%3Aa&s=SRR29413843,SRR29413844,SRR29413845,SRR29413846,SRR29413847,SRR29413848,SRR29413849,SRR29413850,SRR29413851,SRR29413852,SRR29413853,SRR29413854,SRR29413855,SRR29413856,SRR29413857,SRR29413858,SRR29413859,SRR29728988,SRR29728989,SRR29728990,SRR29728991,SRR29728992,SRR29728993,SRR29728994,SRR29728995,SRR29728996,SRR29728997,SRR29728998,SRR29728999,SRR29729000,SRR29729001,SRR29729002,SRR29729003,SRR29729004,SRR29729005,SRR29729006


/home/project/11003581/Tools/sratoolkit.3.1.1-ubuntu64/bin/fastq-dump SRR29728988


import subprocess

# samples correspond to the project
sra_numbers = [
    "SRR29413843", "SRR29413844", "SRR29413845", "SRR29413846", "SRR29413847", "SRR29413848", "SRR29413849", "SRR29413850", "SRR29413851", "SRR29413852", "SRR29413853", "SRR29413854", "SRR29413855", "SRR29413856", "SRR29413857", "SRR29413858", "SRR29413859", "SRR29728988", "SRR29728989", "SRR29728990", "SRR29728991", "SRR29728992", "SRR29728993", "SRR29728994", "SRR29728995", "SRR29728996", "SRR29728997", "SRR29728998", "SRR29728999", "SRR29729000", "SRR29729001", "SRR29729002", "SRR29729003", "SRR29729004", "SRR29729005", "SRR29729006"
]


# this will download the .sra files to ~/ncbi/public/sra/ (will create directory if not present)
for sra_id in sra_numbers:
    print ("Currently downloading: " + sra_id)
    prefetch = "/home/project/11003581/Tools/sratoolkit.3.1.1-ubuntu64/bin/prefetch " + sra_id
    print ("The command used was: " + prefetch)
    subprocess.call(prefetch, shell=True)

# this will extract the .sra files from above into a folder named 'fastq'
for sra_id in sra_numbers:
    print ("Generating fastq for: " + sra_id)
    fastq_dump = "/home/project/11003581/Tools/sratoolkit.3.1.1-ubuntu64/bin/fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip ~/ncbi/public/sra/" + sra_id + ".sra"
    print ("The command used was: " + fastq_dump)
    subprocess.call(fastq_dump, shell=True)


