module load gcc
module load python/3.12.1-gcc11

####### Download PBMC 
aws s3 sync s3://ashok-lab-nccs-data/PBMC/01.RawData/ /home/users/nus/ash.ps/scratch/NCCS-MASH/blood-fastqs/ --exact-timestamps

##### Download fastqs where partial of them are standard and some are deep ######

#csv created as follows- saved as files.csv
prefix	filename	size_gb	storage_class	status
WHL383	WHL383_HS005-PE-R00056_L007_R1.fastq.gz	23.566448	STANDARD	Available
WHL383	WHL383_HS005-PE-R00056_L007_R1.fastq.gz.md5	0	STANDARD	Available
WHT190	WHT190_HS005-PE-R00098_HJN53ALXX_L004_R1_001.fastq.gz	4.172587	DEEP_ARCHIVE	Deep Archive Retrieval Requested

#script used to download the files
#!/bin/bash

# File path for the CSV
CSV_FILE="files.csv"

# Bucket name (replace with your S3 bucket name)
BUCKET_NAME="s3://ashok-lab-nccs-data"

# Temporary file to store filenames for downloading
FILES_TO_DOWNLOAD="files_to_download.txt"
FAILED_FILES="failed_files.txt"

# Extract filenames using column names
echo "Extracting filenames from $CSV_FILE..."
header=$(head -n 1 "$CSV_FILE")
filename_col=$(echo "$header" | awk -F, '{for (i=1; i<=NF; i++) if ($i == "filename") print i}')

if [[ -z $filename_col ]]; then
  echo "Error: 'filename' column not found in the CSV."
  exit 1
fi

awk -F, -v fname_col="$filename_col" 'NR > 1 {print $fname_col}' "$CSV_FILE" > "$FILES_TO_DOWNLOAD"

# Check if any files exist in the list
if [[ ! -s $FILES_TO_DOWNLOAD ]]; then
  echo "No files found in the CSV to download."
  exit 0
fi

# Initialize the failed files list
> "$FAILED_FILES"

# Iterate over each file and try downloading
echo "Downloading files from S3..."
while IFS= read -r FILENAME; do
  # Download the file from S3
  aws s3 cp "$BUCKET_NAME/$FILENAME" . --quiet

  # Check if the file download was successful
  if [[ $? -ne 0 ]]; then
    echo "❌ Failed to download '$FILENAME'. Adding to failed files list."
    echo "$FILENAME" >> "$FAILED_FILES"
  else
    echo "✅ Successfully downloaded '$FILENAME'."
  fi
done < "$FILES_TO_DOWNLOAD"

# Summary of failed downloads
if [[ -s $FAILED_FILES ]]; then
  echo "Some files failed to download. See the list of failed files in '$FAILED_FILES'."
else
  echo "All files downloaded successfully!"
fi

# Cleanup
rm -f "$FILES_TO_DOWNLOAD"

echo "Script execution completed."