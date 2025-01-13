###### Demultiplexing only the read 1 ######

#!/bin/bash

#PBS -l select=1:ngpus=1:ncpus=64
#PBS -l walltime=15:00:00
#PBS -P 11003581
#PBS -N demultiplex
#PBS -j oe

module load miniforge3
conda activate /home/project/11003581/conda-envs/mimseq


input_dir="/home/project/11003581/Data/final-JQQ-analysis/fastqs"
output_dir="/home/project/11003581/Data/final-JQQ-analysis/demultiplexed"
barcode_file="/home/project/11003581/Data/final-JQQ-analysis/barcodes.fa"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

for i in "$input_dir"/*.gz; do
    fn=$(basename "$i" .gz)
    cutadapt --no-indels -q 30,30 --trimmed-only -j 10 -a file:$barcode_file -m 10 \
    -o "${output_dir}/${fn}_{name}.fq.gz" "$i" \
    1> "${output_dir}/${fn}_log.txt"
done

######### Change the name of the files after demultiplexing ##########
#!/bin/bash

# Define the input directory
input_dir="/home/project/11003581/Data/final-JQQ-analysis/demultiplexed"

# Change to the input directory
cd "$input_dir" || exit

# Loop through the specified file patterns
for file in *_WT.fq.gz *_DOX5d.fq.gz *_DOX20d.fq.gz; do
    # Check if the file exists to avoid errors
    [ -e "$file" ] || continue

    # Extract the Rep number and the suffix (_DOX5d or _DOX20d)
    if [[ $file =~ (Rep[1-5])_.*(_DOX5d|_DOX20d|_WT)\.fq\.gz ]]; then
        rep="${BASH_REMATCH[1]}"
        suffix="${BASH_REMATCH[2]}"
        
        # Create the new filename
        new_file="${rep}${suffix}.fq.gz"

        # Rename the file
        mv -- "$file" "$new_file"
        echo "Renamed: $file -> $new_file"
    fi
done



################### Shell command to generate the sample.txt #######################
ls /home/project/11003581/Data/final-JQQ-analysis/demultiplexed/*.fq.gz | \
awk '{
    file=$0;
    if (file ~ /WT/) {
        type="WT";
    } else if (file ~ /DOX20d/) {
        type="DOX20d";
    } else if (file ~ /DOX5d/) {
        type="DOX5d";
    } else {
        type="Unknown";
    }
    print file "\t" type;
}' > sample.txt

###############################################################################################
################### MIM SEQ script  #######################
###############################################################################################
#!/bin/bash

#PBS -l select=1:ngpus=1:ncpus=128
#PBS -l walltime=15:00:00
#PBS -P 11003581
#PBS -N run-mimseq
#PBS -j oe

module load miniforge3
conda activate /home/project/11003581/conda-envs/mimseq
module load r/4.2.0

cd /home/project/11003581/Data/final-JQQ-analysis

mimseq --species Hsap --cluster-id 0.97 --threads 128 --min-cov 0.0005 \
    --max-mismatches 0.075 --control-condition WT \
    -n WTvsDox --out-dir WTvsDox \
    --max-multi 4 --remap --remap-mismatches 0.05 sample.txt


