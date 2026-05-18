#!/bin/bash
#PBS -N md5_compare
#PBS -P 11003581
#PBS -l select=1:ncpus=128:mem=128g
#PBS -l walltime=6:00:00
#PBS -j oe

set -euo pipefail

BASE_DIR="/home/project/11003581/Data/sep-function/WT_12May/X401SC26043443-Z01-F001"

cd "${BASE_DIR}"

OUTPUT="/home/project/11003581/Data/sep-function/WT_12May/md5_comparison.tsv"

echo -e "filename\tmd5_source\tmd5_destination\tmatch" > ${OUTPUT}

while read -r expected_md5 filepath
do
    actual_md5=$(md5sum "${filepath}" | awk '{print $1}')

    if [[ "${expected_md5}" == "${actual_md5}" ]]; then
        status="MATCH"
    else
        status="FAILED"
    fi

    echo -e "${filepath}\t${expected_md5}\t${actual_md5}\t${status}" >> ${OUTPUT}

done < MD5.txt

echo "====================================="
echo "MD5 verification completed"
date
echo "====================================="

echo ""
echo "Summary:"
grep -c "MATCH" ${OUTPUT} | awk '{print "Matched files:", $1-1}'
grep -c "FAILED" ${OUTPUT} | awk '{print "Failed files:", $1}'

echo ""
echo "Output table:"
echo "${BASE_DIR}/${OUTPUT}"