'''

cd /home/project/11003581/Tools/AnnotSV/
make PREFIX=. install
make PREFIX=. install-human-annotation
make PREFIX=. install-mouse-annotation


'''

#!/bin/bash

#PBS -l select=1:mem=128gb
#PBS -l walltime=1:00:00
#PBS -P 11003581
#PBS -N annot-sv
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load bedtools
module load bcftools

/home/project/11003581/Tools/AnnotSV/bin/AnnotSV -SVinputFile /home/users/nus/ash.ps/scratch/YS-analysis/manta-vcfs/T01_vs_N01.manta.somatic_sv.vcf.gz \
    -annotationsDir /home/project/11003581/Tools/AnnotSV/annotation_dir/Annotations_Human/ \
    -annotationMode full \
    -genomeBuild GRCh38 \
    -outputDir /home/users/nus/ash.ps/scratch/YS-analysis/sv-annot/ \
    -REreport 1 \
    -snvIndelPASS 1 \


