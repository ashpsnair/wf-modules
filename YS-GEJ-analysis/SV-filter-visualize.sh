'''

cd /home/project/11003581/Tools/AnnotSV/
make PREFIX=. install
make PREFIX=. install-human-annotation
make PREFIX=. install-mouse-annotation

git clone https://github.com/mobidic/knotAnnotSV.git
cd knotAnnotSV


cpan YAML::XS Sort::Key::Natural Excel::Writer::XLSX

module load miniforge3
conda activate /home/project/11003581/conda-envs/mut-profile

'''

#!/bin/bash

#PBS -l select=1:mem=128gb
#PBS -l walltime=1:00:00
#PBS -P 11003581
#PBS -N annot-sv
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load bedtools/2.30.0
module load bcftools/1.15.1

/home/project/11003581/Tools/AnnotSV/bin/AnnotSV -SVinputFile /home/users/nus/ash.ps/scratch/YS-analysis/manta-vcfs/T01_vs_N01.manta.somatic_sv.vcf.gz \
    -annotationMode full \
    -genomeBuild GRCh38 \
    -outputDir /home/users/nus/ash.ps/scratch/YS-analysis/sv-annot/ \
    -REreport 1 \
    -snvIndelPASS 1


########### vcf2circos

'''
https://hgdownload.soe.ucsc.edu/gbdb/
bgzip -d ncbirefseqfile
sort -k1,1V -k2,2n file > sortedfile

'''


################################ Running knotsv

#!/bin/bash

#PBS -l select=1:mem=128gb
#PBS -l walltime=1:00:00
#PBS -P 11003581
#PBS -N knot
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR


perl /home/project/11003581/Tools/knotAnnotSV/knotAnnotSV.pl \
    --annotSVfile /home/users/nus/ash.ps/scratch/YS-analysis/sv-annot/T01_vs_N01.manta.somatic_sv.annotated.tsv \
    --outDir /home/users/nus/ash.ps/scratch/YS-analysis/sv-annot/knot-out/ \
    --outPrefix T01\
    --genomeBuild hg38 \
    --datatableDir /home/project/11003581/Tools/knotAnnotSV/DataTables
  