##################### --------------------------------------------- #####################
# BRCA1_R133C
##################### --------------------------------------------- #####################

#samplesheet
patient,sex,status,sample,lane,fastq_1,fastq_2
BRCA1_R133C,XX,1,B_68_01,1,/home/users/nus/ash.ps/scratch/sep-function/B_68_01/B_68_01_DKDN250027426-1A_22NT5GLT4_L7_1.fq.gz,/home/users/nus/ash.ps/scratch/sep-function/B_68_01/B_68_01_DKDN250027426-1A_22NT5GLT4_L7_2.fq.gz
BRCA1_R133C,XX,1,B_68_02,1,/home/users/nus/ash.ps/scratch/sep-function/B_68_02/B_68_02_DKDN250027427-1A_22NT5GLT4_L7_1.fq.gz,/home/users/nus/ash.ps/scratch/sep-function/B_68_02/B_68_02_DKDN250027427-1A_22NT5GLT4_L7_2.fq.gz
BRCA1_R133C,XX,1,B_68_03,1,/home/users/nus/ash.ps/scratch/sep-function/B_68_03/B_68_03_DKDN250027428-1A_22NT5GLT4_L7_1.fq.gz,/home/users/nus/ash.ps/scratch/sep-function/B_68_03/B_68_03_DKDN250027428-1A_22NT5GLT4_L7_2.fq.gz
BRCA1_R133C,XX,1,B_68_04,1,/home/users/nus/ash.ps/scratch/sep-function/B_68_04/B_68_04_DKDN250027429-1A_22NT5GLT4_L7_1.fq.gz,/home/users/nus/ash.ps/scratch/sep-function/B_68_04/B_68_04_DKDN250027429-1A_22NT5GLT4_L7_2.fq.gz
BRCA1_R133C,XX,1,B_68_05,1,/home/users/nus/ash.ps/scratch/sep-function/B_68_05/B_68_05_DKDN250027430-1A_22NT5GLT4_L8_1.fq.gz,/home/users/nus/ash.ps/scratch/sep-function/B_68_05/B_68_05_DKDN250027430-1A_22NT5GLT4_L8_2.fq.gz
BRCA1_R133C,XX,1,B_68_06,1,/home/users/nus/ash.ps/scratch/sep-function/B_68_06/B_68_06_DKDN250027431-1A_22NT5GLT4_L8_1.fq.gz,/home/users/nus/ash.ps/scratch/sep-function/B_68_06/B_68_06_DKDN250027431-1A_22NT5GLT4_L8_2.fq.gz
BRCA1_R133C,XX,1,B_68_07,1,/home/users/nus/ash.ps/scratch/sep-function/B_68_07/B_68_07_DKDN250027432-1A_22NT5GLT4_L8_1.fq.gz,/home/users/nus/ash.ps/scratch/sep-function/B_68_07/B_68_07_DKDN250027432-1A_22NT5GLT4_L8_2.fq.gz
BRCA1_R133C,XX,1,B_68_08,1,/home/users/nus/ash.ps/scratch/sep-function/B_68_08/B_68_08_DKDN250027433-1A_22NT5GLT4_L8_1.fq.gz,/home/users/nus/ash.ps/scratch/sep-function/B_68_08/B_68_08_DKDN250027433-1A_22NT5GLT4_L8_2.fq.gz
BRCA1_R133C,XX,1,B_68_09,1,/home/users/nus/ash.ps/scratch/sep-function/B_68_09/B_68_09_DKDN250027434-1A_22NT5GLT4_L8_1.fq.gz,/home/users/nus/ash.ps/scratch/sep-function/B_68_09/B_68_09_DKDN250027434-1A_22NT5GLT4_L8_2.fq.gz
BRCA1_R133C,XX,1,B_68_10,1,/home/users/nus/ash.ps/scratch/sep-function/B_68_10/B_68_10_DKDN250027435-1A_22NT5GLT4_L8_1.fq.gz,/home/users/nus/ash.ps/scratch/sep-function/B_68_10/B_68_10_DKDN250027435-1A_22NT5GLT4_L8_2.fq.gz

### nf-core sarek code
#!/bin/bash

#PBS -l select=1:ncpus=128:mem=256g
#PBS -l walltime=24:00:00
#PBS -P 11003581
#PBS -N run-BRCA1_R133C
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load nextflow/24.10.5
module load miniforge3

export NXF_CONDA_CACHEDIR=/home/project/11003581/Tools/nf-conda-cache
export NXF_CONDA_USE_MAMBA=true

WORKDIR=/home/users/nus/ash.ps/scratch/sep-fun-analysis/BRCA1_R133C

nextflow run nf-core/sarek -r 3.5.1 \
   -profile conda \
   --input $WORKDIR/samplesheet.csv \
   --outdir $WORKDIR \
   --genome GATK.GRCh38



##################### --------------------------------------------- #####################
# H2A_KO
##################### --------------------------------------------- #####################

#samplesheet
patient,sex,status,sample,lane,fastq_1,fastq_2
H2A_KO,XX,1,H_16_01,1,/home/users/nus/ash.ps/scratch/sep-function/H_16_01/H_16_01_DKDN250027406-1A_22NT5GLT4_L4_1.fq.gz,/home/users/nus/ash.ps/scratch/sep-function/H_16_01/H_16_01_DKDN250027406-1A_22NT5GLT4_L4_2.fq.gz
H2A_KO,XX,1,H_16_02,1,/home/users/nus/ash.ps/scratch/sep-function/H_16_02/H_16_02_DKDN250027407-1A_22NT5GLT4_L4_1.fq.gz,/home/users/nus/ash.ps/scratch/sep-function/H_16_02/H_16_02_DKDN250027407-1A_22NT5GLT4_L4_2.fq.gz
H2A_KO,XX,1,H_16_03,1,/home/users/nus/ash.ps/scratch/sep-function/H_16_03/H_16_03_DKDN250027408-1A_22NT5GLT4_L4_1.fq.gz,/home/users/nus/ash.ps/scratch/sep-function/H_16_03/H_16_03_DKDN250027408-1A_22NT5GLT4_L4_2.fq.gz
H2A_KO,XX,1,H_16_04,1,/home/users/nus/ash.ps/scratch/sep-function/H_16_04/H_16_04_DKDN250027409-1A_22NT5GLT4_L5_1.fq.gz,/home/users/nus/ash.ps/scratch/sep-function/H_16_04/H_16_04_DKDN250027409-1A_22NT5GLT4_L5_2.fq.gz
H2A_KO,XX,1,H_16_05,1,/home/users/nus/ash.ps/scratch/sep-function/H_16_05/H_16_05_DKDN250027410-1A_22NT5GLT4_L5_1.fq.gz,/home/users/nus/ash.ps/scratch/sep-function/H_16_05/H_16_05_DKDN250027410-1A_22NT5GLT4_L5_2.fq.gz
H2A_KO,XX,1,H_16_06,1,/home/users/nus/ash.ps/scratch/sep-function/H_16_06/H_16_06_DKDN250027411-1A_22NT5GLT4_L5_1.fq.gz,/home/users/nus/ash.ps/scratch/sep-function/H_16_06/H_16_06_DKDN250027411-1A_22NT5GLT4_L5_2.fq.gz
H2A_KO,XX,1,H_16_07,1,/home/users/nus/ash.ps/scratch/sep-function/H_16_07/H_16_07_DKDN250027412-1A_22NT5GLT4_L5_1.fq.gz,/home/users/nus/ash.ps/scratch/sep-function/H_16_07/H_16_07_DKDN250027412-1A_22NT5GLT4_L5_2.fq.gz
H2A_KO,XX,1,H_16_08,1,/home/users/nus/ash.ps/scratch/sep-function/H_16_08/H_16_08_DKDN250027413-1A_22NT5GLT4_L5_1.fq.gz,/home/users/nus/ash.ps/scratch/sep-function/H_16_08/H_16_08_DKDN250027413-1A_22NT5GLT4_L5_2.fq.gz
H2A_KO,XX,1,H_16_09,1,/home/users/nus/ash.ps/scratch/sep-function/H_16_09/H_16_09_DKDN250027414-1A_22NT5GLT4_L5_1.fq.gz,/home/users/nus/ash.ps/scratch/sep-function/H_16_09/H_16_09_DKDN250027414-1A_22NT5GLT4_L5_2.fq.gz
H2A_KO,XX,1,H_16_10,1,/home/users/nus/ash.ps/scratch/sep-function/H_16_10/H_16_10_DKDN250027415-1A_22NT5GLT4_L5_1.fq.gz,/home/users/nus/ash.ps/scratch/sep-function/H_16_10/H_16_10_DKDN250027415-1A_22NT5GLT4_L5_2.fq.gz

### nf-core sarek code
#!/bin/bash

#PBS -l select=1:ncpus=128:mem=256g
#PBS -l walltime=24:00:00
#PBS -P 11003581
#PBS -N run-H2A-KO
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load nextflow/24.10.5
module load miniforge3

export NXF_CONDA_CACHEDIR=/home/project/11003581/Tools/nf-conda-cache
export NXF_CONDA_USE_MAMBA=true

WORKDIR=/home/users/nus/ash.ps/scratch/sep-fun-analysis/H2A_KO

nextflow run nf-core/sarek -r 3.5.1 \
   -profile conda \
   --input $WORKDIR/samplesheet.csv \
   --outdir $WORKDIR \
   --genome GATK.GRCh38


##################### --------------------------------------------- #####################
# BRCA2_P3292L
##################### --------------------------------------------- #####################

#samplesheet
patient,sex,status,sample,lane,fastq_1,fastq_2
BRCA2_P3292L,XX,1,P_15_01,1,/home/users/nus/ash.ps/scratch/sep-function/P_15_01/P_15_01_DKDN250027396-1A_22NT5GLT4_L3_1.fq.gz,/home/users/nus/ash.ps/scratch/sep-function/P_15_01/P_15_01_DKDN250027396-1A_22NT5GLT4_L3_2.fq.gz
BRCA2_P3292L,XX,1,P_15_02,1,/home/users/nus/ash.ps/scratch/sep-function/P_15_02/P_15_02_DKDN250027397-1A_22NT5GLT4_L3_1.fq.gz,/home/users/nus/ash.ps/scratch/sep-function/P_15_02/P_15_02_DKDN250027397-1A_22NT5GLT4_L3_2.fq.gz
BRCA2_P3292L,XX,1,P_15_03,1,/home/users/nus/ash.ps/scratch/sep-function/P_15_03/P_15_03_DKDN250027398-1A_22NT5GLT4_L3_1.fq.gz,/home/users/nus/ash.ps/scratch/sep-function/P_15_03/P_15_03_DKDN250027398-1A_22NT5GLT4_L3_2.fq.gz
BRCA2_P3292L,XX,1,P_15_04,1,/home/users/nus/ash.ps/scratch/sep-function/P_15_04/P_15_04_DKDN250027399-1A_22NT5GLT4_L3_1.fq.gz,/home/users/nus/ash.ps/scratch/sep-function/P_15_04/P_15_04_DKDN250027399-1A_22NT5GLT4_L3_2.fq.gz
BRCA2_P3292L,XX,1,P_15_05,1,/home/users/nus/ash.ps/scratch/sep-function/P_15_05/P_15_05_DKDN250027400-1A_22NT5GLT4_L3_1.fq.gz,/home/users/nus/ash.ps/scratch/sep-function/P_15_05/P_15_05_DKDN250027400-1A_22NT5GLT4_L3_2.fq.gz
BRCA2_P3292L,XX,1,P_15_06,1,/home/users/nus/ash.ps/scratch/sep-function/P_15_06/P_15_06_DKDN250027401-1A_22NT5GLT4_L3_1.fq.gz,/home/users/nus/ash.ps/scratch/sep-function/P_15_06/P_15_06_DKDN250027401-1A_22NT5GLT4_L3_2.fq.gz
BRCA2_P3292L,XX,1,P_15_07,1,/home/users/nus/ash.ps/scratch/sep-function/P_15_07/P_15_07_DKDN250027402-1A_22NT5GLT4_L4_1.fq.gz,/home/users/nus/ash.ps/scratch/sep-function/P_15_07/P_15_07_DKDN250027402-1A_22NT5GLT4_L4_2.fq.gz
BRCA2_P3292L,XX,1,P_15_08,1,/home/users/nus/ash.ps/scratch/sep-function/P_15_08/P_15_08_DKDN250027403-1A_22NT5GLT4_L4_1.fq.gz,/home/users/nus/ash.ps/scratch/sep-function/P_15_08/P_15_08_DKDN250027403-1A_22NT5GLT4_L4_2.fq.gz
BRCA2_P3292L,XX,1,P_15_09,1,/home/users/nus/ash.ps/scratch/sep-function/P_15_09/P_15_09_DKDN250027404-1A_22NT5GLT4_L4_1.fq.gz,/home/users/nus/ash.ps/scratch/sep-function/P_15_09/P_15_09_DKDN250027404-1A_22NT5GLT4_L4_2.fq.gz
BRCA2_P3292L,XX,1,P_15_10,1,/home/users/nus/ash.ps/scratch/sep-function/P_15_10/P_15_10_DKDN250027405-1A_22NT5GLT4_L4_1.fq.gz,/home/users/nus/ash.ps/scratch/sep-function/P_15_10/P_15_10_DKDN250027405-1A_22NT5GLT4_L4_2.fq.gz


### nf-core sarek code
#!/bin/bash

#PBS -l select=1:ncpus=128:mem=256g
#PBS -l walltime=24:00:00
#PBS -P 11003581
#PBS -N run-H2A-KO
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load nextflow/24.10.5
module load miniforge3

export NXF_CONDA_CACHEDIR=/home/project/11003581/Tools/nf-conda-cache
export NXF_CONDA_USE_MAMBA=true

WORKDIR=/home/users/nus/ash.ps/scratch/sep-fun-analysis/BRCA2_P3292L

nextflow run nf-core/sarek -r 3.5.1 \
   -profile conda \
   --input $WORKDIR/samplesheet.csv \
   --outdir $WORKDIR \
   --genome GATK.GRCh38


##################### --------------------------------------------- #####################
# RAD51_S181P
##################### --------------------------------------------- #####################

#samplesheet
patient,sex,status,sample,lane,fastq_1,fastq_2
RAD51_S181P,XX,1,R_21_01,1,/home/users/nus/ash.ps/scratch/sep-function/R_21_01/R_21_01_DKDN250027416-1A_22NT5GLT4_L6_1.fq.gz,/home/users/nus/ash.ps/scratch/sep-function/R_21_01/R_21_01_DKDN250027416-1A_22NT5GLT4_L6_2.fq.gz
RAD51_S181P,XX,1,R_21_02,1,/home/users/nus/ash.ps/scratch/sep-function/R_21_02/R_21_02_DKDN250027417-1A_22NT5GLT4_L6_1.fq.gz,/home/users/nus/ash.ps/scratch/sep-function/R_21_02/R_21_02_DKDN250027417-1A_22NT5GLT4_L6_2.fq.gz
RAD51_S181P,XX,1,R_21_03,1,/home/users/nus/ash.ps/scratch/sep-function/R_21_03/R_21_03_DKDN250027418-1A_22NT5GLT4_L6_1.fq.gz,/home/users/nus/ash.ps/scratch/sep-function/R_21_03/R_21_03_DKDN250027418-1A_22NT5GLT4_L6_2.fq.gz
RAD51_S181P,XX,1,R_21_04,1,/home/users/nus/ash.ps/scratch/sep-function/R_21_04/R_21_04_DKDN250027419-1A_22NT5GLT4_L6_1.fq.gz,/home/users/nus/ash.ps/scratch/sep-function/R_21_04/R_21_04_DKDN250027419-1A_22NT5GLT4_L6_2.fq.gz
RAD51_S181P,XX,1,R_21_05,1,/home/users/nus/ash.ps/scratch/sep-function/R_21_05/R_21_05_DKDN250027420-1A_22NT5GLT4_L6_1.fq.gz,/home/users/nus/ash.ps/scratch/sep-function/R_21_05/R_21_05_DKDN250027420-1A_22NT5GLT4_L6_2.fq.gz
RAD51_S181P,XX,1,R_21_06,1,/home/users/nus/ash.ps/scratch/sep-function/R_21_06/R_21_06_DKDN250027421-1A_22NT5GLT4_L6_1.fq.gz,/home/users/nus/ash.ps/scratch/sep-function/R_21_06/R_21_06_DKDN250027421-1A_22NT5GLT4_L6_2.fq.gz
RAD51_S181P,XX,1,R_21_07,1,/home/users/nus/ash.ps/scratch/sep-function/R_21_07/R_21_07_DKDN250027422-1A_22NT5GLT4_L6_1.fq.gz,/home/users/nus/ash.ps/scratch/sep-function/R_21_07/R_21_07_DKDN250027422-1A_22NT5GLT4_L6_2.fq.gz
RAD51_S181P,XX,1,R_21_08,1,/home/users/nus/ash.ps/scratch/sep-function/R_21_08/R_21_08_DKDN250027423-1A_22NYGMLT4_L7_1.fq.gz,/home/users/nus/ash.ps/scratch/sep-function/R_21_08/R_21_08_DKDN250027423-1A_22NYGMLT4_L7_2.fq.gz
RAD51_S181P,XX,1,R_21_08,2,/home/users/nus/ash.ps/scratch/sep-function/R_21_08/R_21_08_DKDN250027423-1A_22YYTKLT3_L3_1.fq.gz,/home/users/nus/ash.ps/scratch/sep-function/R_21_08/R_21_08_DKDN250027423-1A_22YYTKLT3_L3_2.fq.gz
RAD51_S181P,XX,1,R_21_09,1,/home/users/nus/ash.ps/scratch/sep-function/R_21_09/R_21_09_DKDN250027424-1A_22NT5GLT4_L7_1.fq.gz,/home/users/nus/ash.ps/scratch/sep-function/R_21_09/R_21_09_DKDN250027424-1A_22NT5GLT4_L7_2.fq.gz
RAD51_S181P,XX,1,R_21_10,1,/home/users/nus/ash.ps/scratch/sep-function/R_21_10/R_21_10_DKDN250027425-1A_22NT5GLT4_L7_1.fq.gz,/home/users/nus/ash.ps/scratch/sep-function/R_21_10/R_21_10_DKDN250027425-1A_22NT5GLT4_L7_2.fq.gz

#!/bin/bash

#PBS -l select=1:ncpus=128:mem=256g
#PBS -l walltime=24:00:00
#PBS -P 11003581
#PBS -N run-RAD51_S181P
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

module load nextflow/24.10.5
module load miniforge3

export NXF_CONDA_CACHEDIR=/home/project/11003581/Tools/nf-conda-cache
export NXF_CONDA_USE_MAMBA=true

WORKDIR=/home/users/nus/ash.ps/scratch/sep-fun-analysis/RAD51_S181P

nextflow run nf-core/sarek -r 3.5.1 \
   -profile conda \
   --input $WORKDIR/samplesheet.csv \
   --outdir $WORKDIR \
   --genome GATK.GRCh38













