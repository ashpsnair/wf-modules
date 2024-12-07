'''
module load miniforge3
conda create --prefix /home/project/11003581/conda-envs/mutsig python=3.7

conda activate /home/project/11003581/conda-envs/mutsig
pip install SigProfilerExtractor

'''



from SigProfilerExtractor import sigpro as sig

#extract the signatures by running the following command
sig.sigProfilerExtractor("vcf", "results", "path/to/21BRCA_vcf", genome_build="GRCh38", maximum_signatures=10)

