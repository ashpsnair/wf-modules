'''
#Installation
git clone https://github.com/AlexandrovLab/SigProfilerExtractor.git
module load gcc
module load python/3.12.1-gcc11
cd SigProfilerExtractor

pip install --upgrade pip setuptools
pip install .



python3
from SigProfilerMatrixGenerator import install as genInstall
genInstall.install('GRCh38')




'''