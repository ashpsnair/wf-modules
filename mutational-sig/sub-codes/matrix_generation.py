import sys
import traceback
import multiprocessing
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen

def main(base_dir, sample_name):
    multiprocessing.set_start_method('spawn', force=True)

    print("Starting matrix generation...")
    matrices = matGen.SigProfilerMatrixGeneratorFunc(
        sample_name, "GRCh38", base_dir, plot=True, exome=False,
        bed_file=None, chrom_based=False, tsb_stat=False, seqInfo=False,
        cushion=100
    )
    print("Matrix generation completed.")

if __name__ == "__main__":
    try:
        base_dir = sys.argv[1]
        sample_name = sys.argv[2]
        main(base_dir, sample_name)
    except Exception as e:
        print("An error occurred:")
        print(traceback.format_exc())
        sys.exit(1)
