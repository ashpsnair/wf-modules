import os
import sys
import traceback
import multiprocessing
from SigProfilerExtractor import sigpro as sig

def main(base_dir, sample_name):
    multiprocessing.set_start_method('spawn', force=True)

    input_matrix = os.path.join(base_dir, "output", "SBS", f"{sample_name}.SBS96.all")

    print("Starting signature extraction...")
    sig.sigProfilerExtractor(
        "matrix", os.path.join(base_dir, "sigpro_extract"), input_matrix,
        reference_genome="GRCh38", opportunity_genome="GRCh38",
        minimum_signatures=1, maximum_signatures=10, nmf_replicates=100,
        cpu=128
    )
    print("Signature extraction completed.")

if __name__ == "__main__":
    try:
        base_dir = sys.argv[1]
        sample_name = sys.argv[2]
        main(base_dir, sample_name)
    except Exception as e:
        print("An error occurred:")
        print(traceback.format_exc())
        sys.exit(1)
