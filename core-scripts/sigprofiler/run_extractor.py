#!/usr/bin/env python3

from SigProfilerExtractor import sigpro
import os
import sys
import time
from datetime import datetime


def log(msg):
    print(
        f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {msg}",
        flush=True
    )


def main():

    project = os.environ["PROJECT_NAME"]
    vcf_dir = os.environ["VCF_DIR"]
    output_dir = os.environ["OUTPUT_DIR"]
    mutation_type = os.environ["MUTATION_TYPE"]

    matrix_map = {
        "SBS": "SBS96",
        "DBS": "DBS78",
        "ID": "ID83",
    }

    if mutation_type not in matrix_map:
        sys.exit(f"Unknown mutation type: {mutation_type}")

    matrix_file = os.path.join(
        vcf_dir,
        "output",
        mutation_type,
        f"{project}.{matrix_map[mutation_type]}.all",
    )

    if not os.path.exists(matrix_file):
        log(f"{mutation_type} matrix not found.")
        sys.exit(0)

    log(f"Running {mutation_type}")

    start = time.time()

    try:

        sigpro.sigProfilerExtractor(
            "matrix",
            os.path.join(output_dir, f"SigProfiler_{mutation_type}"),
            matrix_file,
            reference_genome="GRCh38",
            opportunity_genome="GRCh38",
            minimum_signatures=1,
            maximum_signatures=10,
            nmf_replicates=100,
            get_all_signature_matrices=True,
        )

        elapsed = (time.time() - start) / 60

        log(f"{mutation_type} completed ({elapsed:.2f} min)")

    except Exception as e:

        log(f"{mutation_type} failed")

        log(str(e))

        sys.exit(1)


if __name__ == "__main__":
    main()
