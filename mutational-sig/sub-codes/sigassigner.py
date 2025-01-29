import sys
from SigProfilerAssignment import Analyzer as Analyze

def main(base_dir, sample_name):
    input_matrix = f"{base_dir}/output/SBS/{sample_name}.SBS96.all"
    signatures = f"{base_dir}/sigminer/denovo_sig.tsv"
    activities = f"{base_dir}/sigminer/denovo_exposure.tsv"
    output = f"{base_dir}/sigminer/"

    Analyze.decompose_fit(
        samples=input_matrix,
        output=output,
        signatures=signatures,
        genome_build="GRCh38",
        signature_database='/Users/ash/opt/anaconda3/lib/python3.12/site-packages/SigProfilerAssignment/data/Reference_Signatures/GRCh38/COSMIC_v3.3_SBS_GRCh38.txt'
    )

    print("Signature assignment completed.")

if __name__ == "__main__":
    base_dir = sys.argv[1]
    sample_name = sys.argv[2]
    main(base_dir, sample_name)
