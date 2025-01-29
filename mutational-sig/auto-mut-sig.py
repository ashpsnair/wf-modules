import os
import sys
import traceback
import multiprocessing
import subprocess
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
from SigProfilerExtractor import sigpro as sig
from SigProfilerAssignment import Analyzer as Analyze

def run_matrix_generation(base_dir, sample_name):
    print("Starting matrix generation...")
    matrices = matGen.SigProfilerMatrixGeneratorFunc(
        sample_name, "GRCh38", base_dir, plot=True, exome=False,
        bed_file=None, chrom_based=False, tsb_stat=False, seqInfo=False,
        cushion=100
    )
    print("Matrix generation completed.")

def run_signature_extraction(base_dir, sample_name):
    input_matrix = os.path.join(base_dir, "output", "SBS", f"{sample_name}.SBS96.all")
    print("Starting signature extraction...")
    sig.sigProfilerExtractor(
        "matrix", os.path.join(base_dir, "sigpro_extract"), input_matrix,
        reference_genome="GRCh38", opportunity_genome="GRCh38",
        minimum_signatures=1, maximum_signatures=10, nmf_replicates=100,
        cpu=128
    )
    print("Signature extraction completed.")

def run_r_analysis(base_dir, sample_name):
    r_script = f"""
    library(sigminer)
    library(tibble)
    library(pheatmap)

    base_dir <- "{base_dir}"
    sample_name <- "{sample_name}"

    input_matrix <- file.path(base_dir, "output", "SBS", paste0(sample_name, ".SBS96.all"))
    output <- file.path(base_dir, "sigminer")

    dir.create(output, showWarnings = FALSE)

    matrix <- read.delim(input_matrix, sep = "\t")
    matrix <- as.matrix(setNames(data.frame(t(matrix[,-1])), matrix[,1]))

    model <- sig_unify_extract(
      matrix,
      range = 1:10,
      approach = "bayes_nmf",
      nrun = 10
    )

    signatures <- as.data.frame(model$Signature.norm)
    signatures <- cbind(MutationsType = rownames(signatures), signatures)
    write.table(signatures, file.path(output, "denovo_sig.tsv"), sep = "\t", quote = FALSE, row.name = FALSE)

    exposure <- tibble::rownames_to_column(as.data.frame(t(model$Exposure)), "Sample")
    write.table(exposure, file.path(output, "denovo_exposure.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

    write.table(model$Raw$summary_run, file.path(output, "run.log"), sep = "\t", quote = FALSE, row.name = FALSE)

    sim <- get_sig_similarity(model, sig_db = "latest_SBS_GRCh38")

    png(file.path(output, "similarity_heatmap.png"), width = 1200, height = 350)   

    pheatmap(sim$similarity,
             color = colorRampPalette(c("#348ABD", "white", "#E24A33"))(70),
             border_color = FALSE)

    dev.off()

    print("R analysis completed.")
    """
    with open("temp_r_script.R", "w") as f:
        f.write(r_script)
    
    subprocess.run(["Rscript", "temp_r_script.R"], check=True)
    os.remove("temp_r_script.R")

def run_signature_assignment(base_dir, sample_name):
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

def main(base_dir, sample_name):
    multiprocessing.set_start_method('spawn', force=True)

    run_matrix_generation(base_dir, sample_name)

    # Run signature extraction and R analysis in parallel
    extraction_process = multiprocessing.Process(target=run_signature_extraction, args=(base_dir, sample_name))
    r_analysis_process = multiprocessing.Process(target=run_r_analysis, args=(base_dir, sample_name))

    extraction_process.start()
    r_analysis_process.start()

    extraction_process.join()
    r_analysis_process.join()

    run_signature_assignment(base_dir, sample_name)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <base_dir> <sample_name>")
        sys.exit(1)

    base_dir = sys.argv[1]
    sample_name = sys.argv[2]

    try:
        main(base_dir, sample_name)
    except Exception as e:
        print("An error occurred:")
        print(traceback.format_exc())
        sys.exit(1)
