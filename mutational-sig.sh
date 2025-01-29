from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
from SigProfilerExtractor import sigpro as sig
from SigProfilerAssignment import Analyzer as Analyze

sample_name= "PD"
base_dir="/Users/ash/NUS Dropbox/Aiswarya PS/Prof/T2DM/PD-DM/pop-filter-vcfs"

input_matrix= base_dir+ "/output/SBS/"+ sample_name+ ".SBS96.all"

matrices = matGen.SigProfilerMatrixGeneratorFunc(sample_name, "GRCh38", base_dir, plot=True, exome=False, bed_file=None, chrom_based=False, tsb_stat=False, seqInfo=False, cushion=100)

sig.sigProfilerExtractor("matrix", base_dir+"/sigpro_extract", input_matrix, reference_genome="GRCh38", opportunity_genome="GRCh38", minimum_signatures=1, maximum_signatures=10, nmf_replicates=100)

######## SIGMINER #########


library(sigminer)
library(tibble)

sample_name= "PD"
base_dir="/Users/ash/NUS Dropbox/Aiswarya PS/Prof/T2DM/PD/pop-filter-vcfs"

input_matrix= paste0(base_dir, "/output/SBS/", sample_name, ".SBS96.all")

matrix = read.delim(input_matrix, sep = "\t")
output = paste0(base_dir,"/sigminer/")

dir.create(output)

# Formats the matrix for sigminer if necessary
matrix <- as.matrix(setNames(data.frame(t(matrix[,-1])), matrix[,1]))

# Signature extraction using Bayesian NMF from K=1 to K=10
model <- sig_unify_extract(
  matrix,
  range = 1:10,
  approach = "bayes_nmf",
  nrun = 10
)


# Save H and W matrices
signatures <- as.data.frame(model$Signature.norm)
signatures <- cbind(MutationsType = rownames(signatures), signatures)
write.table(signatures, paste0(output, "denovo_sig.tsv"), sep = "\t", quote = FALSE, row.name = FALSE)
exposure <- tibble::rownames_to_column(as.data.frame(t(model$Exposure)), "Sample")
write.table(exposure, paste0(output, "denovo_exposure.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Save run logs
write.table(model$Raw$summary_run, paste0(output, "run.log"), sep = "\t", quote = FALSE, row.name = FALSE)

# [optional] Visualize de novo signatures similarity to COSMIC
library(pheatmap)
library(repr)

sim <- get_sig_similarity(model, sig_db = "latest_SBS_GRCh38")

# Save heatmap as SVG without defining aspect ratio
png(paste0(output, "similarity_heatmap.png"), width = 1200, height = 350)   

pheatmap(sim$similarity,
         color = colorRampPalette(c("#348ABD", "white", "#E24A33"))(70),
         border_color = FALSE)

dev.off()  # Close the png device


############# Sigprofile Assigned ###########

signatures = base_dir+ '/sigminer/' +'/denovo_sig.tsv'
activities = base_dir+ '/sigminer/' + '/denovo_exposure.tsv'
output = base_dir+ '/sigminer/'

Analyze.decompose_fit(samples=input_matrix,output=output,signatures=signatures,genome_build="GRCh38",signature_database='/home/pawan/anaconda3/lib/python3.8/site-packages/SigProfilerAssignment/data/Reference_Signatures/GRCh38/COSMIC_v3.3_SBS_GRCh38.txt')
