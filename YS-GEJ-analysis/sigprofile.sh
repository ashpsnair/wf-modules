# Set the path to the directory
path <- "/home/users/nus/ash.ps/scratch/YS-analysis/filter-annovar/Tumor/"

# List all directories in the specified path
folder_names <- list.dirs(path, full.names = FALSE, recursive = FALSE)
# Assuming folder_names has been defined as per your previous request

for (i in folder_names) {
  # Print the processing message
  print(paste0("Processing ", i))
  
  # Construct the file path for reading the matrix
  matrix <- read.delim(paste0("~/Downloads/MASH-vcfs/patient-wise-vaf20/", i, "/sigprofile/output/SBS/", i, ".SBS96.all"), sep = "\t")
  
  # Define output directory
  output <- paste0("~/Downloads/MASH-vcfs/patient-wise-vaf20/", i, "/sigprofile/sigminer/")
  
  # Create output directory if it doesn't exist
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  
  # Format the matrix for sigminer if necessary
  matrix <- as.matrix(setNames(data.frame(t(matrix[,-1])), matrix[,1]))
  
  # Signature extraction using Bayesian NMF from K=1 to K=10
  model <- sig_unify_extract(
    matrix,
    range = 1:10,
    approach = "bayes_nmf",
    nrun = 10
  )
  
  # Open a PNG graphics device
  png(filename = paste0(output, "sig-survey.png"), width = 1000, height = 1000)  # Adjust width and height as needed
  
  # Generate the plot
  show_sig_number_survey2(model$survey)
  
  # Close the graphics device to save the file
  dev.off()
  
  # Save H and W matrices
  signatures <- as.data.frame(model$Signature.norm)
  signatures <- cbind(MutationsType = rownames(signatures), signatures)
  write.table(signatures, paste0(output, "denovo_sig.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  
  exposure <- tibble::rownames_to_column(as.data.frame(t(model$Exposure)), "Sample")
  write.table(exposure, paste0(output, "denovo_exposure.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Save run logs
  write.table(model$Raw$summary_run, paste0(output, "run.log"), sep = "\t", quote = FALSE, row.names = FALSE)
  
  # [optional] Visualize de novo signatures similarity to COSMIC
  library(pheatmap)
  
  sim <- get_sig_similarity(model, sig_db = "latest_SBS_GRCh38")
  
  # Save heatmap as PNG without defining aspect ratio
  png(paste0(output, "similarity_heatmap.png"), width = 1200, height = 350)   
  
  pheatmap(sim$similarity,
           color = colorRampPalette(c("#348ABD", "white", "#E24A33"))(70),
           border_color = FALSE)
  
  dev.off()  # Close the png device
}