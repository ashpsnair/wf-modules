library(dplyr)
library(stringr)
library(tidyr)

############################################################
# INPUT
############################################################

dir <- "~/Downloads/"

input_file  <- file.path(dir, "checkSize.xls")
output_file <- file.path(dir, "samplesheet.txt")

############################################################
# READ FILE
############################################################

df <- read.table(
  input_file,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
)

colnames(df) <- c("size", "path")

############################################################
# PARSE FASTQ INFO
############################################################

parsed <- df %>%
  mutate(
    
    # sample name
    sample = str_extract(path, "(?<=01.RawData/)[^/]+"),
    
    # patient name
    patient = str_extract(sample, "^[^_]+"),
    
    # lane
    lane = str_extract(path, "(?<=_L)\\d+"),
    
    # read pair
    read = str_extract(path, "(?<=_)\\d(?=\\.fq\\.gz)"),
    
    # flowcell
    flowcell = str_extract(path, "(?<=-1A_)[A-Z0-9]+(?=_L)"),
    
    # clean path
    fastq = str_remove(path, "^\\./"),
    
    sex = "XX",
    status = 1
  )

############################################################
# CHECK FOR MULTIPLE FLOWCELLS
############################################################

flowcell_check <- parsed %>%
  distinct(sample, lane, flowcell) %>%
  group_by(sample, lane) %>%
  summarise(
    n_flowcells = n(),
    flowcells = paste(flowcell, collapse = ", "),
    .groups = "drop"
  ) %>%
  filter(n_flowcells > 1)

if(nrow(flowcell_check) > 0){
  
  cat("\n====================================================\n")
  cat("WARNING: MULTIPLE FLOWCELLS DETECTED FOR SAME LANE\n")
  cat("====================================================\n\n")
  
  print(flowcell_check)
  
  cat("\nThese samples will be split into separate rows\n")
  cat("to prevent incorrect R1/R2 pairing in nf-core/sarek.\n\n")
}

############################################################
# CREATE SAMPLE SHEET
############################################################

samplesheet <- parsed %>%
  
  pivot_wider(
    id_cols = c(
      patient,
      sex,
      status,
      sample,
      lane,
      flowcell
    ),
    
    names_from = read,
    values_from = fastq,
    names_prefix = "fastq_"
  ) %>%
  
  filter(
    !is.na(fastq_1),
    !is.na(fastq_2)
  ) %>%
  
  arrange(
    sample,
    as.numeric(lane),
    flowcell
  ) %>%
  
  # flowcell removed from final output
  select(
    patient,
    sex,
    status,
    sample,
    lane,
    fastq_1,
    fastq_2
  )

############################################################
# WRITE OUTPUT
############################################################

write.table(
  samplesheet,
  file = output_file,
  sep = ",",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

############################################################
# DONE
############################################################

cat("\nSamplesheet written to:\n")
cat(output_file, "\n\n")

print(samplesheet)
