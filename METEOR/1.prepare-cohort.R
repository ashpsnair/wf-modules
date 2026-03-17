library(data.table)
library(dplyr)
library(stringr)
library(tibble)

dir.create("results", showWarnings = FALSE)
dir.create("temp", showWarnings = FALSE)

# ----------------------------
# INPUT PATHS
# ----------------------------
expr_file <- "data/OV_star_tpm.tsv"
seg_file  <- "data/OV_ascat3_segments.tsv"
maf_file  <- "data/OV_somatic_mutation.maf"
pheno_file <- "data/OV_phenotype.tsv"
surv_file  <- "data/OV_survival.tsv"

# ----------------------------
# HELPERS
# ----------------------------
tcga_patient <- function(x) substr(x, 1, 12)
tcga_sample  <- function(x) substr(x, 1, 16)
sample_type_code <- function(x) substr(x, 14, 15)

is_primary_tumor <- function(x) sample_type_code(x) == "01"

# ----------------------------
# 1. LOAD EXPRESSION
# Assumes Xena-style matrix:
# first column = gene, remaining columns = TCGA samples
# ----------------------------
expr_raw <- fread(expr_file)

gene_col <- colnames(expr_raw)[1]
expr_df <- as.data.frame(expr_raw)
rownames(expr_df) <- expr_df[[gene_col]]
expr_df[[gene_col]] <- NULL

# keep TCGA columns only
expr_df <- expr_df[, grepl("^TCGA-", colnames(expr_df)), drop = FALSE]

# keep primary tumors only
tumor_cols <- colnames(expr_df)[is_primary_tumor(colnames(expr_df))]
expr_df <- expr_df[, tumor_cols, drop = FALSE]

# remove duplicated genes by mean
expr_df <- expr_df %>%
  rownames_to_column("gene") %>%
  group_by(gene) %>%
  summarise(across(everything(), ~ mean(as.numeric(.x), na.rm = TRUE)), .groups = "drop") %>%
  as.data.frame()

rownames(expr_df) <- expr_df$gene
expr_df$gene <- NULL

expr_patients <- unique(tcga_patient(colnames(expr_df)))

# ----------------------------
# 2. LOAD SEGMENTS
# ----------------------------
seg <- fread(seg_file)
seg_cols <- colnames(seg)

sample_col_seg <- seg_cols[str_detect(tolower(seg_cols), "sample")][1]
if (is.na(sample_col_seg)) stop("Could not detect sample column in ASCAT3 segment file.")

seg <- seg %>%
  mutate(patient_id = tcga_patient(.data[[sample_col_seg]]))

seg_patients <- unique(seg$patient_id)

# ----------------------------
# 3. LOAD MAF
# ----------------------------
maf <- fread(maf_file)

maf_cols <- colnames(maf)
sample_col_maf <- maf_cols[str_detect(tolower(maf_cols), "tumor_sample_barcode|sample")][1]
gene_col_maf   <- maf_cols[str_detect(tolower(maf_cols), "hugo_symbol|gene")][1]
class_col_maf  <- maf_cols[str_detect(tolower(maf_cols), "variant_classification")][1]

if (is.na(sample_col_maf) | is.na(gene_col_maf)) {
  stop("Could not detect sample/gene columns in MAF file.")
}

maf <- maf %>%
  mutate(patient_id = tcga_patient(.data[[sample_col_maf]]))

maf_patients <- unique(maf$patient_id)

# ----------------------------
# 4. PHENOTYPE / SURVIVAL
# ----------------------------
pheno <- fread(pheno_file)
surv  <- fread(surv_file)

# Try to detect patient barcode column
pheno_id_col <- colnames(pheno)[str_detect(tolower(colnames(pheno)), "submitter|sample|barcode|patient")][1]
surv_id_col  <- colnames(surv)[str_detect(tolower(colnames(surv)), "submitter|sample|barcode|patient")][1]

if (!is.na(pheno_id_col)) {
  pheno <- pheno %>% mutate(patient_id = tcga_patient(.data[[pheno_id_col]]))
}
if (!is.na(surv_id_col)) {
  surv <- surv %>% mutate(patient_id = tcga_patient(.data[[surv_id_col]]))
}

# ----------------------------
# 5. INTERSECTION
# ----------------------------
common_patients <- Reduce(intersect, list(expr_patients, seg_patients, maf_patients))

cat("Matched patients across expression + segments + maf:", length(common_patients), "\n")

writeLines(common_patients, "results/common_patients.txt")

# subset expression to matched patients
matched_expr_cols <- colnames(expr_df)[tcga_patient(colnames(expr_df)) %in% common_patients]
expr_matched <- expr_df[, matched_expr_cols, drop = FALSE]

saveRDS(expr_matched, "results/expr_tpm_matched.rds")
saveRDS(seg %>% filter(patient_id %in% common_patients), "results/ascat3_segments_matched.rds")
saveRDS(maf %>% filter(patient_id %in% common_patients), "results/maf_matched.rds")
saveRDS(pheno %>% filter(patient_id %in% common_patients), "results/pheno_matched.rds")
saveRDS(surv %>% filter(patient_id %in% common_patients), "results/surv_matched.rds")