library(data.table)
library(dplyr)
library(ggplot2)
library(stringr)

expr <- readRDS("results/expr_tpm_matched.rds")
maf  <- readRDS("results/maf_matched.rds")
pheno <- readRDS("results/pheno_matched.rds")
surv  <- readRDS("results/surv_matched.rds")
hrd  <- fread("results/hrd/hrd_scores.tsv")

# ----------------------------
# choose the HRD-sum column
# ----------------------------
score_col <- colnames(hrd)[str_detect(tolower(colnames(hrd)), "hrd.sum|hrd_score|hrdsum|hrd\\.sum")][1]
if (is.na(score_col)) {
  stop("Could not find HRD-sum column in scarHRD output. Inspect colnames(hrd).")
}

hrd <- hrd %>%
  rename(HRD_sum = !!score_col) %>%
  mutate(patient_id = patient_id %||% .data[[colnames(hrd)[str_detect(tolower(colnames(hrd)), "patient|sample")][1]]])

# ----------------------------
# BRCA1/2 mutation label
# ----------------------------
maf_cols <- colnames(maf)
sample_col_maf <- maf_cols[str_detect(tolower(maf_cols), "tumor_sample_barcode|sample")][1]
gene_col_maf   <- maf_cols[str_detect(tolower(maf_cols), "hugo_symbol|gene")][1]
class_col_maf  <- maf_cols[str_detect(tolower(maf_cols), "variant_classification")][1]

non_silent <- c(
  "Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins",
  "Splice_Site", "Translation_Start_Site", "Nonstop_Mutation", "In_Frame_Del",
  "In_Frame_Ins", "De_novo_Start_OutOfFrame", "De_novo_Start_InFrame"
)

brca_patients <- maf %>%
  filter(.data[[gene_col_maf]] %in% c("BRCA1", "BRCA2")) %>%
  filter(is.na(class_col_maf) | .data[[class_col_maf]] %in% non_silent) %>%
  mutate(patient_id = substr(.data[[sample_col_maf]], 1, 12)) %>%
  distinct(patient_id) %>%
  pull(patient_id)

label_df <- hrd %>%
  mutate(
    BRCA_status = ifelse(patient_id %in% brca_patients, "BRCA_mut", "BRCA_WT")
  )

# pilot-friendly HRD grouping: top vs bottom quartile
q1 <- quantile(label_df$HRD_sum, 0.25, na.rm = TRUE)
q3 <- quantile(label_df$HRD_sum, 0.75, na.rm = TRUE)

label_df <- label_df %>%
  mutate(
    HRD_group = case_when(
      HRD_sum >= q3 ~ "HRD_high",
      HRD_sum <= q1 ~ "HRD_low",
      TRUE ~ "Intermediate"
    )
  )

write.table(label_df, "results/hrd/hrd_labels.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

# quick QC plot
p <- ggplot(label_df, aes(BRCA_status, HRD_sum, fill = BRCA_status)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  theme_bw(base_size = 12) +
  labs(title = "HRD score by BRCA status in TCGA-OV")

ggsave("results/hrd/hrd_by_brca_violin.png", p, width = 6, height = 4, dpi = 300)