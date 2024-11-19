#source:https://bioconductor.org/packages/devel/bioc/vignettes/MutationalPatterns/inst/doc/Introduction_to_MutationalPatterns.html#:~:text=Mutational%20signatures%20are%20thought%20to,negative%20matrix%20factorization%20(NMF).
'''
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("MutationalPatterns")

'''
library(MutationalPatterns)
library(BSgenome.Hsapiens.UCSC.hg38)
ref_genome = BSgenome.Hsapiens.UCSC.hg38


vcf_files <- list.files(system.file("extdata", package = "MutationalPatterns"),
  pattern = "sample.vcf", full.names = TRUE
)

#You also need to define corresponding sample names for the VCF files:

sample_names <- c(
  "colon1", "colon2", "colon3",
  "intestine1", "intestine2", "intestine3",
  "liver1", "liver2", "liver3"
)

## ----read_vcfs_as_granges, message=FALSE------------------------------------------------------
grl <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)

## ----store tissue variable--------------------------------------------------------------------
tissue <- c(rep("colon", 3), rep("intestine", 3), rep("liver", 3))

## ----locate blood data------------------------------------------------------------------------
blood_vcf_fnames <- list.files(
  system.file("extdata", package = "MutationalPatterns"), 
  pattern = "blood.*vcf", full.names = TRUE)

## ----set blood data sample names--------------------------------------------------------------
blood_sample_names <- c("blood1", "blood2", "blood3")

## ----read blood data--------------------------------------------------------------------------
blood_grl <- read_vcfs_as_granges(blood_vcf_fnames, blood_sample_names, 
                                  ref_genome, type = "all")

## ----get mut types----------------------------------------------------------------------------
snv_grl <- get_mut_type(blood_grl, type = "snv")
indel_grl <- get_mut_type(blood_grl, type = "indel")
dbs_grl <- get_mut_type(blood_grl, type = "dbs")
mbs_grl <- get_mut_type(blood_grl, type = "mbs")

## ----read indels, eval=FALSE------------------------------------------------------------------
#  indel_grl <- read_vcfs_as_granges(blood_vcf_fnames, blood_sample_names,
#                                    ref_genome, type = "indel")

## ----read predefined dbs, eval=FALSE----------------------------------------------------------
#  predefined_dbs_grl <- read_vcfs_as_granges(blood_vcf_fnames, blood_sample_names,
#                                    ref_genome, type = "dbs",
#                                    predefined_dbs_mbs = TRUE)

## ----mutations_from_vcf-----------------------------------------------------------------------
muts <- mutations_from_vcf(grl[[1]])
head(muts, 12)

## ----mut_type---------------------------------------------------------------------------------
types <- mut_type(grl[[1]])
head(types, 12)

## ----mut_context------------------------------------------------------------------------------
context <- mut_context(grl[[1]], ref_genome)
head(context, 12)

## ----type_context-----------------------------------------------------------------------------
type_context <- type_context(grl[[1]], ref_genome)
lapply(type_context, head, 12)

## ----mut_type_occurrences---------------------------------------------------------------------
type_occurrences <- mut_type_occurrences(grl, ref_genome)
type_occurrences

## ----plot_spectrum----------------------------------------------------------------------------
p1 <- plot_spectrum(type_occurrences)

## ----plot_spectrum_2--------------------------------------------------------------------------
p2 <- plot_spectrum(type_occurrences, CT = TRUE)

## ----plot_spectrum_3--------------------------------------------------------------------------
p3 <- plot_spectrum(type_occurrences, CT = TRUE, 
                    indv_points = TRUE, legend = FALSE)

## ----combine_plot_spectrum, eval=TRUE, fig.wide = TRUE, message=FALSE-------------------------
library("gridExtra")
grid.arrange(p1, p2, p3, ncol = 3, widths = c(3, 3, 1.75))

## ----plot_spectrum_4--------------------------------------------------------------------------
p4 <- plot_spectrum(type_occurrences, by = tissue, CT = TRUE, legend = TRUE)

## ----plot_spectrum_5--------------------------------------------------------------------------
p5 <- plot_spectrum(type_occurrences, CT = TRUE, 
                    legend = TRUE, error_bars = "stdev")

## ----combine_plot_spectrum_2, fig.wide=TRUE, message=FALSE------------------------------------
grid.arrange(p4, p5, ncol = 2, widths = c(4, 2.3))

## ----mut_matrix-------------------------------------------------------------------------------
mut_mat <- mut_matrix(vcf_list = grl, ref_genome = ref_genome)
head(mut_mat)

## ----plot_96_profile, fig.wide=TRUE-----------------------------------------------------------
plot_96_profile(mut_mat[, c(1, 7)])

## ----mut_mat_extendend_context----------------------------------------------------------------
mut_mat_ext_context <- mut_matrix(grl, ref_genome, extension = 2)
head(mut_mat_ext_context)

## ----plot_profile_heatmap, fig.wide=TRUE------------------------------------------------------
plot_profile_heatmap(mut_mat_ext_context, by = tissue)

## ----riverplot, fig.wide=TRUE-----------------------------------------------------------------
plot_river(mut_mat_ext_context[,c(1,4)])

## ----get_indel_context------------------------------------------------------------------------
indel_grl <- get_indel_context(indel_grl, ref_genome)
head(indel_grl[[1]], n = 3)

## ----count_indel_contexts---------------------------------------------------------------------
indel_counts <- count_indel_contexts(indel_grl)
head(indel_counts)

## ----fig.wide=TRUE----------------------------------------------------------------------------
plot_indel_contexts(indel_counts, condensed = TRUE)

## ----fig.wide=TRUE----------------------------------------------------------------------------
plot_main_indel_contexts(indel_counts)

## ----get_DBS_contexts-------------------------------------------------------------------------
head(dbs_grl[[1]])
dbs_grl <- get_dbs_context(dbs_grl)
head(dbs_grl[[1]])

## ----count_DBS_contexts-----------------------------------------------------------------------
dbs_counts <- count_dbs_contexts(dbs_grl)

## ----plot_DBS_contexts, fig.wide=TRUE---------------------------------------------------------
plot_dbs_contexts(dbs_counts, same_y = TRUE)

## ----plot_main_DBS_contexts, fig.wide=TRUE----------------------------------------------------
plot_main_dbs_contexts(dbs_counts, same_y = TRUE)

## ----count_MBS_contexts-----------------------------------------------------------------------
mbs_counts <- count_mbs_contexts(mbs_grl)

## ----plot_MBS_contexts------------------------------------------------------------------------
plot_mbs_contexts(mbs_counts, same_y = TRUE)

## ----pool_mut_mat-----------------------------------------------------------------------------
pooled_mut_mat <- pool_mut_mat(mut_mat, grouping = tissue)
head(pooled_mut_mat)

## ----psuedo_count-----------------------------------------------------------------------------
mut_mat <- mut_mat + 0.0001

## ----use_nmf----------------------------------------------------------------------------------
library("NMF")
estimate <- nmf(mut_mat, rank = 2:5, method = "brunet", 
                nrun = 10, seed = 123456, .opt = "v-p")

## ----estimate_rank, fig.wide=TRUE-------------------------------------------------------------
plot(estimate)

## ----extract_signatures-----------------------------------------------------------------------
nmf_res <- extract_signatures(mut_mat, rank = 2, nrun = 10, single_core = TRUE)

## ----extract_signatures_snv_indel-------------------------------------------------------------
combi_mat = rbind(indel_counts, dbs_counts)
nmf_res_combi <- extract_signatures(combi_mat, rank = 2, nrun = 10, single_core = TRUE)

## ----use_bayesnmf, message=FALSE--------------------------------------------------------------
# BiocManager::install("ccfindR")
library("ccfindR")
sc <- scNMFSet(count = mut_mat)
set.seed(4)
estimate_bayes <- vb_factorize(sc, ranks = 1:3, nrun = 1, 
                               progress.bar = FALSE, verbose = 0)
plot(estimate_bayes)

## ----extract_signatures_bayes-----------------------------------------------------------------
nmf_res_bayes <- extract_signatures(mut_mat, rank = 2, nrun = 10, 
                                    nmf_type = "variational_bayes")

## ----add_column_names-------------------------------------------------------------------------
colnames(nmf_res$signatures) <- c("Signature A", "Signature B")
rownames(nmf_res$contribution) <- c("Signature A", "Signature B")

## ----read_signatures--------------------------------------------------------------------------
signatures = get_known_signatures()

## ----change_names_NMF_sigs--------------------------------------------------------------------
nmf_res <- rename_nmf_signatures(nmf_res, signatures, cutoff = 0.85)
colnames(nmf_res$signatures)

## ----plot_96_profile_signatures, fig.wide=TRUE------------------------------------------------
plot_96_profile(nmf_res$signatures, condensed = TRUE)

## ----plot_contribution------------------------------------------------------------------------
plot_contribution(nmf_res$contribution, nmf_res$signature,
  mode = "relative"
)

## ----plot_contribution_heatmap_clust, fig.wide=TRUE-------------------------------------------
plot_contribution_heatmap(nmf_res$contribution, 
                          cluster_samples = TRUE, 
                          cluster_sigs = TRUE)

## ----cluster signatures-----------------------------------------------------------------------
hclust_signatures <- cluster_signatures(nmf_res$signatures, method = "average")
signatures_order <- colnames(nmf_res$signatures)[hclust_signatures$order]
signatures_order

## ----cluster samples--------------------------------------------------------------------------
hclust_samples <- cluster_signatures(mut_mat, method = "average")
samples_order <- colnames(mut_mat)[hclust_samples$order]
samples_order

## ----plot_contribution_heatmap_order, fig.wide=TRUE-------------------------------------------
plot_contribution_heatmap(nmf_res$contribution,
  sig_order = signatures_order, sample_order = samples_order,
  cluster_sigs = FALSE, cluster_samples = FALSE
)

## ----plot_compare_profiles--------------------------------------------------------------------
plot_compare_profiles(mut_mat[, 1],
  nmf_res$reconstructed[, 1],
  profile_names = c("Original", "Reconstructed"),
  condensed = TRUE
)

## ----plot_ori_vs_rec--------------------------------------------------------------------------
plot_original_vs_reconstructed(mut_mat, nmf_res$reconstructed, 
                               y_intercept = 0.95)

## ----fit_to_signatures------------------------------------------------------------------------
fit_res <- fit_to_signatures(mut_mat, signatures)

## ----plot_contribution_refit------------------------------------------------------------------
plot_contribution(fit_res$contribution,
  coord_flip = FALSE,
  mode = "absolute"
)

## ----plot_ori_vs_rec_fit----------------------------------------------------------------------
plot_original_vs_reconstructed(mut_mat, fit_res$reconstructed, 
                               y_intercept = 0.95)

## ----fit_to_signatures_strict-----------------------------------------------------------------
strict_refit <- fit_to_signatures_strict(mut_mat, signatures, max_delta = 0.004)

## ----strict_refit_decay, fig.wide=TRUE--------------------------------------------------------
fig_list <- strict_refit$sim_decay_fig
fig_list[[1]]

## ----plot_contribution_refit_strict-----------------------------------------------------------
fit_res_strict <- strict_refit$fit_res
plot_contribution(fit_res_strict$contribution,
  coord_flip = FALSE,
  mode = "absolute"
)

## ----fit_to_signatures_strict_best_subset-----------------------------------------------------
best_subset_refit <- fit_to_signatures_strict(mut_mat, signatures[,1:5], max_delta = 0.002, method = "best_subset")

## ----Merge_signatures-------------------------------------------------------------------------
merged_signatures <- merge_signatures(signatures, cos_sim_cutoff = 0.8)

## ----fit_to_signatures_bootstrapped, message=FALSE--------------------------------------------
contri_boots <- fit_to_signatures_bootstrapped(mut_mat[, c(3, 7)],
  signatures,
  n_boots = 50,
  method = "strict"
)

## ----plot_bootstrapped_refit------------------------------------------------------------------
plot_bootstrapped_contribution(contri_boots)

## ----plot_bootstrapped_refit_rel_boxtplot-----------------------------------------------------
plot_bootstrapped_contribution(contri_boots, 
                               mode = "relative", 
                               plot_type = "dotplot")

## ----plot_correlation_signatures--------------------------------------------------------------
fig_list <- plot_correlation_bootstrap(contri_boots)
fig_list[[2]]

## ----Calculate_cossim_single------------------------------------------------------------------
cos_sim(mut_mat[, 1], signatures[, 1])

## ----Calculate_cossim-------------------------------------------------------------------------
cos_sim_samples_signatures <- cos_sim_matrix(mut_mat, signatures)
cos_sim_samples_signatures[1:3, 1:3]

## ----Plot_cosine_heatmap----------------------------------------------------------------------
plot_cosine_heatmap(cos_sim_samples_signatures, 
                    cluster_rows = TRUE, cluster_cols = TRUE)

## ----plot_cosine_heatmap_samples--------------------------------------------------------------
cos_sim_samples <- cos_sim_matrix(mut_mat, mut_mat)
plot_cosine_heatmap(cos_sim_samples, cluster_rows = TRUE, cluster_cols = TRUE)

## ----get_transcription_annotation_ojbect------------------------------------------------------
# For example get known genes table from UCSC for hg19 using
# BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
# BiocManager::install("AnnotationDbi")
# BiocManager::install("GenomicFeatures")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

## ----Choose_gene_ids--------------------------------------------------------------------------
gene_ids <- c(7157, 3845, 4893, 673, 675, 1029, 8289, 5728, 7015)

## ----context_potential_damage_analysis--------------------------------------------------------
contexts <- rownames(mut_mat)
context_mismatches <- context_potential_damage_analysis(contexts, txdb, ref_genome, gene_ids)
head(context_mismatches)

## ---------------------------------------------------------------------------------------------
sig_damage <- signature_potential_damage_analysis(signatures, contexts, context_mismatches)
head(sig_damage)

## ----load_signatures_indels-------------------------------------------------------------------
signatures_indel = get_known_signatures(muttype = "indel")
signatures_indel[1:5, 1:5]

## ----load_signatures_artifacts----------------------------------------------------------------
signatures_artifacts = get_known_signatures(incl_poss_artifacts = TRUE)
dim(signatures_artifacts)

## ----load_signatures_grch38-------------------------------------------------------------------
signatures_GRCh38 = get_known_signatures(genome = "GRCh38")
dim(signatures_GRCh38)

## ----load_signatures_signal-------------------------------------------------------------------
signatures_signal = get_known_signatures(source = "SIGNAL")
signatures_signal[1:5, 1:5]

## ----load_signatures_exposure-----------------------------------------------------------------
signatures_exposure = get_known_signatures(source = "SIGNAL", sig_type = "exposure")
signatures_exposure[1:5, 1:5]

## ----load_signatures_tissue-------------------------------------------------------------------
signatures_stomach = get_known_signatures(source = "SIGNAL", sig_type = "tissue", tissue_type = "Stomach")
signatures_stomach[1:5, 1:5]

## ----load_signatures_tissue_error, eval=FALSE-------------------------------------------------
#  get_known_signatures(source = "SIGNAL", sig_type = "tissue", tissue_type = "?")

## ----convert_signatures_to_ref_1--------------------------------------------------------------
fit_res_tissue <- fit_to_signatures(mut_mat, signatures_stomach)
fit_res_tissue$contribution[1:5, 1:5]

## ----convert_signatures_to_ref_2--------------------------------------------------------------
fit_res_tissue <- convert_sigs_to_ref(fit_res_tissue)
fit_res_tissue$contribution[1:5, 1:5]

## ----get_genes--------------------------------------------------------------------------------
genes_hg19 <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
genes_hg19

## ----mut_strand-------------------------------------------------------------------------------
strand <- mut_strand(grl[[1]], genes_hg19)
head(strand, 10)

## ----mut_matrix_stranded----------------------------------------------------------------------
mut_mat_s <- mut_matrix_stranded(grl, ref_genome, genes_hg19)
mut_mat_s[1:5, 1:5]

## ----plot_192_spectrum------------------------------------------------------------------------
plot_192_profile(mut_mat_s[, 1:2])

## ----strand_occurrences-----------------------------------------------------------------------
strand_counts <- strand_occurrences(mut_mat_s, by = tissue)
head(strand_counts)

## ----strand_bias_test-------------------------------------------------------------------------
strand_bias <- strand_bias_test(strand_counts)
head(strand_bias)

## ----plot_strand------------------------------------------------------------------------------
ps1 <- plot_strand(strand_counts, mode = "relative")

## ----plot_strand_bias-------------------------------------------------------------------------
ps2 <- plot_strand_bias(strand_bias, sig_type = "p")

## ----plot_strand_bias_combi, fig.wide=TRUE----------------------------------------------------
grid.arrange(ps1, ps2)

## ----strand_bias_test_notstrict---------------------------------------------------------------
strand_bias_notstrict <- strand_bias_test(strand_counts,
  p_cutoffs = c(0.5, 0.1, 0.05),
  fdr_cutoffs = 0.5
)
plot_strand_bias(strand_bias_notstrict, sig_type = "p")

## ----repli_file-------------------------------------------------------------------------------
repli_file <- system.file("extdata/ReplicationDirectionRegions.bed",
  package = "MutationalPatterns"
)
repli_strand <- read.table(repli_file, header = TRUE)
# Store in GRanges object
repli_strand_granges <- GRanges(
  seqnames = repli_strand$Chr,
  ranges = IRanges(
    start = repli_strand$Start + 1,
    end = repli_strand$Stop
  ),
  strand_info = repli_strand$Class
)
# UCSC seqlevelsstyle
seqlevelsStyle(repli_strand_granges) <- "UCSC"
repli_strand_granges

## ----specify_levels---------------------------------------------------------------------------
repli_strand_granges$strand_info <- factor(repli_strand_granges$strand_info,
  levels = c("right", "left")
)

## ----mut_matrix_stranded_rep------------------------------------------------------------------
mut_mat_s_rep <- mut_matrix_stranded(grl, ref_genome, repli_strand_granges,
  mode = "replication"
)
strand_counts_rep <- strand_occurrences(mut_mat_s_rep, by = tissue)
strand_bias_rep <- strand_bias_test(strand_counts_rep)

## ----plot_strand_rep, fig.wide=TRUE-----------------------------------------------------------
ps1 <- plot_strand(strand_counts_rep, mode = "relative")
ps2 <- plot_strand_bias(strand_bias_rep)
grid.arrange(ps1, ps2)

## ----extract_signatures_bias------------------------------------------------------------------
nmf_res_strand <- extract_signatures(mut_mat_s, rank = 2, single_core = TRUE)
colnames(nmf_res_strand$signatures) <- c("Signature A", "Signature B")

## ----plot_rainfall, fig.wide=TRUE-------------------------------------------------------------
# Define autosomal chromosomes
chromosomes <- seqnames(get(ref_genome))[1:22]

# Make a rainfall plot
plot_rainfall(grl[[1]],
  title = names(grl[1]),
  chromosomes = chromosomes, cex = 1.5, ylim = 1e+09
)

## ----load_biomart-----------------------------------------------------------------------------
library(biomaRt)

## ----download_using_biomaRt-------------------------------------------------------------------
# regulatory <- useEnsembl(biomart="regulation",
#                          dataset="hsapiens_regulatory_feature",
#                          GRCh = 37)

## Download the regulatory CTCF binding sites and convert them to
## a GRanges object.
# CTCF <- getBM(attributes = c('chromosome_name',
#                             'chromosome_start',
#                             'chromosome_end',
#                             'feature_type_name'),
#              filters = "regulatory_feature_type_name",
#              values = "CTCF Binding Site",
#              mart = regulatory)
#
# CTCF_g <- reduce(GRanges(CTCF$chromosome_name,
#                 IRanges(CTCF$chromosome_start,
#                 CTCF$chromosome_end)))

CTCF_g <- readRDS(system.file("states/CTCF_g_data.rds",
  package = "MutationalPatterns"
))

## Download the promoter regions and convert them to a GRanges object.

# promoter = getBM(attributes = c('chromosome_name', 'chromosome_start',
#                                 'chromosome_end', 'feature_type_name'),
#                  filters = "regulatory_feature_type_name",
#                  values = "Promoter",
#                  mart = regulatory)
# promoter_g = reduce(GRanges(promoter$chromosome_name,
#                     IRanges(promoter$chromosome_start,
#                             promoter$chromosome_end)))

promoter_g <- readRDS(system.file("states/promoter_g_data.rds",
  package = "MutationalPatterns"
))

## Download the promoter flanking regions and convert them to a GRanges object.

# flanking = getBM(attributes = c('chromosome_name',
#                                 'chromosome_start',
#                                 'chromosome_end',
#                                 'feature_type_name'),
#                  filters = "regulatory_feature_type_name",
#                  values = "Promoter Flanking Region",
#                  mart = regulatory)
# flanking_g = reduce(GRanges(
#                        flanking$chromosome_name,
#                        IRanges(flanking$chromosome_start,
#                        flanking$chromosome_end)))

flanking_g <- readRDS(system.file("states/promoter_flanking_g_data.rds",
  package = "MutationalPatterns"
))

## ----combine_genomic_regions------------------------------------------------------------------
regions <- GRangesList(promoter_g, flanking_g, CTCF_g)

names(regions) <- c("Promoter", "Promoter flanking", "CTCF")

## ----combine_genomic_regions_2----------------------------------------------------------------
seqlevelsStyle(regions) <- "UCSC"

## ----download_bed_data------------------------------------------------------------------------
## Get the filename with surveyed/callable regions
surveyed_file <- system.file("extdata/callableloci-sample.bed",
  package = "MutationalPatterns"
)

## Import the file using rtracklayer and use the UCSC naming standard
library(rtracklayer)
surveyed <- import(surveyed_file)
seqlevelsStyle(surveyed) <- "UCSC"

## For this example we use the same surveyed file for each sample.
surveyed_list <- rep(list(surveyed), 9)

## ----genomic_distribution---------------------------------------------------------------------
distr <- genomic_distribution(grl, surveyed_list, regions)

## ----enrichment_depletion_test----------------------------------------------------------------
distr_test <- enrichment_depletion_test(distr, by = tissue)
head(distr_test)

## ----plot_enrichment_depletion, fig.wide=TRUE-------------------------------------------------
plot_enrichment_depletion(distr_test, sig_type = "p")

## ----split_muts_region------------------------------------------------------------------------
grl_region <- split_muts_region(grl, regions)
names(grl_region)

## ----extract_signatures_regions, fig.wide=TRUE------------------------------------------------
mut_mat_region <- mut_matrix(grl_region, ref_genome)
nmf_res_region <- extract_signatures(mut_mat_region, rank = 2, nrun = 10, single_core = TRUE)
nmf_res_region <- rename_nmf_signatures(nmf_res_region, 
                                        signatures, 
                                        cutoff = 0.85)
plot_contribution_heatmap(nmf_res_region$contribution, 
                          cluster_samples = TRUE, 
                          cluster_sigs = TRUE)

## ----plot_spectrum_region, fig.wide=TRUE------------------------------------------------------
type_occurrences_region <- mut_type_occurrences(grl_region, ref_genome)
plot_spectrum_region(type_occurrences_region)

## ----plot_spectrum_region_yaxis, fig.wide=TRUE------------------------------------------------
plot_spectrum_region(type_occurrences_region, mode = "relative_sample")

## ----make_long_profile------------------------------------------------------------------------
mut_mat_region <- mut_matrix(grl_region, ref_genome)
mut_mat_long <- lengthen_mut_matrix(mut_mat_region)
mut_mat_long[1:5, 1:5]

## ----plot_profile_region, fig.wide=TRUE-------------------------------------------------------
plot_profile_region(mut_mat_long[, c(1, 4, 7)])

## ----bin_mutation_density---------------------------------------------------------------------
regions_density <- bin_mutation_density(grl, ref_genome, nrbins = 3)
names(regions_density) <- c("Low", "Medium", "High")

## ----split_muts_region_density----------------------------------------------------------------
grl_region_dens <- split_muts_region(grl, regions_density, include_other = FALSE)

## ----combine substitutions--------------------------------------------------------------------
gr = unlist(grl)

## ----determine_regional_similarity------------------------------------------------------------
regional_sims <- determine_regional_similarity(gr,
                              ref_genome, 
                              chromosomes = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6"),
                              window_size = 40,
                              stepsize = 10,
                              max_window_size_gen = 40000000
)

## ----plot_regional_similarity-----------------------------------------------------------------
plot_regional_similarity(regional_sims)

## ----plot_lesion_segregation------------------------------------------------------------------
plot_lesion_segregation(grl[1:2])

## ----calculate_lesion_segregation-------------------------------------------------------------
lesion_segretation <- calculate_lesion_segregation(grl, sample_names)
head(lesion_segretation)

## ----calculate_lesion_segregation_split-------------------------------------------------------
lesion_seg_type <- calculate_lesion_segregation(grl,
                                                sample_names,
                                                split_by_type = TRUE,
                                                ref_genome = ref_genome)

## ----calculate_lesion_segregation_wald--------------------------------------------------------
lesion_segretation_wald <- calculate_lesion_segregation(grl, sample_names, 
                                                        test = "wald-wolfowitz")
head(lesion_segretation_wald)

## ----calculate_lesion_segregation_rl20--------------------------------------------------------
lesion_segretation_rl20 <- calculate_lesion_segregation(grl,
  sample_names,
  test = "rl20",
  ref_genome = ref_genome,
  chromosomes = chromosomes
)
head(lesion_segretation_rl20)

## ----plot_spectrum_xaxis----------------------------------------------------------------------
p <- plot_spectrum(type_occurrences, legend = FALSE)
p_axis <- p +
  theme(axis.text.y = element_text(size = 14, angle = 90))

## ----plot_spectrum_theme----------------------------------------------------------------------
p_theme <- p +
  theme_classic()

## ----combine_plot_spectrum2, eval=TRUE, fig.wide = TRUE, message=FALSE------------------------
grid.arrange(p, p_axis, p_theme, ncol = 3, widths = c(3, 3, 3))

## ----sessioninfo------------------------------------------------------------------------------
sessionInfo()


