#!/usr/bin/env Rscript
###############################################################################
# script name: isd_crete_differential_abundance.R
# developed by: Savvas Paragkamian
# framework: ISD Crete
###############################################################################
# GOAL:
# Aim of this script is to normalise the reads of each sample biodiversity using
# the ANCOMBC methodology
#
###############################################################################
# 75 in minutes
# OUTPUT: RDS files of each run
###############################################################################
# usage:./scripts/isd_crete_differential_abundance.R
###############################################################################

library(magrittr)
library(dplyr)
library(tibble)
library(readr)
library(tidyr)
library(mia)
library(ANCOMBC)
library(phyloseq)

community_matrix_l <- read_delim("results/community_matrix_l.tsv",delim="\t") 
community_matrix <- readRDS("results/community_matrix.RDS") |> as.matrix() |> t()
tax_tab <- readRDS("results/tax_tab.RDS")
results <- readRDS("results/ancombc2_results.RDS")
metadata <- read_delim("results/sample_metadata.tsv", delim="\t")

tax_tab1 <- community_matrix_l %>%
    distinct(scientificName, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
    as.matrix()

tax_tab <- tax_tab1[,-1]
rownames(tax_tab) <- tax_tab1[,1]

metadata <- metadata %>% filter(ENA_RUN %in% colnames(community_matrix))
# ancombc analysis

assays = S4Vectors::SimpleList(counts = community_matrix)
smd = S4Vectors::DataFrame(metadata)
tax_tab = S4Vectors::DataFrame(tax_tab)

#otu_mat = matrix(sample(1:100, 100, replace = TRUE), nrow = 10, ncol = 10)
#rownames(otu_mat) = paste0("taxon", 1:nrow(otu_mat))
#colnames(otu_mat) = paste0("sample", 1:ncol(otu_mat))
#assays = SimpleList(counts = otu_mat)
## create TSE
tse = TreeSummarizedExperiment::TreeSummarizedExperiment(assays = assays, colData = smd)
set.seed(123)

print("Starting ancombc2 for Label2")

output_label2 = ancombc2(data = tse, assay_name = "counts", tax_level = NULL,
                  fix_formula = "LABEL2", rand_formula = NULL,
                  p_adj_method = "holm", pseudo_sens = TRUE,
                  prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                  group = "LABEL2", struc_zero = FALSE, neg_lb = FALSE,
                  alpha = 0.05, n_cl = 2, verbose = TRUE,
                  global = FALSE, pairwise = TRUE, 
                  dunnet = FALSE, trend = FALSE,
                  iter_control = list(tol = 1e-5, max_iter = 20, 
                                      verbose = FALSE),
                  em_control = list(tol = 1e-5, max_iter = 100),
                  lme_control = NULL, 
                  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100), 
                  trend_control = NULL)

saveRDS(output_label2, "results/ancombc2_label2_results.RDS")
#print("Starting ancombc2 for bioclim")
#output_clim = ancombc2(data = tse, assay_name = "counts", tax_level = NULL,
#                  fix_formula = "elevation_bin + bio_1 + bio_12 + elevation", rand_formula = NULL,
#                  p_adj_method = "holm", pseudo_sens = TRUE,
#                  prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
#                  group = "elevation_bin", struc_zero = FALSE, neg_lb = FALSE,
#                  alpha = 0.05, n_cl = 2, verbose = TRUE,
#                  global = FALSE, pairwise = TRUE, 
#                  dunnet = FALSE, trend = FALSE,
#                  iter_control = list(tol = 1e-5, max_iter = 20, 
#                                      verbose = FALSE),
#                  em_control = list(tol = 1e-5, max_iter = 100),
#                  lme_control = NULL, 
#                  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100), 
#                  trend_control = NULL)
#
#saveRDS(output_clim, "results/ancombc2_results_clim.RDS")
#
#print("Starting ancombc2 for metadata")
#output_metadata = ancombc2(data = tse, assay_name = "counts", tax_level = NULL,
#                  fix_formula = "elevation_bin + total_nitrogen + water_content + total_organic_carbon + carbon_nitrogen_ratio", rand_formula = NULL,
#                  p_adj_method = "holm", pseudo_sens = TRUE,
#                  prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
#                  group = "elevation_bin", struc_zero = FALSE, neg_lb = FALSE,
#                  alpha = 0.05, n_cl = 2, verbose = TRUE,
#                  global = FALSE, pairwise = TRUE, 
#                  dunnet = FALSE, trend = FALSE,
#                  iter_control = list(tol = 1e-5, max_iter = 20, 
#                                      verbose = FALSE),
#                  em_control = list(tol = 1e-5, max_iter = 100),
#                  lme_control = NULL, 
#                  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100), 
#                  trend_control = NULL)
#
#saveRDS(output_metadata, "results/ancombc2_results_metadata.RDS")
#
### geology
#print("Starting ancombc2 for geology")
#
#output_geology = ancombc2(data = tse, assay_name = "counts", tax_level = NULL,
#                  fix_formula = "geology_na + total_nitrogen", rand_formula = NULL,
#                  p_adj_method = "holm", pseudo_sens = TRUE,
#                  prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
#                  group = "geology_na", struc_zero = FALSE, neg_lb = FALSE,
#                  alpha = 0.05, n_cl = 2, verbose = TRUE,
#                  global = FALSE, pairwise = TRUE, 
#                  dunnet = FALSE, trend = FALSE,
#                  iter_control = list(tol = 1e-5, max_iter = 20, 
#                                      verbose = FALSE),
#                  em_control = list(tol = 1e-5, max_iter = 100),
#                  lme_control = NULL, 
#                  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100), 
#                  trend_control = NULL)
#
#saveRDS(output_geology, "results/ancombc2_geology_results.RDS")
#
############################## results ##################################
#
#output_label2 <- readRDS( "results/ancombc2_label2_results.RDS")
#
significant_label2 <- filter(output_label2$res,if_any(where(is.logical),~ . =="TRUE"))
write_delim(output_label2$res, "results/ancombc2_label2_results_only.tsv", delim="\t")
write_delim(significant_label2, "results/ancombc2_significant_label2.tsv", delim="\t")
#
#output_geology<- readRDS( "results/ancombc2_geology_results.RDS")
#significant_geology <- filter(output_geology$res,if_any(where(is.logical),~ . =="TRUE"))
#write_delim(output_geology$res, "results/ancombc2_geology_results_only.tsv", delim="\t")
#
#write_delim(significant_geology, "results/ancombc2_significant_geology.tsv", delim="\t")

#output_metadata <- readRDS("results/ancombc2_results_metadata.RDS")
#significant_metadata <- filter(output_metadata$res,if_any(where(is.logical),~ . =="TRUE"))
#write_delim(significant_metadata, "results/ancombc2_significant_metadata.tsv", delim="\t")

#output_clim <- readRDS("results/ancombc2_results_clim.RDS")
#significant_clim <- filter(output_clim$res,if_any(where(is.logical),~ . =="TRUE"))
#write_delim(significant_clim, "results/ancombc2_significant_clim.tsv", delim="\t")


#output_label2 = ancombc2(data = tse, assay_name = "counts", tax_level = NULL,
#                  fix_formula = "LABEL2", rand_formula = NULL,
#                  p_adj_method = "holm", pseudo_sens = TRUE,
#                  prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
#                  group = "LABEL2", struc_zero = FALSE, neg_lb = FALSE,
#                  alpha = 0.05, n_cl = 2, verbose = TRUE,
#                  global = FALSE, pairwise = TRUE, 
#                  dunnet = FALSE, trend = FALSE,
#                  iter_control = list(tol = 1e-5, max_iter = 20, 
#                                      verbose = FALSE),
#                  em_control = list(tol = 1e-5, max_iter = 100),
#                  lme_control = NULL, 
#                  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100), 
#                  trend_control = NULL)
#
#saveRDS(output_label2, "results/ancombc2_label2_results.RDS")

#output_label3 = ancombc2(data = tse, assay_name = "counts", tax_level = NULL,
#                  fix_formula = "LABEL3", rand_formula = NULL,
#                  p_adj_method = "holm", pseudo_sens = TRUE,
#                  prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
#                  group = "LABEL3", struc_zero = FALSE, neg_lb = FALSE,
#                  alpha = 0.05, n_cl = 2, verbose = TRUE,
#                  global = FALSE, pairwise = TRUE, 
#                  dunnet = FALSE, trend = FALSE,
#                  iter_control = list(tol = 1e-5, max_iter = 20, 
#                                      verbose = FALSE),
#                  em_control = list(tol = 1e-5, max_iter = 100),
#                  lme_control = NULL, 
#                  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100), 
#                  trend_control = NULL)
#
#saveRDS(output_label3, "results/ancombc2_label3_results.RDS")
#print("Starting ancombc2 for bioclim")
#output_clim = ancombc2(data = tse, assay_name = "counts", tax_level = NULL,
#                  fix_formula = "elevation_bin + bio_1 + bio_12 + elevation", rand_formula = NULL,
#                  p_adj_method = "holm", pseudo_sens = TRUE,
#                  prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
#                  group = "elevation_bin", struc_zero = FALSE, neg_lb = FALSE,
#                  alpha = 0.05, n_cl = 2, verbose = TRUE,
#                  global = FALSE, pairwise = TRUE, 
#                  dunnet = FALSE, trend = FALSE,
#                  iter_control = list(tol = 1e-5, max_iter = 20, 
#                                      verbose = FALSE),
#                  em_control = list(tol = 1e-5, max_iter = 100),
#                  lme_control = NULL, 
#                  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100), 
#                  trend_control = NULL)
#
#saveRDS(output_clim, "results/ancombc2_results_clim.RDS")
#
#print("Starting ancombc2 for metadata")
#output_metadata = ancombc2(data = tse, assay_name = "counts", tax_level = NULL,
#                  fix_formula = "elevation_bin + total_nitrogen + water_content + total_organic_carbon + carbon_nitrogen_ratio", rand_formula = NULL,
#                  p_adj_method = "holm", pseudo_sens = TRUE,
#                  prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
#                  group = "elevation_bin", struc_zero = FALSE, neg_lb = FALSE,
#                  alpha = 0.05, n_cl = 2, verbose = TRUE,
#                  global = FALSE, pairwise = TRUE, 
#                  dunnet = FALSE, trend = FALSE,
#                  iter_control = list(tol = 1e-5, max_iter = 20, 
#                                      verbose = FALSE),
#                  em_control = list(tol = 1e-5, max_iter = 100),
#                  lme_control = NULL, 
#                  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100), 
#                  trend_control = NULL)
#
#saveRDS(output_metadata, "results/ancombc2_results_metadata.RDS")
#
### geology
#print("Starting ancombc2 for geology")
#
#output_geology = ancombc2(data = tse, assay_name = "counts", tax_level = NULL,
#                  fix_formula = "geology_na + total_nitrogen", rand_formula = NULL,
#                  p_adj_method = "holm", pseudo_sens = TRUE,
#                  prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
#                  group = "geology_na", struc_zero = FALSE, neg_lb = FALSE,
#                  alpha = 0.05, n_cl = 2, verbose = TRUE,
#                  global = FALSE, pairwise = TRUE, 
#                  dunnet = FALSE, trend = FALSE,
#                  iter_control = list(tol = 1e-5, max_iter = 20, 
#                                      verbose = FALSE),
#                  em_control = list(tol = 1e-5, max_iter = 100),
#                  lme_control = NULL, 
#                  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100), 
#                  trend_control = NULL)
#
#saveRDS(output_geology, "results/ancombc2_geology_results.RDS")
#
# Aridity
print("aridity and desertification")


output_esa = ancombc2(data = tse, assay_name = "counts", tax_level = NULL,
                  fix_formula = "ESA_12CL", rand_formula = NULL,
                  p_adj_method = "holm", pseudo_sens = TRUE,
                  prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                  group = "ESA_12CL", struc_zero = FALSE, neg_lb = FALSE,
                  alpha = 0.05, n_cl = 2, verbose = TRUE,
                  global = FALSE, pairwise = TRUE, 
                  dunnet = FALSE, trend = FALSE,
                  iter_control = list(tol = 1e-5, max_iter = 20, 
                                      verbose = FALSE),
                  em_control = list(tol = 1e-5, max_iter = 100),
                  lme_control = NULL, 
                  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100), 
                  trend_control = NULL)

saveRDS(output_esa, "results/ancombc2_esa_results.RDS")

output_arid = ancombc2(data = tse, assay_name = "counts", tax_level = NULL,
                  fix_formula = "aridity_class", rand_formula = NULL,
                  p_adj_method = "holm", pseudo_sens = TRUE,
                  prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                  group = "aridity_class", struc_zero = FALSE, neg_lb = FALSE,
                  alpha = 0.05, n_cl = 2, verbose = TRUE,
                  global = FALSE, pairwise = TRUE, 
                  dunnet = FALSE, trend = FALSE,
                  iter_control = list(tol = 1e-5, max_iter = 20, 
                                      verbose = FALSE),
                  em_control = list(tol = 1e-5, max_iter = 100),
                  lme_control = NULL, 
                  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100), 
                  trend_control = NULL)

saveRDS(output_arid, "results/ancombc2_aridity_class_results.RDS")
############################## results ##################################
#
#output_label2 <- readRDS( "results/ancombc2_label2_results.RDS")
#
#significant_label2 <- filter(output_label2$res,if_any(where(is.logical),~ . =="TRUE"))
#write_delim(output_label2$res, "results/ancombc2_label2_results_only.tsv", delim="\t")
#write_delim(significant_label2, "results/ancombc2_significant_label2.tsv", delim="\t")
#significant_label3 <- filter(output_label3$res,if_any(where(is.logical),~ . =="TRUE"))
#write_delim(output_label3$res, "results/ancombc2_label3_results_only.tsv", delim="\t")
#write_delim(significant_label3, "results/ancombc2_significant_label3.tsv", delim="\t")
#
#output_geology<- readRDS( "results/ancombc2_geology_results.RDS")
#significant_geology <- filter(output_geology$res,if_any(where(is.logical),~ . =="TRUE"))
#write_delim(output_geology$res, "results/ancombc2_geology_results_only.tsv", delim="\t")
#
#write_delim(significant_geology, "results/ancombc2_significant_geology.tsv", delim="\t")

#output_metadata <- readRDS("results/ancombc2_results_metadata.RDS")
#significant_metadata <- filter(output_metadata$res,if_any(where(is.logical),~ . =="TRUE"))
#write_delim(significant_metadata, "results/ancombc2_significant_metadata.tsv", delim="\t")

#output_clim <- readRDS("results/ancombc2_results_clim.RDS")
#significant_clim <- filter(output_clim$res,if_any(where(is.logical),~ . =="TRUE"))
#write_delim(significant_clim, "results/ancombc2_significant_clim.tsv", delim="\t")
