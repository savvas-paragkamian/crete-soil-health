#!/usr/bin/env Rscript
###############################################################################
# script name: isd_crete_compositionality.R
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
# usage:./scripts/isd_crete_compositionality.R
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

#output_label2 = ancombc2(data = tse, assay_name = "counts", tax_level = NULL,
#                  fix_formula = "LABEL2 + elevation_bin", rand_formula = NULL,
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
#saveRDS(output_label2, "results/ancombc2_results.RDS")

### geology
print("Starting ancombc2 for geology")

output_geology = ancombc2(data = tse, assay_name = "counts", tax_level = NULL,
                  fix_formula = "geology_na + total_nitrogen", rand_formula = NULL,
                  p_adj_method = "holm", pseudo_sens = TRUE,
                  prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                  group = "geology_na", struc_zero = FALSE, neg_lb = FALSE,
                  alpha = 0.05, n_cl = 2, verbose = TRUE,
                  global = FALSE, pairwise = TRUE, 
                  dunnet = FALSE, trend = FALSE,
                  iter_control = list(tol = 1e-5, max_iter = 20, 
                                      verbose = FALSE),
                  em_control = list(tol = 1e-5, max_iter = 100),
                  lme_control = NULL, 
                  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100), 
                  trend_control = NULL)

saveRDS(output_geology, "results/ancombc2_geology_results.RDS")
