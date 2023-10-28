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
# OUTPUT: NOT WORKING!!!
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

crete_biodiversity_matrix <- readRDS("results/crete_biodiversity_matrix.RDS")
biodiversity_srs <- readRDS("results/biodiversity_srs.RDS")
tax_tab <- readRDS("results/tax_tab.RDS")
metadata <- read_delim("results/sample_metadata.tsv", delim="\t")
# ancombc analysis

assays = SimpleList(counts = crete_biodiversity_matrix)
smd = DataFrame(metadata)
tax_tab = DataFrame(tax_tab)

## create TSE
tse = TreeSummarizedExperiment(assays = assays,
                               colData = smd,
                               rowData = tax_tab)
## convert TSE to phyloseq
pseq = makePhyloseqFromTreeSummarizedExperiment(tse)

total = median(sample_sums(pseq))
standf = function(x, t=total) round(t * (x / sum(x)))

pseq_std  = transform_sample_counts(pseq,standf )
gpsf = filter_taxa(pseq_std, function(x) sd(x)/mean(x) > 3.0, TRUE)

gpsfb = subset_taxa(gpsf, Order=="Frankiales")

jpeg(file="saving_plot1.jpeg")
plot_bar(gpsfb, "LABEL1", "Abundance", "Family", title="test")
dev.off()


GP = pseq_std
#wh0 = genefilter_sample(GP, filterfun_sample(function(x) x > 5), A=0.5*nsamples(GP))
#GP1 = prune_taxa(wh0, GP)

#GP1 = transform_sample_counts(GP, function(x) 1E6 * x/sum(x))

GP.ord <- ordinate(pseq, "NMDS", "bray")

jpeg(file="saving_plot2.jpeg")
p1 = plot_ordination(GP, GP.ord, type="samples", color="LABEL1", title="taxa")
p1
dev.off()


###
 otu_mat = matrix(sample(1:100, 100, replace = TRUE), nrow = 10, ncol = 10)
     rownames(otu_mat) = paste0("taxon", 1:nrow(otu_mat))
     colnames(otu_mat) = paste0("sample", 1:ncol(otu_mat))
     assays = SimpleList(counts = otu_mat)

     # sample metadata
     smd = data.frame(group = sample(LETTERS[1:4], size = 10, replace = TRUE),
                      row.names = paste0("sample", 1:ncol(otu_mat)),
                      stringsAsFactors = FALSE)
     smd = DataFrame(smd)

     # taxonomy table
     tax_tab = matrix(sample(letters, 70, replace = TRUE),
                      nrow = nrow(otu_mat), ncol = 7)             
     rownames(tax_tab) = rownames(otu_mat)                                                                                                                                                                            
     colnames(tax_tab) = c("Kingdom", "Phylum", "Class", "Order",                                                                                                                                                     
                           "Family", "Genus", "Species")                                                                                                                                                              
     tax_tab = DataFrame(tax_tab)                                                                                                                                                                                     
                                                                                                                                                                                                                      
     # create TSE                                                                                                                                                                                                     
     tse = TreeSummarizedExperiment(assays = assays,                                                                                                                                                                  
                                    colData = smd,                                                                                                                                                                    
                                    rowData = tax_tab)   

pseq = makePhyloseqFromTreeSummarizedExperiment(tse)

data(dietswap, package = "microbiome")                                                                                                                                                                           
     tse = mia::makeTreeSummarizedExperimentFromPhyloseq(dietswap)                                                                                                                                                    
                                                                                                                                                                                                                      
     colData(tse)$bmi_group = factor(colData(tse)$bmi_group,                                                                                                                                                          
                                     levels = c("obese",                                                                                                                                                              
                                                "overweight",                                                                                                                                                         
                                                "lean"))                                                                                                                                                              
     # Note that setting pseudo_sens = FALSE, max_iter = 1, and B = 1 is                                                                                                                                              
     # only for the sake of speed                                                                                                                                                                                     
     # Set pseudo_sens = TRUE, and use default or larger values for max_iter and B                                                                                                                                    
     # for better performance                   
     #
     #
     #
     #
  out = ancombc2(data = tse, assay_name = "counts", tax_level = "Phylum",
                 fix_formula = "nationality + timepoint + bmi_group",
                 rand_formula = "(timepoint | subject)",
                 p_adj_method = "holm", pseudo = 0, pseudo_sens = FALSE,
                 prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                 group = "bmi_group", struc_zero = TRUE, neg_lb = TRUE,
                 alpha = 0.05, n_cl = 1, verbose = TRUE,
                 global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = TRUE,
                 iter_control = list(tol = 1e-2, max_iter = 1, verbose = TRUE),
                 em_control = list(tol = 1e-5, max_iter = 1),
                 lme_control = lme4::lmerControl(),
                 mdfdr_control = list(fwer_ctrl_method = "holm", B = 1),
                 trend_control = list(contrast =list(matrix(c(1, 0, -1, 1),nrow = 2,byrow = TRUE)), 
                                      node = list(2),
                                      solver = "ECOS",B = 1))

res_prim = out$res

res_global = out$res_global

res_pair = out$res_pair

res_dunn = out$res_dunn

res_trend = out$res_trend

set.seed(123)
data("GlobalPatterns", package = "mia")  
se <- GlobalPatterns

ancombc2_out <- ancombc2(data = se,
                         tax_level = "Family",
                         group = "pairwise")


output = ancombc2(data = tse, assay_name = "counts", tax_level = "Family",
                  fix_formula = NULL, rand_formula = NULL,
                  p_adj_method = "holm", 
                  prv_cut = 0.10, lib_cut = 100, s0_perc = 0.05,
                  group = NULL, struc_zero = FALSE, neg_lb = FALSE,
                  alpha = 0.05, n_cl = 2, verbose = TRUE,
                  global = FALSE, pairwise = FALSE, 
                  dunnet = FALSE, trend = FALSE,
                  iter_control = list(tol = 1e-5, max_iter = 20, 
                                      verbose = FALSE),
                  em_control = list(tol = 1e-5, max_iter = 100),
                  lme_control = NULL, 
                  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100), 
                  trend_control = NULL)


output = ancombc2(data = tse, assay_name = "counts", tax_level = "Family",
                  fix_formula = NULL, rand_formula = NULL,
                  p_adj_method = "holm", 
                  prv_cut = 0.10, lib_cut = 100, s0_perc = 0.05,
                  group = NULL, struc_zero = FALSE, neg_lb = FALSE,
                  alpha = 0.05, n_cl = 2, verbose = TRUE,
                  global = FALSE, pairwise = FALSE, 
                  dunnet = FALSE, trend = FALSE,
                  iter_control = list(tol = 1e-5, max_iter = 20, 
                                      verbose = FALSE),
                  em_control = list(tol = 1e-5, max_iter = 100),
                  lme_control = NULL, 
                  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100), 
                  trend_control = NULL)

data(dietswap, package = "microbiome")
    tse = mia::makeTreeSummarizedExperimentFromPhyloseq(dietswap)
    colData(tse)$bmi_group = factor(colData(tse)$bmi_group,
                                    levels = c("obese",
                                               "overweight",
                                               "lean"))
set.seed(123)
# Note that setting max_iter = 1 and B = 1 is only for the sake of speed
# Use default or larger values for max_iter and B for better performance
out = ancombc2(data = tse, assay_name = "counts", tax_level = "Phylum",
               fix_formula = "nationality + timepoint + bmi_group",
               rand_formula = "(timepoint | subject)",
               p_adj_method = "holm",
               prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
               group = "bmi_group", struc_zero = TRUE, neg_lb = TRUE,
               alpha = 0.05, n_cl = 1, verbose = TRUE, global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = TRUE,
               iter_control = list(tol = 1e-2, max_iter = 1, verbose = TRUE),
               em_control = list(tol = 1e-5, max_iter = 1),
               lme_control = lme4::lmerControl(),
               mdfdr_control = list(fwer_ctrl_method = "holm", B = 1),
               trend_control = list(contrast =list(matrix(c(1, 0, -1, 1),nrow = 2, byrow = TRUE)),node = list(2), solver = "ECOS", B = 1))
res_prim = out$res
res_global = out$res_global
res_pair = out$res_pair
res_dunn = out$res_dunn
res_trend = out$res_trend






