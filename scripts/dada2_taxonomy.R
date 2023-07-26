#!/usr/bin/env Rscript

###############################################################################
# script name: dada2_taxonomy.R
# developed by: Johanna Holms, Savvas Paragkamian
# framework: ISD Crete
###############################################################################
# GOAL:
# Aim of this script is to assign taxonomy to ASVs with DADA2 of ISD Crete 2016
###############################################################################
# usage:./dada2_taxonomy.R
###############################################################################

library(dada2,lib.loc="/home1/s.paragkamian/software/R/4.1.1")

#faf
seqtab.nochim<-readRDS("/home1/s.paragkamian/isd-crete/Crete/all_runs_dada2_abundance_table.rds")

#seqtab.nochim<- seqtab.nochim[,1:1000]

print("test data initiated.")
taxa <- assignTaxonomy(seqtab.nochim,
                       "/home1/s.paragkamian/databases/SILVA_138_SSU/silva_nr99_v138.1_wSpecies_train_set.fa",
                       multithread=20,
                       tryRC = TRUE, 
                       verbose = T)

print("assignTaxonomy done.")
taxa <- addSpecies(taxa,
                   "/home1/s.paragkamian/databases/SILVA_138_SSU/silva_species_assignment_v138.1.fa",
                   tryRC = TRUE)

print("addSpeciesdone.")

print("begin saving data.")
saveRDS(taxa, "/home1/s.paragkamian/isd-crete/output/dada2_taxa.RDS")
print("data saved. \n")
