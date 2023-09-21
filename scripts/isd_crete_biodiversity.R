#!/usr/bin/env Rscript

###############################################################################
# script name: isd_crete_biodiversity.R
# developed by: Savvas Paragkamian, Johanna Holms
# framework: ISD Crete
###############################################################################
# GOAL:
# Aim of this script is to use ASVs and the taxonomy to filter and normalise
# the sample biodiversity
#
###############################################################################
# OUTPUT:
# the main output of this script is the crete_biodiversity_asv.tsv
# which contains the sample - asv occurrences with abundance along with 
# taxonomic information.
#
# The minimum information of this file is 
# ENA-RUN, asv, abundance
# Taxonomy is also included.
# 
# Other 5 files are produced
# crete_biodiversity_matrix.RDS, a matrix of abundances
# tax_tab.RDS, taxonomy table with the remaining asvs
# sample_stats.tsv
# sample_stats_total.tsv
# asv_stas.tsv
###############################################################################
# usage:./isd_crete_biodiversity.R
###############################################################################
library(dplyr)
library(tibble)
library(readr)
library(magrittr)
library(tidyr)
# Load data old repo
abundance_asv_old <- readRDS("Crete/all_runs_dada2_abundance_table.rds")
master_metadata_old <- read.delim("Crete/Composite_MetaData_from_master.csv", sep=",")

# Load data
abundance_asv_long <- read_delim("dada2_output/taxonomy/seqtab_nochim_long.tsv", delim="\t") %>% 
    mutate(ENA_RUN=gsub("_1_filt.fastq.gz", "", file, perl=T))
#abundance_asv <- readRDS("dada2_output/taxonomy/seqtab_nochim.RDS") # use this one if data are RDS
asv_fasta <- read_delim("dada2_output/taxonomy/asv_fasta_ids.tsv", delim="\t")
taxa_asv <- readRDS("dada2_output/taxonomy/dada2_taxonomy.RDS")
species_asv <- readRDS("dada2_output/taxonomy/dada2_taxa_species.RDS")

# transform the matrices to long tables for stats

## this is for the species exact matching results from DADA2 addSpecies function
species <- data.frame(genus = species_asv[,6], species = species_asv[,8]) %>% 
    rownames_to_column("asv") %>%
    mutate(speciesname = paste(genus, species,sep=" ")) %>%
    group_by(speciesname) %>%
    summarise(n=n())

## the abundance matrix has all the biodiverstity information
## sampleID, ASV and abundance. NOTE that contains many zeros!!!
## UNCOMMENT if data are from RDS!
#abundance_asv_l <- abundance_asv %>% 
#    data.frame() %>%
#    rownames_to_column("sampleID") %>%
#    pivot_longer(!sampleID, names_to = "asv", values_to = "abundance")

## the taxonomy matrix has all the taxonomic information of each ASV
taxa_asv_l <- taxa_asv %>% 
    data.frame() %>%
    rownames_to_column("asv") %>%
    mutate(Species = ifelse(is.na(Species),
                                 NA,
                                 paste(Genus, Species,sep=" "))) %>%
    pivot_longer(!asv,
                 names_to = "higherClassification",
                 values_to = "scientificName") %>%
    na.omit("scientificName") %>%
    ungroup() %>%
    left_join(asv_fasta, by=c("asv"="asv"))


taxa_asv_s <- taxa_asv_l %>%
    group_by(asv) %>%
    summarise(higherClassification = paste0(higherClassification,
                                            collapse = "|"), 
              scientificName_all = paste0(scientificName,
                                       collapse = "|")) %>%
    ungroup()

# this function keeps the last occurrence of a string separated by |
keep_last <- function(x) tail(strsplit(x, split="\\|")[[1]],1)
# keep the lowest taxonomic inforfation per ASV
taxa_asv_s$scientificName <- sapply(taxa_asv_s$scientificName_all,
                                           keep_last,
                                           simplify = T, USE.NAMES=F)

taxa_asv_s$classification <- sapply(taxa_asv_s$higherClassification,
                                           keep_last,
                                           simplify = T, USE.NAMES=F)

taxa_asv_all <- taxa_asv_l %>% 
    pivot_wider(names_from=higherClassification, 
                values_from=scientificName) %>% 
    left_join(taxa_asv_s, by=c("asv"="asv"))

# merge abundance and taxonomy

crete_biodiversity_all <- abundance_asv_long %>%
    left_join(taxa_asv_all,by=c("asv_id"="asv_id")) %>%
    filter(abundance>0)

# Decontamination
## very important step, we need the samples controls.

# This is the MASTER dataset of Crete biodiversity and taxonomy
### remove the asvs that don't have a taxonomy

crete_biodiversity <- crete_biodiversity_all %>%
    filter(!is.na(classification), abundance < 100)

write_delim(crete_biodiversity,"results/crete_biodiversity_asv.tsv",delim="\t")

## create a abundance matrix
crete_biodiversity_m <- crete_biodiversity %>%
    select(ENA_RUN, asv_id, abundance) %>%
    pivot_wider(names_from=ENA_RUN, values_from=abundance, values_fill = 0) %>%
    as.matrix()

crete_biodiversity_matrix <- crete_biodiversity_m[,-1]
crete_biodiversity_matrix <- apply(crete_biodiversity_matrix, 2, as.numeric)
rownames(crete_biodiversity_matrix) <- crete_biodiversity_m[,1]
saveRDS(crete_biodiversity_matrix, "results/crete_biodiversity_matrix.RDS")

# Create taxonomy table of the remaining asvs
tax_tab1 <- crete_biodiversity %>%
    distinct(asv_id, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
    as.matrix()

tax_tab <- tax_tab1[,-1]
rownames(tax_tab) <- tax_tab1[,1]
saveRDS(tax_tab, "results/tax_tab.RDS")

## total ASVs and taxonomy
asv_no_taxonomy <- crete_biodiversity_all %>% 
    distinct(asv,classification) %>%
    filter(is.na(classification))

print(paste0("ASVs without taxonomy: ", nrow(asv_no_taxonomy)))

## total ASV summary
asv_stats <- crete_biodiversity %>% 
    group_by(asv_id, classification, scientificName) %>% 
    summarise(n_samples=n(),
              reads=sum(abundance), .groups="keep") %>%
    ungroup()

write_delim(asv_stats,"results/asv_stats.tsv",delim="\t")
## singletons
singletons <- asv_stats %>% 
    filter(reads==1) %>%
    nrow()

print(paste0("there are ",singletons," singletons asvs"))

## asv and samples distribution
asv_sample_dist <- asv_stats %>%
    group_by(n_samples) %>%
    summarise(n_asv=n(), reads=sum(reads))

## taxonomic, asv and read diversity per sample
sample_reads <- crete_biodiversity %>%
    group_by(ENA_RUN) %>%
    summarise(asvs=n(),reads=sum(abundance))

summary(sample_reads)

sample_stats <- crete_biodiversity %>% 
    group_by(ENA_RUN, classification, scientificName) %>% 
    summarise(asvs=n(), reads=sum(abundance), .groups="keep") %>%
    group_by(ENA_RUN,classification) %>%
    summarise(taxa=n(),reads=sum(reads), asvs=sum(asvs), .groups="keep") %>%
#    pivot_wider(names_from=classification,values_from=n_taxa) %>%
    ungroup()

write_delim(sample_stats,
            "results/sample_stats.tsv",
            delim="\t")

sample_stats_total <- sample_stats %>%
    group_by(ENA_RUN) %>%
    summarise(taxa=sum(taxa),reads=sum(reads), asvs=sum(asvs))

write_delim(sample_stats_total,
            "results/sample_stats_total.tsv",
            delim="\t")

