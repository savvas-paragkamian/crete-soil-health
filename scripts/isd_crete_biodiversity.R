#!/usr/bin/Rscript

library(vegan)
library(tidyverse)

# Load data
abundance_asv <- readRDS("../dada2_output/all_runs_dada2_abundance_table.rds")
taxa_asv <- readRDS("../dada2_output/2023-07-27-dada2-taxa-silva-v138-1.RDS")
species_asv <- readRDS("../dada2_output/2023-07-31-dada2_taxa_species.RDS")
master_metadata <- read.delim("../dada2_output/Composite_MetaData_from_master.csv", sep=",")

# summary

# transform the matrices to long tables for stats

## this is for the species exact matching results from DADA2 addSpecies function
species <- data.frame(genus = species_asv[,6], species = species_asv[,8]) %>% 
    rownames_to_column("asv") %>%
    mutate(speciesname = paste(genus, species,sep=" "))# %>%
    group_by(speciesname) %>%
    summarise(n=n())

## the abundance matrix has all the biodiverstity information
## sampleID, ASV and abundance. NOTE that contains many zeros!!!
abundance_asv_l <- abundance_asv %>% 
    data.frame() %>%
    rownames_to_column("sampleID") %>%
    pivot_longer(!sampleID, names_to = "asv", values_to = "abundance")

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
    group_by(asv) %>% 
    ungroup()

taxa_asv_l_filter <- taxa_asv_l %>%
    group_by(asv) %>%
    summarise(higherClassification = paste0(higherClassification,
                                            collapse = "|"), 
              scientificName_all = paste0(scientificName,
                                       collapse = "|")) %>%
    ungroup()

# this function keeps the last occurrence of a string separated by |
keep_last <- function(x) tail(strsplit(x, split="\\|")[[1]],1)
# keep the lowest taxonomic inforfation per ASV
taxa_asv_l_filter$scientificName <- sapply(taxa_asv_l_filter$scientificName_all,
                                           keep_last,
                                           simplify = T, USE.NAMES=F)

taxa_asv_l_filter$classification <- sapply(taxa_asv_l_filter$higherClassification,
                                           keep_last,
                                           simplify = T, USE.NAMES=F)

# merge abundance and taxonomy

crete_biodiversity_all <- abundance_asv_l %>%
    left_join(taxa_asv_l_filter,by=c("asv"="asv")) %>%
    filter(abundance>0)

### remove the asvs that don't have a taxonomy
crete_biodiversity <- crete_biodiversity_all %>%
    filter(!is.na(classification)) 

write_delim(crete_biodiversity,"../results/crete_biodiversity_asv.tsv",delim="\t")
# Descriptives

## total ASVs and taxonomy
asv_no_taxonomy <- crete_biodiversity_all %>% 
    distinct(asv,classification) %>%
    filter(is.na(classification))

print(paste0("ASVs without taxonomy: ", nrow(asv_no_taxonomy)))

## total ASV summary
crete_biodiversity_asv <- crete_biodiversity %>% 
    group_by(asv) %>% 
    summarise(n_samples=n(),
              total_abundance=sum(abundance)) %>%
    ungroup()

asv_sample_dist <- crete_biodiversity_asv %>%
    group_by(n_samples) %>%
    summarise(n_asv=n(), sum_abundance=sum())

write_delim(asv_sample_dist,"../results/asv_sample_dist.tsv",delim="\t")

## taxonomic diversity per sample

crete_biodiversity_s <- crete_biodiversity %>% 
    filter(abundance>10, !is.na(classification)) %>%
    group_by(sampleID, classification) %>% 
    summarise(n_taxa=n(), .groups="keep")
## create bar plots for each sample at family level, class level, Phylum etc

## shannon per sample
## bray curtis per sample
## 
