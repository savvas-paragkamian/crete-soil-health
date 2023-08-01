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
abundance_asv_l <- abundance_asv %>% 
    data.frame() %>%
    rownames_to_column("sampleID") %>%
    pivot_longer(!sampleID, names_to = "asv", values_to = "abundance")

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
#    mutate(scientificName = str_split(scientificName_all, "|", n=1, simplify=T)) %>%
#    mutate(scientificName=tail(strsplit(scientificName_all, split="\\|")[[1]],1)) %>%
#    mutate(scientificName = sub('^.*\\|([[:alnum:]]+)$', '\\1', scientificName_all)) %>%
    ungroup()

#taxa_asv_l_filter$scientificName <- sapply(tail(strsplit(taxa_asv_l_filter$scientificName_all, split="\\|")[[1]],1))

species <- data.frame(genus = species_asv[,6], species = species_asv[,8]) %>% 
    rownames_to_column("asv") %>%
    mutate(speciesname = paste(genus, species,sep=" "))# %>%
    group_by(speciesname) %>%
    summarise(n=n())

##
table(taxa_asv_l$Classification_levels)
table(taxa_asv_l$Classification_levels, taxa_asv_l$higherClassification)
