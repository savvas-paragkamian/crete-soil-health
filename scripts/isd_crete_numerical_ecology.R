#!/usr/bin/Rscript

###############################################################################
# script name: isd_crete_numerical_ecology.R
# developed by: Savvas Paragkamian, Johanna Holms
# framework: ISD Crete
###############################################################################
# GOAL:
# Aim of this script is to use Crete biodiversity data of ASVs, taxonomy and 
# the sample metadata to perform ecological analyses on biodiversity, ordination
# and multivariate comparison.
#
###############################################################################
# OUTPUT:
#
###############################################################################
# usage:./isd_crete_numerical_ecology.R
###############################################################################
library(vegan)
library(dplyr)
library(tibble)
library(readr)
library(magrittr)
library(tidyr)

# Load data
master_metadata_old <- read.delim("Crete/Composite_MetaData_from_master.csv", sep=",")

metadata_long <- read_delim("ena_metadata/ena_isd_2016_attributes.tsv", delim="\t") %>%
    mutate(VALUE=gsub("\\r(?!\\n)","", VALUE, perl=T))

crete_biodiversity <- read_delim("results/crete_biodiversity_asv.tsv",delim="\t")
#
#
#
#
# Metadata

#summary
# metadata to wide format

metadata_wide <- metadata_long %>% 
    dplyr::select(-c(UNITS)) %>%
    mutate(TAG=gsub(" ","_", TAG, perl=T)) %>%
    pivot_wider(names_from=TAG, 
                values_from=VALUE)

metadata_wide$total_nitrogen <- as.numeric(metadata_wide$total_nitrogen)
metadata_wide$water_content <- as.numeric(metadata_wide$water_content)
metadata_wide$total_organic_carbon <- as.numeric(metadata_wide$total_organic_carbon)
# differences of old and new data

master_metadata_old$team_site_location_id[which(!(master_metadata_old$team_site_location_id %in% metadata_wide$source_material_identifiers))]



### highest species biodiversity sample
print("sample with the highest microbial species diversity")
crete_biodiversity_s[which(crete_biodiversity_s$Species==max(crete_biodiversity_s$Species)),]

highest_crete_biodiversity_s <- crete_biodiversity_s %>% 
    left_join(metadata_wide,
              by=c("ENA-RUN"="ENA-RUN")) %>% 
    select(place_name,source_material_identifiers, Species, `ENA-RUN`)

write_delim(highest_crete_biodiversity_s,
            "results/highest_crete_biodiversity_sample_taxonomy.tsv",
            delim="\t")

## Phyla distribution

phyla_dist <- crete_biodiversity %>%
    group_by(sampleID,Phylum) %>%
    summarise(sum_abundance=sum(abundance), .groups="keep") %>%
    na.omit(Phylum)

total_phyla_dist <- phyla_dist %>% 
    group_by(Phylum) %>%
    summarise(sum_abundance=sum(sum_abundance),
              n_samples=n()) %>%
    mutate(sampleID="All samples") %>%
    arrange(desc(n_samples))


## create bar plots for each sample at family level, class level, Phylum etc
## shannon per sample
## bray curtis per sample
## 
