#!/usr/bin/Rscript

###############################################################################
# script name: isd_crete_spatial.R
# developed by: Savvas Paragkamian
# framework: ISD Crete
###############################################################################
# GOAL:
# Aim of this script is to use spatial data to enrich the metadata of the 
# ISD Crete project
#
###############################################################################
# OUTPUT:
#
###############################################################################
# usage:./isd_crete_spatial.R
###############################################################################
library(sf)
library(dplyr)
library(tibble)
library(readr)
library(magrittr)
library(tidyr)


# Metadata

master_metadata_old <- read.delim("Crete/Composite_MetaData_from_master.csv", sep=",")

metadata_long <- read_delim("ena_metadata/ena_isd_2016_attributes.tsv", delim="\t") %>%
    mutate(VALUE=gsub("\\r(?!\\n)","", VALUE, perl=T))
# metadata to wide format

metadata_wide <- metadata_long %>% 
    dplyr::select(-c(UNITS)) %>%
    mutate(TAG=gsub(" ","_", TAG, perl=T)) %>%
    pivot_wider(names_from=TAG, 
                values_from=VALUE)

metadata_wide$total_nitrogen <- as.numeric(metadata_wide$total_nitrogen)
metadata_wide$water_content <- as.numeric(metadata_wide$water_content)
metadata_wide$total_organic_carbon <- as.numeric(metadata_wide$total_organic_carbon)
metadata_wide$sample_volume_or_weight_for_DNA_extraction <- as.numeric(metadata_wide$sample_volume_or_weight_for_DNA_extraction)
metadata_wide$DNA_concentration <- as.numeric(metadata_wide$DNA_concentration)
metadata_wide$latitude <- as.numeric(metadata_wide$`geographic_location_(latitude)`)
metadata_wide$longitude <- as.numeric(metadata_wide$`geographic_location_(longitude)`)
metadata_wide$elevation <- as.numeric(metadata_wide$`geographic_location_(elevation)`)
metadata_wide$amount_or_size_of_sample_collected <- as.numeric(metadata_wide$amount_or_size_of_sample_collected)

metadata <- metadata_wide %>%
    select(`ENA-RUN`, source_material_identifiers, total_nitrogen, water_content,total_organic_carbon,sample_volume_or_weight_for_DNA_extraction,DNA_concentration,latitude, longitude,elevation, amount_or_size_of_sample_collected, vegetation_zone) %>%
    arrange(`ENA-RUN`) %>%
    as.data.frame()

rownames(metadata) <- metadata$`ENA-RUN`
metadata <- metadata[,-1]
