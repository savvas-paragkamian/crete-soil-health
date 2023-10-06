#!/usr/bin/env Rscript

###############################################################################
# script name: isd_crete_biodiversity.R
# developed by: Savvas Paragkamian, Johanna Holms
# framework: ISD Crete
###############################################################################
# GOAL:
# Aim of this script is provide a streamline of files with all the information
# required for further analyses. That is to use ASVs and the taxonomy to filter
# and normalise the reads of each sample biodiversity and to combine sample 
# metadata.
#
###############################################################################
# OUTPUT:
#
# the main output of this script is 
# 1. crete_biodiversity_asv.tsv
# which contains the sample - asv occurrences with abundance along with 
# taxonomic information.
# 2. sample_metadata.tsv
# 3. asv_metadata.tsv
# 
# Other 3 files are produced
# crete_biodiversity_matrix.RDS, a matrix of abundances
# tax_tab.RDS, taxonomy table with the remaining asvs
# sample_stats.tsv
###############################################################################
# RUNNING TIME: 9 minutes
###############################################################################
# usage:./scripts/isd_crete_biodiversity.R
###############################################################################
library(magrittr)
library(dplyr)
library(tibble)
library(readr)
library(tidyr)
library(vegan)
library(SRS) # normalisation
#library(ANCOMBC)
################################## functions ##################################
# this function keeps the last occurrence of a string separated by |
keep_last <- function(x) tail(strsplit(x, split="\\|")[[1]],1)

dist_long <- function(x,method){
    method <- method
    df <- as.data.frame(as.matrix(x)) %>%
    rownames_to_column() %>%
    pivot_longer(-rowname,
                 values_to=method,
                 names_to="colname")

    return(df)
}
################################## Load data ##################################
# Load data old repo
#abundance_asv_old <- readRDS("Crete/all_runs_dada2_abundance_table.rds")
#master_metadata_old <- read.delim("Crete/Composite_MetaData_from_master.csv", sep=",")

# Load data
abundance_asv_long <- read_delim("dada2_output/taxonomy/seqtab_nochim_long.tsv", delim="\t") %>% 
    mutate(ENA_RUN=gsub("_1_filt.fastq.gz", "", file, perl=T))
#abundance_asv <- readRDS("dada2_output/taxonomy/seqtab_nochim.RDS") # use this one if data are RDS
asv_fasta <- read_delim("dada2_output/taxonomy/asv_fasta_ids.tsv", delim="\t")
taxa_asv <- readRDS("dada2_output/taxonomy/dada2_taxonomy.RDS")
species_asv <- readRDS("dada2_output/taxonomy/dada2_taxa_species.RDS")
metadata <- read_delim("results/metadata_spatial.tsv", delim="\t")

## location pairs of each site
metadata$sites <- gsub("_loc_*.","", metadata$source_material_identifiers)
metadata$location <- gsub(".*(loc_.*)","\\1", metadata$source_material_identifiers)

########################## transform the matrices to long tables ##########################

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

####################### merge abundance and taxonomy ##########################

crete_biodiversity_all <- abundance_asv_long %>%
    left_join(taxa_asv_all,by=c("asv_id"="asv_id")) %>%
    filter(abundance>0)

######################## Decontamination #####################################
## very important step, we need the samples controls.

# This is the MASTER dataset of Crete biodiversity and taxonomy


########################## Normalisation Compositionality #####################
## total ASVs and taxonomy
asv_no_taxonomy <- crete_biodiversity_all %>% 
    distinct(asv,classification) %>%
    filter(is.na(classification))

print(paste0("ASVs without taxonomy: ", nrow(asv_no_taxonomy)))
## create a abundance matrix
crete_biodiversity_m <- crete_biodiversity_all %>%
    filter(!is.na(classification)) %>%
    select(ENA_RUN, asv_id, abundance) %>%
    pivot_wider(names_from=ENA_RUN, values_from=abundance, values_fill = 0) %>%
    as.matrix()

crete_biodiversity_matrix <- crete_biodiversity_m[,-1]
crete_biodiversity_matrix <- apply(crete_biodiversity_matrix, 2, as.numeric)
rownames(crete_biodiversity_matrix) <- crete_biodiversity_m[,1]
saveRDS(crete_biodiversity_matrix, "results/crete_biodiversity_matrix.RDS")

# SRS compositional preparation

Cmin <- min(colSums(crete_biodiversity_matrix))
Cmax <- max(colSums(crete_biodiversity_matrix))
summary(colSums(crete_biodiversity_matrix))
table(colSums(crete_biodiversity_matrix) > 10000)

## SRS curve
jpeg(file="results/isd_crete_srs_curve.jpeg")

SRScurve(crete_biodiversity_matrix,
         metric = "richness",
         step = 500,
         xlim=c(0,Cmax),
         xlab = "sample size",
         ylab = "richness",
         label = F,
         col = c("#5A4A6F", "#E47250",  "#EBB261"))

dev.off()

crete_biodiversity_matrix_df <- as.data.frame(crete_biodiversity_matrix)
rownames(crete_biodiversity_matrix_df) <- rownames(crete_biodiversity_matrix)

biodiversity_srs <- SRS(crete_biodiversity_matrix_df,
                        10000,
                        set_seed = TRUE,
                        seed = 1)

rownames(biodiversity_srs) <- rownames(crete_biodiversity_matrix_df)

# how many samples don't have ASVs
length(which(colSums(biodiversity_srs)==0))
# how many ASVs have 0 abundance after the SRS
length(which(rowSums(biodiversity_srs)==0))
dim(biodiversity_srs)
## filter empty
biodiversity_srs <- biodiversity_srs[-which(rowSums(biodiversity_srs)==0) ,]

#### Save ######
saveRDS(biodiversity_srs, "results/biodiversity_srs.RDS")

### remove the asvs that don't have a taxonomy
############################# Merge SRS ################################
biodiversity_srs_l <- dist_long(biodiversity_srs,"srs_abundance")

crete_biodiversity <- crete_biodiversity_all %>%
    filter(!is.na(classification)) %>% 
    left_join(biodiversity_srs_l, by=c("asv_id"="rowname", "ENA_RUN"="colname"))

write_delim(crete_biodiversity,"results/crete_biodiversity_asv.tsv",delim="\t")

########################## Samples diversity #########################
################################# Indices ###################################
# total ASV per sample
samples_srs_asvs <- biodiversity_srs_l %>%
    group_by(colname) %>%
    filter(srs_abundance>0) %>%
    summarise(n_srs_asv=n())

# Explore alpha-metrics
biodiversity_srs_t <- t(biodiversity_srs)

shannon <- data.frame(shannon = as.matrix(diversity(biodiversity_srs_t, index = "shannon")))
observed <- data.frame(t(estimateR(biodiversity_srs_t)))

biodiversity_index <- cbind(shannon,observed)

biodiversity_index$ENA_RUN <- rownames(biodiversity_index)
rownames(biodiversity_index) <- NULL
## taxonomic, asv and read diversity per sample

sample_taxonomy_stats <- crete_biodiversity %>% 
    group_by(ENA_RUN, classification, scientificName) %>% 
    summarise(asvs=n(),
              reads=sum(abundance),
              reads_srs=sum(srs_abundance,na.rm=T),
              .groups="keep") %>%
    group_by(ENA_RUN,classification) %>%
    summarise(taxa=n(),
              reads=sum(reads),
              asvs=sum(asvs),
              reads_srs=sum(reads_srs),
              .groups="keep") %>%
#    pivot_wider(names_from=classification,values_from=n_taxa) %>%
    dplyr::ungroup()

write_delim(sample_taxonomy_stats,
            "results/sample_taxonomy_stats.tsv",
            delim="\t")

sample_stats_total <- sample_taxonomy_stats %>%
    group_by(ENA_RUN) %>%
    summarise(taxa=sum(taxa),
              reads=sum(reads),
              asvs=sum(asvs),
              reads_srs=sum(reads_srs))

## keep only the sample metadata after filtering
## filter also metadata and taxonomy
metadata_all <-  metadata %>%
    left_join(biodiversity_index, by=c("ENA_RUN"="ENA_RUN")) %>%
    left_join(sample_stats_total, by=c("ENA_RUN"="ENA_RUN")) %>%
    left_join(samples_srs_asvs, by=c("ENA_RUN"="colname"))

write_delim(metadata_all,
            "results/sample_metadata.tsv",
            delim="\t")

########################## ASV summary ###############################
# Create taxonomy table of the remaining asvs
tax_tab1 <- crete_biodiversity %>%
    distinct(asv_id, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
    as.matrix()

tax_tab <- tax_tab1[,-1]
rownames(tax_tab) <- tax_tab1[,1]
saveRDS(tax_tab, "results/tax_tab.RDS")

asv_stats <- crete_biodiversity %>% 
    group_by(asv_id, classification, scientificName) %>% 
    summarise(n_samples=n(),
              reads=sum(abundance),
              reads_srs=sum(srs_abundance, na.rm=T),
              reads_srs_mean=mean(srs_abundance, na.rm=T),
              reads_srs_sd=sd(srs_abundance, na.rm=T),
              .groups="keep") %>%
    dplyr::ungroup()

write_delim(asv_stats,"results/asv_metadata.tsv",delim="\t")

## singletons
singletons <- asv_stats %>% 
    filter(reads==1) %>%
    nrow()

print(paste0("there are ",singletons," singletons asvs"))

## asv and samples distribution
asv_sample_dist <- asv_stats %>%
    group_by(n_samples) %>%
    summarise(n_asv=n(),
              reads=sum(reads),
              reads_srs=sum(reads_srs, na.rm=T))

### highest species biodiversity sample

crete_biodiversity_s <- crete_biodiversity %>%
    distinct(ENA_RUN, Species) %>% 
    group_by(ENA_RUN) %>%
    summarise(Species=n())

print("sample with the highest microbial species diversity")
crete_biodiversity_s[which(crete_biodiversity_s$Species==max(crete_biodiversity_s$Species)),]

