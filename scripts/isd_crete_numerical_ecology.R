#!/usr/bin/env Rscript

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
#library(mia)
library(phyloseq)
library(vegan)
library(dplyr)
library(tibble)
library(readr)
library(magrittr)
library(tidyr)
library(ggplot2)

################################## Load data ##################################
crete_biodiversity <- read_delim("results/crete_biodiversity_asv.tsv",delim="\t")
asv_metadata <- read_delim("results/asv_metadata.tsv",delim="\t")
phyla_samples_summary <- read_delim("results/phyla_samples_summary.tsv",delim="\t")
#crete_biodiversity_matrix <- readRDS("results/crete_biodiversity_matrix.RDS")
biodiversity_srs <- readRDS("results/biodiversity_srs.RDS")
tax_tab <- readRDS("results/tax_tab.RDS")

# Metadata

master_metadata_old <- read.delim("Crete/Composite_MetaData_from_master.csv", sep=",")

metadata <- read_delim("results/sample_metadata.tsv", delim="\t")


# differences of old and new data

master_metadata_old$team_site_location_id[which(!(master_metadata_old$team_site_location_id %in% metadata$source_material_identifiers))]

################################## functions ##################################
dist_long <- function(x,method){
    method <- method
    df <- as.data.frame(as.matrix(x)) %>%
    rownames_to_column() %>%
    pivot_longer(-rowname,
                 values_to=method,
                 names_to="colname")

    return(df)
}


####################### Destriptors Statistics ###############################
pw <- pairwise.wilcox.test(metadata$Observed, metadata_all$vegetation_zone, p.adjust.method="BH")
pw_e <- pairwise.wilcox.test(metadata$Observed, metadata_all$elevation_bin, p.adjust.method="BH")
pw_s <- pairwise.wilcox.test(metadata$Shannon, metadata_all$LABEL1, p.adjust.method="BH")
print(pw)

# correlations of diversity and other numerical metadata
with(metadata, cor(Shannon, water_content))


metadata_n <- metadata
rownames(metadata_n) <- metadata$ENA_RUN
nums <- unlist(lapply(metadata_n, is.numeric), use.names = FALSE)
metadata_n <- metadata_n[,c(nums)]
cc <- cor(metadata_n)

cc_sp <- cor(metadata_n, method="spearman")

write.table(cc_sp,
            "results/metadata_sprearman.tsv",
            sep="\t",
            row.names=T,
            col.names=NA)

############################## Community analysis ###########################
###################### Co-occurrence of samples and ASVs ####################

biodiversity_m <- biodiversity_srs
biodiversity_m[biodiversity_m > 0 ] <- 1
biodiversity_m <- as.matrix(biodiversity_m)

## matrix multiplication takes up a lot of memory and CPU, I had an error
## Error: vector memory exhausted (limit reached?)
## cd ~ ; touch .Renviron 
## echo R_MAX_VSIZE=100Gb >> .Renviron

asv_cooccur <- biodiversity_m %*% t(biodiversity_m)
sample_cooccur <- t(biodiversity_m) %*% biodiversity_m
isSymmetric(sample_cooccur) # is true so we can remove the lower triangle
sample_cooccur[lower.tri(sample_cooccur)] <- NA

sample_cooccur_l <- dist_long(sample_cooccur,"cooccurrence") %>%
    filter(rowname!=colname) %>%
    na.omit()


############################## Dissimilarity ###########################
# use the vegan package, the matrix must be transposed
biodiversity_srs_t <- t(biodiversity_srs)

bray <- vegdist(biodiversity_srs_t,
                method="bray")

jaccard <- vegdist(biodiversity_srs_t,
                method="jaccard",
                binary=TRUE)

aitchison <- vegdist(biodiversity_srs_t,
                method="robust.aitchison")


bray_l <- dist_long(bray, "bray")
jaccard_l <- dist_long(jaccard, "jaccard")
aitchison_l <- dist_long(aitchison, "robust.aitchison")

########################## Phylum level ########################
## Phyla distribution, average relative abundance and ubiquity
## Biogeography of soil bacteria and archaea across France

phyla_dist_samples <- phyla_samples_summary %>% 
    group_by(Phylum) %>%
    summarise(n_samples=n(),
              total_asvs=sum(asvs),
              total_reads_srs=sum(reads_srs_sum),
              average_relative=mean(relative_srs)) %>%
    arrange(desc(average_relative)) %>% 
    as.data.frame()

rownames(phyla_dist_samples) <- phyla_dist_samples$Phylum
phyla_dist_samples <- phyla_dist_samples[,-1]


######################## Site locations comparison #################
dissi_loc <- samples_locations %>%
    left_join(sample_cooccur_l,
              by=c("loc_1"="rowname", "loc_2"="colname")) %>%
    left_join(bray_l,
              by=c("loc_1"="rowname", "loc_2"="colname")) %>%
    left_join(jaccard_l,
              by=c("loc_1"="rowname", "loc_2"="colname")) %>%
    left_join(aitchison_l,
              by=c("loc_1"="rowname", "loc_2"="colname"))


summary(dissi_loc)

#####################################################################
z <- betadiver(biodiversity_srs_t, "z")
mod <- with(metadata_all, betadisper(z, LABEL1))
#sac <- specaccum(biodiversity_srs_t)

# Ordination
nmds <- vegan::metaMDS(biodiversity_srs_t,
                       k=2,
                       distance = "bray",
                       trymax=100)
stressplot(nmds)
ordiplot(nmds,display="sites", cex=1.25)
ordisurf(nmds,metadata_all$dem,main="",col="forestgreen")
ordihull(nmds,display="sites",label=T,  groups=metadata_all$LABEL1, cex=1.25)
ordihull(nmds,display="sites",label=T,  groups=metadata_all$elevation, cex=1.25)


tse <- runUMAP(tse,
               name = "UMAP",
               assay.type = "counts",
               ncomponents = 3)
plotReducedDim(tse, "UMAP",
               colour_by = "Group",
               ncomponents = c(1:3))
#### tests 

nmds <- vegan::metaMDS(t(crete_biodiversity_matrix),
                       k=2,
                       distance = "bray",
                       trymax=100)

stressplot(nmds)
ordiplot(nmds,display="sites", cex=1.25)

