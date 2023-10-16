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
#library(phyloseq)
library(vegan)
library(ape)
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
genera_phyla_samples <- read_delim("results/genera_phyla_samples.tsv",delim="\t")
genera_phyla_stats <- read_delim("results/genera_phyla_stats.tsv",delim="\t")
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

######
###
metadata <- metadata %>% filter(ENA_RUN %in% colnames(biodiversity_srs))
####################### Destriptors Statistics ###############################
print("descriptors stats")
pw <- pairwise.wilcox.test(metadata$S.obs, metadata$vegetation_zone, p.adjust.method="BH")
pw_e <- pairwise.wilcox.test(metadata$S.obs, metadata$elevation_bin, p.adjust.method="BH")
pw_s <- pairwise.wilcox.test(metadata$shannon, metadata$LABEL1, p.adjust.method="BH")
print(pw)

# correlations of diversity and other numerical metadata
with(metadata, cor(shannon, water_content))


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


############################## Dissimilarity ###########################
# use the vegan package, the matrix must be transposed
print("(dis)similarities")
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

#####################################################################
#################### ASV community matrix ###################
#z <- betadiver(biodiversity_srs_t, "z")
#mod <- with(metadata, betadisper(z, LABEL1))
#sac <- specaccum(biodiversity_srs_t)

# Ordination
#nmds <- vegan::metaMDS(biodiversity_srs_t,
#                       k=2,
#                       distance = "bray",
#                       trymax=100)
#stressplot(nmds)
#ordiplot(nmds,display="sites", cex=1.25)
#ordisurf(nmds,metadata$dem,main="",col="forestgreen")
#ordihull(nmds,display="sites",label=T,  groups=metadata$LABEL1, cex=1.25)
#ordihull(nmds,display="sites",label=T,  groups=metadata$elevation, cex=1.25)


#tse <- runUMAP(tse,
#               name = "UMAP",
#               assay.type = "counts",
#               ncomponents = 3)
#plotReducedDim(tse, "UMAP",
#               colour_by = "Group",
#               ncomponents = c(1:3))
#### tests 

#nmds <- vegan::metaMDS(t(crete_biodiversity_matrix),
#                       k=2,
#                       distance = "bray",
#                       trymax=100)

#stressplot(nmds)
#ordiplot(nmds,display="sites", cex=1.25)

############### genera as taxa in the community matrix################
print("genera as community")
genera_samples_m <- genera_phyla_samples %>%
    dplyr::select(ENA_RUN,relative_srs,Genus) %>%
    pivot_wider(names_from=Genus,
                values_from=relative_srs,
                values_fill=0) %>%
    as.data.frame()

write_delim(genera_samples_m,"results/genera_samples_matrix.tsv", delim="\t")

community_matrix <- genera_samples_m
rownames(community_matrix) <- community_matrix[,1]
community_matrix <- community_matrix[,-1]

genera_tax <- genera_phyla_samples %>% distinct(Phylum, Genus)

#### Community matrix with genera

bray <- vegdist(community_matrix,
                method="bray")

png(file="figures/clustering_bray_hclust_samples.png",
    width = 50,
    height = 30,
    res=300,
    units = "cm",
    bg="white")
plot(hclust(bray))
dev.off()

bray_tax <- vegdist(t(community_matrix),method="bray")

png(file="figures/clustering_hclust_taxa.png",
    width = 50,
    height = 50,
    res=300,
    units = "cm",
    bg="white")
plot(hclust(bray_tax))
dev.off()

bray_samples <- vegdist(community_matrix,method="bray")
#homoscedasticity_s <- betadisper(bray_samples, metadata$LABEL1, type = c("median","centroid"), bias.adjust = FALSE)

bray_l <- dist_long(bray, "bray")

z <- betadiver(community_matrix, "z")
#mod <- with(metadata, betadisper(z, LABEL1))
#sac <- specaccum(biodiversity_srs_t)

######################### Ordination ############################
####################### PCoA #########################
print("starting PCoA")

#### sites
pcoa_bray <- ape::pcoa(bray)
pcoa_bray_m <- pcoa_bray$vectors %>% as.data.frame() %>% rownames_to_column("ENA_RUN")

write_delim(pcoa_bray_m,"results/ordination_pcoa_bray_sites.tsv", delim="\t")


####################### nMDS #########################
print("starting nMDS")
nmds_isd <- vegan::metaMDS(community_matrix,
                       k=2,
                       distance = "bray",
                       trymax=100)

# fit environmental numerical vectors
env_isd <- metadata %>%
    filter(ENA_RUN %in% rownames(community_matrix)) %>% 
    column_to_rownames(var="ENA_RUN")# %>%

print("starting envfit")
envfit_isd <- envfit(nmds_isd, env_isd, permutations = 999, na.rm=T) 
env_scores_isd <- as.data.frame(scores(envfit_isd, display = "vectors"))
write_delim(env_scores_isd,"results/env_scores_isd.tsv", delim="\t")

# plotting
png(file="figures/ordination_nmds_stressplot.png",
    width = 30,
    height = 30,
    res=300,
    units = "cm",
    bg="white")
stressplot(nmds_isd)
dev.off()

png(file="figures/ordination_nmds_sites_lat.png",
    width = 30,
    height = 30,
    res=300,
    units = "cm",
    bg="white")
ordiplot(nmds_isd,display="sites", cex=1.25)
ordisurf(nmds_isd,env_isd$latitude,main="",col="firebrick") ## interesting
#ordisurf(nmds,metadata$dem,main="",col="orange")
dev.off()

png(file="figures/ordination_nmds_sites_dem.png",
    width = 30,
    height = 30,
    res=300,
    units = "cm",
    bg="white")
ordiplot(nmds_isd,display="sites", cex=1.25)
ordisurf(nmds_isd,env_isd$dem,main="",col="firebrick") ## interesting
dev.off()


nmds_isd_taxa <- as.data.frame(scores(nmds_isd, "species")) %>%
    rownames_to_column("scientificName") %>%
    left_join(genera_tax, by=c("scientificName"="Genus"))

write_delim(nmds_isd_taxa,"results/nmds_isd_taxa.tsv", delim="\t")

nmds_isd_sites <- as.data.frame(scores(nmds_isd,"sites")) %>%
    rownames_to_column("ENA_RUN") %>%
    left_join(metadata[,c("ENA_RUN","elevation_bin", "LABEL1", "LABEL2", "vegetation_zone")],
              by=c("ENA_RUN"="ENA_RUN"))

write_delim(nmds_isd_sites,"results/nmds_isd_sites.tsv", delim="\t")

############################ nmds k3 ###########################

#nmds_isd_k3 <- vegan::metaMDS(community_matrix,
#                       k=3,
#                       distance = "bray",
#                       trymax=100)
#nmds_isd_taxa_k3 <- as.data.frame(scores(nmds_isd_k3,"species"))
#nmds_isd_sites_k3 <- as.data.frame(scores(nmds_isd_k3,"sites"))
############################# dbRDA ############################

#dbrda_isd <- dbrda(community_matrix ~ elevation + latitude + longitude + total_organic_carbon + total_nitrogen + water_content,env_isd, dist="bray")

############################# UMAP ############################
# the python script isd_crete_umap.py
# performs the UMAP algorithm

############################## Community analysis ###########################
###################### Co-occurrence of samples and ASVs ####################
print("starting co-occurrence")

#biodiversity_m <- biodiversity_srs
#biodiversity_m[biodiversity_m > 0 ] <- 1
#biodiversity_m <- as.matrix(biodiversity_m)

## matrix multiplication takes up a lot of memory and CPU, I had an error
## Error: vector memory exhausted (limit reached?)
## cd ~ ; touch .Renviron 
## echo R_MAX_VSIZE=200Gb >> .Renviron

#asv_cooccur <- biodiversity_m %*% t(biodiversity_m)
community_matrix_m <- community_matrix
community_matrix_m[community_matrix_m > 0] <- 1
community_matrix_m <- as.matrix(community_matrix_m)

sample_cooccur <- community_matrix_m %*% t(community_matrix_m)
taxa_cooccur <- t(community_matrix_m) %*% community_matrix_m

isSymmetric(taxa_cooccur) # is true so we can remove the lower triangle
taxa_cooccur[lower.tri(taxa_cooccur)] <- NA

taxa_cooccur_l <- dist_long(taxa_cooccur,"cooccurrence") %>%
    filter(rowname!=colname) %>%
    na.omit()

write_delim(taxa_cooccur_l,"results/taxa_cooccur_l.tsv", delim="\t")

isSymmetric(sample_cooccur) # is true so we can remove the lower triangle
sample_cooccur[lower.tri(sample_cooccur)] <- NA

sample_cooccur_l <- dist_long(sample_cooccur,"cooccurrence") %>%
    filter(rowname!=colname) %>%
    na.omit() %>% 
    left_join(bray_l,
              by=c("rowname"="rowname", "colname"="colname")) %>%
    left_join(jaccard_l,
              by=c("rowname"="rowname", "colname"="colname")) %>%
    left_join(aitchison_l,
              by=c("rowname"="rowname", "colname"="colname"))


write_delim(sample_cooccur_l,"results/sample_cooccur_l.tsv", delim="\t")
######################## Site locations comparison ASV #################
samples_locations <- metadata %>%
    pivot_wider(id_cols=sites,
                names_from=location,
                values_from=ENA_RUN)


dissi_loc <- samples_locations %>%
    left_join(sample_cooccur_l,
              by=c("loc_1"="rowname", "loc_2"="colname"))

summary(dissi_loc)
print("finish")
