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
source("scripts/functions.R")
library(vegan)
library(sf)
library(terra)
library(ape)
library(dplyr)
library(tibble)
library(readr)
library(magrittr)
library(tidyr)
library(ggplot2)
library(dendextend) 

################################## Load data ##################################
crete_biodiversity <- read_delim("results/crete_biodiversity_asv.tsv",delim="\t")
community_matrix_l <- read_delim("results/community_matrix_l.tsv",delim="\t")

community_matrix <- readRDS("results/community_matrix.RDS")
asv_metadata <- read_delim("results/asv_metadata.tsv", delim="\t")

# Metadata
metadata <- read_delim("results/sample_metadata.tsv", delim="\t")

######
metadata <- metadata %>% filter(ENA_RUN %in% rownames(community_matrix))

### 
print("samples with highest values of physicochemical properties")
metadata %>% arrange(desc(total_nitrogen)) %>% head(n=2) # ERR3697708 , ERR3697732
metadata %>% arrange(desc(total_organic_carbon)) %>% head(n=10) # ERR3697655, ERR3697675
metadata %>% arrange(desc(water_content)) %>% head(n=2) ## ERR3697703, ERR3697702 

################################# Metadata correlations #############################
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
######################################community matrix##########################
print("community matrix")

taxa <- community_matrix_l %>% distinct(Kingdom,Phylum,Class,Order,Family,Genus,Species,scientificName,classification)


############################## Dissimilarity ###########################
# use the vegan package, the matrix must be transposed
print("(dis)similarities")

########################## Phylum level ########################
## Phyla distribution, average relative abundance and ubiquity
## Biogeography of soil bacteria and archaea across France

#### Community matrix

bray <- vegdist(community_matrix,
                method="bray")

hc <- hclust(bray)

hc_df <- as.data.frame(cutree(hc,k=6)) |>
    rownames_to_column("ENA_RUN") |>
    as_tibble() 
colnames(hc_df) <- c("ENA_RUN", "cluster")

cluster_cols=c("#D55E00", "#F0E442","#009E73", "#56B4E9", "#BE81A3", "#999999")

cluster_cols=c("1"="#009E73",
              "2"="#56B4E9",
              "3"="#999999",
              "4"="#BE81A3",
              "5"="#F0E442",
              "6"="#D55E00")

dend <- as.dendrogram(hc) |>
    set("labels_col", value = cluster_cols, k=6) |>
    set("branches_k_color", value = cluster_cols, k = 6) 
    #color_branches(k = 6, col=cluster_cols) |>
    #color_labels(k = 6, col=cluster_cols)

png(file=paste0("figures/clustering_bray_hclust_samples.png"),
    width = 55,
    height = 20,
    res=300,
    units = "cm",
    bg="white")
plot(dend)
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

jaccard <- vegdist(community_matrix,
                method="jaccard",
                binary=TRUE)

aitchison <- vegdist(community_matrix,
                method="robust.aitchison")

jaccard_l <- dist_long(jaccard, "jaccard")
aitchison_l <- dist_long(aitchison, "robust.aitchison")

# beta diversity
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
    left_join(taxa, by=c("scientificName"="scientificName"))

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
################################# statistics ##########################

umap_isd_sites <- read_delim("results/umap_samples_2.tsv", delim="\t")
#umap_isd_sites_k1 <- read_delim("results/umap_samples_1.tsv", delim="\t")
#colnames(umap_isd_sites_k1) <- c("id", "UCIE")

metadata <- metadata |>
    left_join(umap_isd_sites, by=c("ENA_RUN"="id")) |>
    left_join(pcoa_bray_m) |>
    left_join(nmds_isd_sites)
#    left_join(umap_isd_sites_k1 ,by=c("ENA_RUN"="id"))
metadata$elevation_bin <- factor(metadata$elevation_bin,
                        levels=unique(metadata$elevation_bin)[order(sort(unique(metadata$elevation_bin)))])

############################# Statistics ###############################
##### regression
## diversity
cor.test(metadata$shannon, metadata$total_nitrogen)
cor.test(metadata$shannon, metadata$carbon_nitrogen_ratio) 

cor.test(metadata$shannon, metadata$elevation)
cor.test(metadata$shannon, metadata$water_content) 

gradient_scatterplot(metadata, "total_organic_carbon","shannon", "elevation_bin") 

####### Drivers numerical
cor.test(metadata$shannon, metadata$total_nitrogen)
cor.test(metadata$shannon, metadata$total_organic_carbon)
cor.test(metadata$shannon, metadata$carbon_nitrogen_ratio)
cor.test(metadata$shannon, metadata$water_content)
cor.test(metadata$shannon, metadata$elevation)
cor.test(metadata$shannon, metadata$aridity)
cor.test(metadata$shannon, metadata$bio_1)
cor.test(metadata$shannon, metadata$bio_12)

cor.test(metadata$shannon, metadata$UMAP1)

lm_s <- lm(metadata$shannon ~ metadata$bio_1 + metadata$geology_na+ metadata$total_organic_carbon)
summary(lm_s)
anova(lm_s)

### drivers of major axis of ordination
# first axis
lm_o <- lm(metadata$UMAP1 ~ metadata$bio_1 + metadata$total_organic_carbon  + metadata$geology_na)
summary(lm_o)
anova(lm_o)
cor.test(metadata$UMAP1, metadata$bio_1)
gradient_scatterplot(metadata, "bio_1","UMAP1", "none") 
gradient_scatterplot(metadata, "bio_12","UMAP1", "none") 
gradient_scatterplot(metadata, "total_organic_carbon","UMAP1", "none") 
gradient_scatterplot(metadata, "total_nitrogen","UMAP1", "none") 
cor.test(metadata$UMAP1, metadata$bio_12)
cor.test(metadata$UMAP1, metadata$total_organic_carbon)
cor.test(metadata$UMAP1, metadata$total_nitrogen)
kruskal.test(UMAP1 ~ LABEL3, data = metadata)  
kruskal.test(UMAP1 ~ geology_na, data = metadata)  

boxplot_single(metadata, "UMAP1", "geology_na", "bio_1")
# second axis
lm_o2 <- lm(metadata$UMAP2 ~ metadata$total_organic_carbon + metadata$water_content)
summary(lm_o2)
anova(lm_o2)
kruskal.test(UMAP2 ~ geology_na, data = metadata)
kruskal.test(UMAP2 ~ elevation_bin, data = metadata)
kruskal.test(UMAP2 ~ geology_na, data = metadata)
kruskal.test(UMAP2 ~ LABEL3, data = metadata)  

cor.test(metadata$UMAP2, metadata$total_organic_carbon)
cor.test(metadata$UMAP2, metadata$total_nitrogen)
cor.test(metadata$UMAP2, metadata$water_content)
gradient_scatterplot(metadata, "water_content","UMAP2", "none") 
gradient_scatterplot(metadata, "total_nitrogen","UMAP2", "none") 
gradient_scatterplot(metadata, "total_organic_carbon","UMAP2", "none") 
boxplot_single(metadata, "UMAP2","LABEL3", "total_organic_carbon")

####### Drivers categorical
kruskal.test(shannon ~ vegetation_zone, data = metadata)
kruskal.test(shannon ~ elevation_bin, data = metadata)
kruskal.test(shannon ~ aridity_class, data = metadata)
kruskal.test(shannon ~ LABEL2, data = metadata)
kruskal.test(shannon ~ LABEL3, data = metadata)

pairwise.wilcox.test(metadata$shannon, metadata$LABEL3, p.adjust.method="BH")

kruskal.test(shannon ~ geology_na, data = metadata)
pairwise.wilcox.test(metadata$shannon, metadata$geology_na, p.adjust.method="BH")

########### community dissimilarity tests #############
# calculate the bray dissimilatiry
bray <- vegdist(community_matrix)

# geology

# multivariate dispersion (variance) for a group of samples is to calculate
# the average distance of group members to the group centroid or spatial
# median (both referred to as 'centroid' from now on unless stated otherwise)
# in multivariate space. 
mod <- betadisper(bray, metadata$geology_na,type="centroid")
png("figures/community_betadisper_geology_box.png",
    res=300,
    width=60,
    height=40,
    unit="cm")
boxplot(mod)
dev.off()

## test to see if there are any significant differences 
anova(mod)
### Pairwise comparisons of group mean dispersions can also be performed using
### permutest.betadisper. An alternative to the classical comparison of group
### dispersions, is to calculate Tukey's Honest Significant Differences between
### groups, via TukeyHSD.betadisper. This is a simple wrapper to TukeyHSD. The
### user is directed to read the help file for TukeyHSD before using this
### function. In particular, note the statement about using the function with unbalanced designs.
permutest(mod, pairwise = TRUE, permutations = 99)
mod.HSD <- TukeyHSD(mod)

png("figures/community_betadisper_geology.png",
    res=300,
    width=60,
    height=40,
    unit="cm")
plot(mod.HSD)
dev.off()

# total nitrogen
mod <- betadisper(bray, metadata$total_nitrogen,type="centroid")
png("figures/community_betadisper_nitrogen_box.png",
    res=300,
    width=60,
    height=40,
    unit="cm")
plot(mod)
dev.off()

anova(mod)

permutest(mod, pairwise = TRUE, permutations = 99)
# label2

mod <- betadisper(bray, metadata$LABEL2,type="centroid")
png("figures/community_betadisper_label2_box.png",
    res=300,
    width=60,
    height=40,
    unit="cm")
boxplot(mod)
dev.off()

anova(mod)

permutest(mod, pairwise = TRUE, permutations = 99)
mod.HSD <- TukeyHSD(mod)
png("figures/community_betadisper_label2.png",
    res=300,
    width=60,
    height=40,
    unit="cm")
plot(mod.HSD)

dev.off()
# label3

mod <- betadisper(bray, metadata$LABEL3,type="centroid")
png("figures/community_betadisper_label3_box.png",
    res=300,
    width=60,
    height=40,
    unit="cm")
boxplot(mod)
dev.off()

anova(mod)

permutest(mod, pairwise = TRUE, permutations = 99)
mod.HSD <- TukeyHSD(mod)
png("figures/community_betadisper_label3.png",
    res=300,
    width=60,
    height=40,
    unit="cm")
plot(mod.HSD)

dev.off()
#plot(mod.HSD)
# elevation
mod <- betadisper(bray, metadata$elevation_bin,type="centroid")

png("figures/community_betadisper_elevation_box.png",
    res=300,
    width=60,
    height=40,
    unit="cm")
boxplot(mod)
dev.off()

anova(mod)

permutest(mod, pairwise = TRUE, permutations = 99)
mod.HSD <- TukeyHSD(mod)
png("figures/community_betadisper_elevation_bin.png",
    res=300,
    width=60,
    height=40,
    unit="cm")
plot(mod.HSD)

dev.off()


#### permanova

#adonis_elevation <- adonis2(community_matrix ~ elevation_bin, data=metadata_f, permutations=99)

adonis_multiple <- adonis2(community_matrix ~ bio_1*bio_12*elevation_bin*total_nitrogen*geology_na*LABEL3*carbon_nitrogen_ratio,
                           data=metadata,
                           permutations=999)


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
#community_matrix_m <- community_matrix
#community_matrix_m[community_matrix_m > 0] <- 1
#community_matrix_m <- as.matrix(community_matrix_m)

#sample_cooccur <- community_matrix_m %*% t(community_matrix_m)
#taxa_cooccur <- t(community_matrix_m) %*% community_matrix_m

#isSymmetric(taxa_cooccur) # is true so we can remove the lower triangle
#taxa_cooccur[lower.tri(taxa_cooccur)] <- NA

#taxa_cooccur_l <- dist_long(taxa_cooccur,"cooccurrence") %>%
#    filter(rowname!=colname) %>%
#    na.omit()

#write_delim(taxa_cooccur_l,"results/taxa_cooccur_l.tsv", delim="\t")

#isSymmetric(sample_cooccur) # is true so we can remove the lower triangle
#sample_cooccur[lower.tri(sample_cooccur)] <- NA

#sample_cooccur_l <- dist_long(sample_cooccur,"cooccurrence") %>%
#    filter(rowname!=colname) %>%
#    na.omit() %>% 
#    left_join(bray_l,
#              by=c("rowname"="rowname", "colname"="colname")) %>%
#    left_join(jaccard_l,
#              by=c("rowname"="rowname", "colname"="colname")) %>%
#    left_join(aitchison_l,
#              by=c("rowname"="rowname", "colname"="colname"))


#write_delim(sample_cooccur_l,"results/sample_cooccur_l.tsv", delim="\t")


######################## Site locations comparison ASV #################
#samples_locations <- metadata %>%
#    pivot_wider(id_cols=sites,
#                names_from=location,
#                values_from=ENA_RUN)
#
#
#dissi_loc <- samples_locations %>%
#    left_join(sample_cooccur_l,
#              by=c("loc_1"="rowname", "loc_2"="colname"))
#
#summary(dissi_loc)
#
########################################################################
###########------- Spatial Summary Table ------------- #################
########################################################################

##  load map data
crete_shp <- sf::st_read("spatial_data/crete/crete.shp")

crete_area <- data.frame(name="Crete",
                         area=sum(units::set_units(st_area(crete_shp), km^2)),map="GDAM")
# ---------------------------#
### protected areas
# ---------------------------#
wdpa_crete <- sf::st_read("spatial_data/wdpa_crete/wdpa_crete.shp") |>
    mutate(DESIG_ENG = gsub("Wildlife Refugee", "Wildlife Refuge", DESIG_ENG))

wdpa_crete_wildlife <- wdpa_crete |>
    filter(DESIG_ENG=="Wildlife Refuge")

wdpa_crete$area <- units::set_units(st_area(wdpa_crete),km^2)

wdpa_crete_all <- data.frame(name="total protected", 
                             area=sum(wdpa_crete$area))

wdpa_crete_combine <- st_union(wdpa_crete) %>%
    st_make_valid() %>%
    st_as_sf() %>%
    filter(st_geometry_type(.) %in% c("MULTIPOLYGON"))

wdpa_crete_combine_area <- data.frame(name="total protected (no overlap)",
                                      area=sum(units::set_units(st_area(wdpa_crete_combine), km^2)))


protected_area <- wdpa_crete |> 
    group_by(DESIG_ENG) |>
    summarise(area=sum(area)) |>
    st_drop_geometry() |>
    dplyr::rename("name"="DESIG_ENG") |>
    bind_rows(wdpa_crete_all, wdpa_crete_combine_area) |>
    arrange(area) |>
    mutate(area=round(area,2)) |>
    mutate(map="wdpa")

# ---------------------------#
# Corine Land Cover Label2
# ---------------------------#

clc_crete_shp <- st_read("spatial_data/clc_crete_shp/clc_crete_shp.shp")

clc_crete_shp$area <- units::set_units(st_area(clc_crete_shp),km^2)

clc_crete_shp_area <- clc_crete_shp |>
    st_drop_geometry() |>
    group_by(LABEL2) |>
    summarise(area=sum(area)) |>
    dplyr::rename("name"="LABEL2") |>
    mutate(map="CLC")

# ---------------------------#
# geology #
# ---------------------------#
crete_geology <- st_read("spatial_data/crete_geology/crete_geology.shp") |>
    st_make_valid()

geol_df <- data.frame(
  symbol = c(
    "J-E",
    "K-E",
    "K.k",
    "K.m",
    "Mk",
    "Mm.I",
    "Ph-T",
    "Q.al",
    "T.br",
    "f",
    "fo",
    "ft",
    "o"
  ),
  name = c(
    "Plattenkalk",
    "Limestone of Pindos",
    "Limestone of Tripolis",
    "Carbonaceous Allochthonous",
    "Neogenic",
    "Neogenic",
    "Phyllites - Chalazites",
    "Quaternary alluvium",
    "Limestone of Tripolis",
    "Schale",
    "Schale of Pindos",
    "Flysch of Tripolis",
    "Ophiolite Complex Allochthonous"
  ),
  stringsAsFactors = FALSE
)
# area calculation
crete_geology$area <- units::set_units(st_area(crete_geology),km^2)

crete_geology_area <- crete_geology |>
    st_drop_geometry() |>
    group_by(geology_na) |>
    summarise(area=sum(area)) |>
    #left_join(geol_df,by=c("geology_na"="symbol")) |>
    na.omit() |>
    mutate(name=geology_na) |>
    #mutate(name=paste0(name," (",geology_na,")")) |>
    dplyr::select(-geology_na) |>
    mutate(map="Geology")

# ---------------------------#
# hilda 1976 - 2016
# ---------------------------#

hilda <- rast("spatial_data/hilda_1976_2016/hilda_1976_2016.tif")
# make 0 as NA 
hilda[hilda == 0] <- NA
# the size of each cell
hilda_a <- cellSize(hilda, unit="km")
# calculation per category
hilda_area <- zonal(hilda_a, hilda, fun="sum", na.rm=F)

hilda_id_names <- read_delim("spatial_data/hilda_1976_2016/hilda_transitions_names.tsv", delim="\t")

hilda_area <- hilda_area |>
    left_join(hilda_id_names, by=c("lyr.1"="hilda_id")) |>
    dplyr::select(-lyr.1) |>
    mutate(area=round(units::set_units(area,km^2),digits=0)) |>
    mutate(map="Hilda+ 1976-2016") |>
    rename("name"="hilda_name")

# ---------------------------#
# raster global aridity index
# ---------------------------#
aridity_crete <- rast("spatial_data/crete_aridity_index.tif")
aridity_crete[aridity_crete[] == 0 ] = NA

# Define the breaks and labels
brks <- c(0, 300, 2000, 5000, 6500, 9000)
labs <- c("Hyper Arid", "Arid", "Semi-Arid", "Dry sub-humid", "Humid")

# Reclassify using cut() through app()
aridity_class <- app(aridity_crete, \(x)
  cut(x,
      breaks = brks,
      include.lowest = TRUE,
      right = FALSE,
      labels = FALSE)
)
# make factor
aridity_class <- as.factor(aridity_class)
# add the labels
levels(aridity_class) <- data.frame(ID = 1:length(labs), label = labs)

# the size of each cell
aridity_class_a <- cellSize(aridity_class, unit="km")
# calculation per category
aridity_class_area <- zonal(aridity_class_a, aridity_class, fun="sum", na.rm=F)

aridity_class_area <- aridity_class_area |>
    mutate(map="Aridity") |>
    mutate(area=round(units::set_units(area,km^2),digits=0)) |>
    rename("name"="label")

# ---------------------------#
# raster desertification and Environmental Sensitive Areas Index
# ---------------------------#
desertification_crete <- rast("spatial_data/crete_desertification_risk/esa3rdp_crete.tif")
# the size of each cell
desertification_crete_a <- cellSize(desertification_crete, unit="km")
# calculation per category
desertification_crete_area <- zonal(desertification_crete_a, desertification_crete, fun="sum", na.rm=F)

desertification_crete_area <- desertification_crete_area |> 
    mutate(map="Desertification") |>
    mutate(area=round(units::set_units(area,km^2),digits=0)) |>
    rename("name"="ESA_12CL")

# ---------------------------#
# soil erosion G2
# ---------------------------#
soil_erosion_g2 <- rast("spatial_data/crete_soil_erosion_g2.tif")
# the size of each cell
soil_erosion_g2_a <- cellSize(soil_erosion_g2, unit="km")
# calculation per category
g2_area <- zonal(soil_erosion_g2_a, soil_erosion_g2, fun="sum", na.rm=F)

g2_area <- g2_area |>
    mutate(map="Soil Erosion G2") |>
    mutate(area=round(units::set_units(area,km^2),digits=0)) |>
    rename("name"="crete_soil_erosion_g2")

# ---------------------------#
# harmonised world soil database v2
# ---------------------------#
hwsd2 <- rast("spatial_data/hwsd2_crete/hwsd2_crete.tif")
# the size of each cell
hwsd2_a <- cellSize(hwsd2, unit="km")
# calculation per category
hwsd2_area <- zonal(hwsd2_a, hwsd2, fun="sum", na.rm=F)

# hswd metadata
# with trimws the leading spaces are removed for the values.
HWSD2_wrb4 <- read_delim("spatial_data/hwsd2_crete/HWSD2_D_WRB4.tsv", delim="\t") |>
    mutate(VALUE=trimws(VALUE)) |>
    distinct(VALUE, CODE) 

HWSD2_SMU <- read_delim("spatial_data/hwsd2_crete/HWSD2_SMU.tsv", delim="\t") |>
    distinct(HWSD2_SMU_ID, WRB4) |>
    left_join(HWSD2_wrb4, by=c("WRB4"="CODE"))

hwsd2_area <-hwsd2_area |>
    left_join(HWSD2_SMU, by=c("HWSD2"="HWSD2_SMU_ID")) |>
    mutate(area=round(units::set_units(area,km^2),digits=0)) |>
    na.omit() |>
    group_by(VALUE) |>
    summarise(area=sum(area)) |>
    mutate(map="HWSD2") |>
    rename("name"="VALUE")

###################### calculate summary ##################
maps <- list(protected_area,
             clc_crete_shp_area,
             crete_geology_area,
             hilda_area,
             aridity_class_area,
             desertification_crete_area,
             g2_area,
             hwsd2_area,
             crete_area)

maps_area <- bind_rows(maps)

write_delim(maps_area, "results/maps_area_summary.tsv",delim="\t")

###################### sample summary ##################
classes <- c("LABEL2", "geology_na", "HWSD2_value","ESA_12CL", "aridity_class","hilda_name","erosion_g2","SPA","SAC","Wildlife_Refuge")

classes_samples <- list()

for (i in seq_along(classes)){
    print(i)
    classes_samples[[i]] <- metadata |> 
        group_by(metadata[[classes[i]]]) |>
        summarise(samples=n(),
                  taxa_richness=sum(taxa),
                  asv_richness=sum(asvs),
                  mean_shannon=mean(shannon),
                  sd_shannon=sd(shannon)) |>
        na.omit()
}

samples_total <- do.call(rbind, classes_samples) 
colnames(samples_total)[1] <- c("name")

area_total <- maps_area |>
    left_join(samples_total) |>
    mutate(area=round(area, digits=1),
           mean_shannon=round(mean_shannon, digits=2),
           sd_shannon=round(sd_shannon, digits=2)) |>
    mutate(across(everything(), ~ replace_na(.x, 0)))

write_delim(area_total, "results/data_cube_summary_table.tsv", delim="\t")

#area_total |> arrange(category, class, area) |> kbl() |> kable_styling(latex_options = "scale_down")

print("finish")
