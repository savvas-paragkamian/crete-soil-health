#!/usr/bin/env Rscript
###############################################################################
# script name: figures.R
# developed by: Savvas Paragkamian
# framework: ISD Crete
###############################################################################
# GOAL:
# Aim of this script is to load the results and with some minor counts and 
# transformations make publication-quality figures. That includes maps, bar 
# plots, scatterplots etc.
###############################################################################
# OUTPUT:
#
###############################################################################
# usage:./scripts/figures.R
###############################################################################

# load packages and functions
#setwd("../")
source("scripts/functions.R")

library(vegan)
library(tidyverse)
library(ggnewscale)
library(ggpubr)
library(pheatmap)
library(sf)
library(jpeg)
library(raster)
library(scales)
library(ucie)

################################## Load data ##################################
## biodiversity
community_matrix_l <- read_delim("results/community_matrix_l.tsv",delim="\t")
crete_biodiversity <- read_delim("results/crete_biodiversity_asv.tsv",delim="\t")
asv_metadata <- read_delim("results/asv_metadata.tsv",delim="\t")
genera_phyla_stats <- read_delim("results/genera_phyla_stats.tsv",delim="\t")
genera_phyla_samples <- read_delim("results/genera_phyla_samples.tsv",delim="\t")

# Metadata
metadata <- read_delim("results/sample_metadata.tsv", delim="\t")

metadata$elevation_bin <- factor(metadata$elevation_bin,
                        levels=unique(metadata$elevation_bin)[order(sort(unique(metadata$elevation_bin)))])

samples_ucie_nmds_genera <- read_delim("results/samples_ucie_nmds_genera.tsv")

## spatial
locations_spatial <- metadata %>%
    st_as_sf(coords=c("longitude", "latitude"),
             remove=F,
             crs="WGS84")

crete_shp <- sf::st_read("spatial_data/crete/crete.shp")
crete_peaks <- read_delim("spatial_data/crete_mountain_peaks.csv", delim=";", col_names=T) %>%
    st_as_sf(coords=c("X", "Y"),
             remove=F,
             crs="WGS84")

clc_crete_shp <- st_read("spatial_data/clc_crete_shp/clc_crete_shp.shp")
crete_geology <- st_read("spatial_data/crete_geology/crete_geology.shp")
natura_crete <- sf::st_read("spatial_data/natura2000/natura2000_crete.shp")
wdpa_crete <- sf::st_read("spatial_data/wdpa_crete/wdpa_crete.shp") %>% filter(DESIG_ENG=="Wildlife Refugee") %>%
    mutate(DESIG_ENG = gsub("Wildlife Refugee", "Wildlife Refuge", DESIG_ENG))

natura_crete_land <- st_intersection(natura_crete, crete_shp)
natura_crete_land_sci <- natura_crete_land %>% filter(SITETYPE=="B")
# raster DEM hangling
dem_crete <- raster("spatial_data/dem_crete/dem_crete.tif")
dem_crete_pixel <- as(dem_crete, "SpatialPixelsDataFrame")
dem_crete_df <- as.data.frame(dem_crete_pixel) %>% filter(dem_crete>0)

####################### UCIE #########################
# UCIE needs 3 axis of ordination
print("starting UCIE")

############################### ucie with umap ########################
# sites
#umap_isd_sites_k3 <- read_delim("results/umap_samples_3.tsv", delim="\t") %>% column_to_rownames("id")

#umap_isd_genera_k3 <- read_delim("results/umap_genera_3.tsv", delim="\t") %>% column_to_rownames("id")

#umap_isd_sites_ucie <- ucie::data2cielab(umap_isd_sites_k3, LAB_coordinates = F)
#colnames(umap_isd_sites_ucie) <- c("ENA_RUN","UCIE")
#write_delim(umap_isd_sites_ucie,"results/umap_isd_sites_ucie.tsv", delim="\t")
# taxa
#umap_isd_taxa_ucie <- ucie::data2cielab(umap_isd_genera_k3)
#colnames(umap_isd_taxa_ucie) <- c("scientificName","UCIE")
#write_delim(umap_isd_taxa_ucie,"results/umap_isd_taxa_ucie.tsv", delim="\t")

############################### ucie with pcoa ########################

pcoa_isd_sites <- read_delim("results/ordination_pcoa_bray_sites.tsv", delim="\t") %>%
    column_to_rownames("ENA_RUN") %>% dplyr::select(Axis.1, Axis.2, Axis.3)

#pcoa_isd_sites_ucie <- ucie::data2cielab(pcoa_isd_sites, LAB_coordinates = F)
#colnames(pcoa_isd_sites_ucie) <- c("ENA_RUN","UCIE")
#write_delim(pcoa_isd_sites_ucie,"results/pcoa_isd_sites_ucie.tsv", delim="\t")
pcoa_isd_sites_ucie <- read_delim("results/pcoa_isd_sites_ucie.tsv")
############################### ucie with pcoa ########################
#nmds_isd_taxa_ucie <- data2cielab(nmds_isd_taxa_k3, Wb=1.2, S=1.6)
#colnames(nmds_isd_taxa_ucie) <- c("scientificName","UCIE")
#write_delim(nmds_isd_taxa_ucie,"results/nmds_isd_taxa_ucie.tsv", delim="\t")

############################# ucie to location data ####################
locations_spatial <- locations_spatial %>%
    left_join(pcoa_isd_sites_ucie,by=c("ENA_RUN"="ENA_RUN"))

locations_spatial$UCIE[is.na(locations_spatial$UCIE)] <- "gray" 
##################################### MAPS #####################################
# Colorblind palette
palette.colors(palette = "Okabe-Ito")
# Crete figures
cols=c("chocolate1","cornflowerblue","darkgoldenrod1", "darkolivegreen4", "darkorchid1", "goldenrod3", "palevioletred3", "peachpuff4", "turquoise","skyblue")

print("printing base maps")

crete_black <- ggplot() +
    geom_sf(crete_shp, mapping=aes()) +
    geom_point(locations_spatial,
            mapping=aes(x=longitude, y=latitude, color=UCIE),
            size=8,
            show.legend=F) +
    geom_jitter(width = 0.25, height = 0.25)+
    scale_color_manual(values=locations_spatial$UCIE, guide="none")+
    coord_sf(crs="wgs84") +
    theme_bw()+
    theme(
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        panel.border = element_blank(),
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        panel.grid.major = element_blank(), #remove major gridlines
        panel.grid.minor = element_blank(), #remove minor gridlines
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent'), #transparent legend panel
        line = element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank(),
        legend.text=element_text(size=8),
        legend.title = element_text(size=8),
        legend.position = "bottom")


ggsave("figures/map_crete_black.png",
       plot=crete_black,
       bg='transparent',
       height = 30,
       width = 60,
       dpi = 300,
       units="cm",
       device="png")


crete_base <- ggplot() +
    geom_sf(crete_shp, mapping=aes()) +
    geom_raster(dem_crete_df, mapping=aes(x=x, y=y, fill=dem_crete))+
    geom_point(locations_spatial,
            mapping=aes(x=longitude, y=latitude, color=UCIE, shape=as.character(route)),
            size=1.7,
            alpha=0.6,
            show.legend=T) +
    geom_jitter(width = 0.25, height = 0.25)+
#    geom_sf(crete_peaks,
#            mapping=aes(),
#            colour=na,
#            size=1,
#            alpha=1,
#            show.legend=f) +
#    geom_label(data = crete_peaks,
#               mapping=aes(x = X, y = Y, label = name),
#               size = 1.8,
#               nudge_x = 0.07,
#               nudge_y=0.07,
#               label.padding = unit(0.1, "lines"))+
    scale_fill_gradientn(guide = guide_colourbar(barwidth = 0.5, barheight = 3.5,
                                  title="elevation",
                                  direction = "vertical",
                                  title.vjust = 0.8),
                        colours = c("snow3","#f0e442","#d55e00","#cc79a7"),
                        breaks = c(100, 800, 1500, 2400),
                        labels = c(100, 800, 1500, 2400))+
    scale_color_manual(values=locations_spatial$UCIE, guide="none")+
    scale_shape_manual(values=c(seq(0,9,1)),name="route")+
    coord_sf(crs="wgs84") +
    theme_bw()+
    theme(axis.title=element_blank(),
          axis.text=element_text(colour="black"),
          legend.text=element_text(size=8),
          legend.title = element_text(size=8),
          legend.position = "bottom",
          legend.box.background = element_blank())


ggsave("figures/map_fig1a.tiff",
       plot=crete_base,
       height = 10,
       width = 20,
       dpi = 300,
       units="cm",
       device="tiff")

ggsave("figures/map_fig1a.png",
       plot=crete_base,
       height = 10,
       width = 20,
       dpi = 300,
       units="cm",
       device="png")

## Crete Corine

crete_corine <- ggplot() +
    geom_sf(crete_shp, mapping=aes()) +
    geom_sf(clc_crete_shp,
            mapping=aes(fill=LABEL1),
            alpha=1,
            colour="transparent",
            show.legend=T) +
    geom_sf(natura_crete_land_sci,
            mapping=aes(color="Natura2000 HSD"),
            linewidth=0.6,
            fill=NA,
            alpha=1,
            show.legend=T) +
#    geom_sf(crete_peaks,
#            mapping=aes(),
#            color = "#D55E00",
#            size=1,
#            alpha=1,
#            show.legend=F) +
#    geom_label(data = crete_peaks, 
#               mapping=aes(x = X, y = Y, label = name),
#               size = 1.5,
#               nudge_x = 0.05,
#               nudge_y=0.05, label.padding = unit(0.1, "lines"))+ 
    scale_fill_manual(values = c("Artificial surfaces"="#000000",
                                 "Agricultural areas"="#E69F00",
                                 "Forest and semi natural areas" = "#009E73",
                                 "Water bodies" = "#0072B2",
                                 "Natura2000 HSD"=NA),
                      guide = "legend") +
    scale_colour_manual(values = c("Natura2000 HSD" = "#56B4E9"),
                        guide = "legend") +
    guides(fill = guide_legend(override.aes = list(color = "transparent", alpha=1) ),
           colour = guide_legend(override.aes = list(alpha=1, fill="transparent") ) )+
    coord_sf(crs="WGS84") +
    theme_bw()+
    theme(axis.title=element_blank(),
          axis.text=element_text(colour="black"),
          legend.title = element_blank(),
          legend.position = "bottom",
          legend.box.background = element_blank(),
          legend.key.size = unit(8, "mm"), 
          legend.text=element_text(size=8))

ggsave("figures/map_fig1b.tiff", 
       plot=crete_corine, 
       height = 10, 
       width = 20,
       dpi = 300, 
       units="cm",
       device="tiff")

ggsave("figures/map_fig1b.png", 
       plot=crete_corine, 
       height = 10, 
       width = 20,
       dpi = 300, 
       units="cm",
       device="png")


fig1 <- ggarrange(crete_base,crete_corine,
          labels = c("A", "B"),
          align = "hv",
          widths = c(1,0.6),
          ncol = 1,
          nrow = 2,
          font.label=list(color="black",size=22),
          legend="bottom") + bgcolor("white")

ggsave("figures/map_fig1.tiff", 
       plot=fig1, 
       height = 30, 
       width = 30,
       dpi = 300, 
       units="cm",
       device="tiff")

ggsave("figures/map_fig1.png", 
       lot=fig1, 
       height = 30, 
       width = 30,
       dpi = 300, 
       units="cm",
       device="png")

ggsave("figures/map_fig1.pdf", 
       plot=fig1, 
       height = 30, 
       width = 30,
       dpi = 300, 
       units="cm",
       device="pdf")

ggsave("figures/map_fig1-small.png", 
       plot=fig1, 
       height = 30, 
       width = 30,
       dpi = 300, 
       units="cm",
       device="png")

#################### geology #######################
cols=c("chocolate1","cornflowerblue","darkgoldenrod1", "darkolivegreen4", "darkorchid1", "goldenrod3", "palevioletred3", "peachpuff4", "turquoise","skyblue")

crete_geology_plot <- ggplot() +
    geom_sf(crete_shp, mapping=aes()) +
    geom_sf(crete_geology,
            mapping=aes(fill=geology_na),
            alpha=1,
            colour="transparent",
            show.legend=T) +
    geom_point(locations_spatial,
            mapping=aes(x=longitude, y=latitude),
            size=1.7,
            color="gray50",
            alpha=0.4,
            show.legend=F) +
    scale_fill_manual(values = c("-"="#000000",
                                 "fo"="chocolate1",
                                 "ft" = "cornflowerblue",
                                 "J-E" = "darkgoldenrod1",
                                 "K-E" = "darkolivegreen4",
                                 "K.k" = "darkorchid1",
                                 "Mk" = "goldenrod3",
                                 "Mm.I" = "palevioletred3",
                                 "o" = "peachpuff4",
                                 "Ph-T" = "turquoise",
                                 "Q.al" = "skyblue",
                                 "T.br" = "#A6CEE3",
                                 "f" = "#1F78B4",
                                 "K.m" = "#B2DF8A"),
                      guide = "legend") +
    guides(fill = guide_legend(override.aes = list(color = "transparent", alpha=1) ),
           colour = "none" )+
    coord_sf(crs="WGS84") +
    theme_bw()+
    theme(axis.title=element_blank(),
          axis.text=element_text(colour="black"),
          legend.title = element_blank(),
          legend.position = "bottom",
          legend.box.background = element_blank(),
          legend.key.size = unit(8, "mm"), 
          legend.text=element_text(size=8))

ggsave("figures/map_crete_geology.tiff", 
       plot=crete_geology_plot, 
       height = 10, 
       width = 20,
       dpi = 300, 
       units="cm",
       device="tiff")

ggsave("figures/map_crete_geology.png", 
       plot=crete_geology_plot, 
       height = 10, 
       width = 20,
       dpi = 300, 
       units="cm",
       device="png")

################################# Taxonomy #####################################
## Phyla distribution, average relative abundance and ubiquity
## Biogeography of soil bacteria and archaea across France
#phyla_samples_summary <- read_delim("results/phyla_samples_summary.tsv",delim="\t")

total_samples <- length(unique(community_matrix_l$ENA_RUN))

phyla_samples_summary <- community_matrix_l %>%
    group_by(ENA_RUN,Phylum) %>%
    summarise(asvs=sum(asvs),
              reads_srs_mean=mean(reads_srs_mean),
              reads_srs_sum=sum(reads_srs_sum), .groups="keep") %>%
    group_by(ENA_RUN) %>%
    mutate(relative_srs=reads_srs_sum/sum(reads_srs_sum)) %>%
    mutate(z_srs=(reads_srs_sum-mean(reads_srs_sum))/sd(reads_srs_sum))
    
phyla_stats <- phyla_samples_summary %>% 
    group_by(Phylum) %>%
    summarise(n_samples=n(),
              total_asvs=sum(asvs),
              total_reads_srs=sum(reads_srs_sum),
              proportion_sample=n_samples/total_samples,
              average_relative=mean(relative_srs)) %>%
    arrange(desc(average_relative))

############################### phyla ratios matter ###########################
n_phyla <- length(unique(phyla_samples_summary$Phylum))
okabe_ito_colors <- palette.colors(palette = "Okabe-Ito")   
fill_colors <- colorRampPalette(okabe_ito_colors)(n_phyla)

phyla_ratios_bar <- ggplot() + 
    geom_col(phyla_samples_summary, mapping=aes(x=ENA_RUN,
                                                y=relative_srs,
                                                fill=Phylum)) +
    scale_fill_manual(values=fill_colors) +
    theme_bw()+
    coord_flip()+
    scale_y_continuous(expand = c(0, 0), name="")+
    theme(axis.text.x = element_text(face="bold",
                                     size = 15),
          axis.ticks.y=element_blank(),axis.title=element_blank(),
          axis.text.y=element_text(size=10, hjust=0, vjust=0)) +
    theme(legend.position="bottom",
            panel.border = element_blank(),
            panel.grid.major = element_blank(), #remove major gridlines
            panel.grid.minor = element_blank())+
    guides(fill=guide_legend(nrow=7,byrow=TRUE))

ggsave("figures/taxonomy_phyla_ratios_samples.png",
       plot=phyla_ratios_bar,
       device="png",
       height = 55,
       width = 38,
       units="cm")

################################# Heatmap ###############################
phyla_samples_w_z <- phyla_samples_summary %>%
    pivot_wider(id_cols=ENA_RUN,
                names_from=Phylum,
                values_from=z_srs,
                values_fill=0) %>%
    as.data.frame() %>% 
    column_to_rownames("ENA_RUN")

phyla_samples_w <- phyla_samples_summary %>%
    pivot_wider(id_cols=ENA_RUN,
                names_from=Phylum,
                values_from=reads_srs_sum,
                values_fill=0) %>%
    as.data.frame() %>%
    column_to_rownames("ENA_RUN")

dcols = vegdist(phyla_samples_w, method="bray")
drows = vegdist(t(phyla_samples_w), method="robust.aitchison")

png("figures/taxonomy_heatmap_phyla_samples.png",
    res=300,
    width=70,
    height=30,
    unit="cm")

pheatmap(t(phyla_samples_w_z),
         clustering_distance_rows = drows,
         clustering_distance_cols = dcols,
         color=colorRampPalette(c("white", "skyblue", "palevioletred3"))(20))

dev.off()

############# Representative phyla of Cretan soils ########################

representative_phyla_rel <- ggplot(phyla_stats,
       aes(x = average_relative, y = reorder(Phylum, average_relative))) +
  geom_point(size = 2) +  # Use a larger dot
  scale_x_log10(labels = label_log()) +
#  scale_x_continuous()+
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey60", linetype = "dashed")
  )

ggsave("figures/taxonomy_representative_phyla_rel.png",
       plot=representative_phyla_rel,
       device="png",
       height = 20,
       width = 23,
       units="cm")

representative_phyla_samples <- ggplot(phyla_stats,
       aes(x = proportion_sample, y = reorder(Phylum, proportion_sample))) +
  geom_point(size = 2) +  # Use a larger dot
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey60", linetype = "dashed")
  )

ggsave("figures/taxonomy_representative_phyla_samples.png",
       plot=representative_phyla_samples,
       device="png",
       height = 20,
       width = 23,
       units="cm")

################################ Phyla and genera  #########################


phyla_genera_bar <- ggplot(genera_phyla_stats, mapping=aes(x=Phylum, y=average_relative,fill=Genus)) +
    geom_col(position="stack") +
    coord_flip()+
    theme_bw()+
    theme(legend.position="none")


ggsave("figures/taxonomy_genera_phyla_stats.png",
       plot=phyla_genera_bar,
       device="png",
       height = 20,
       width = 23,
       units="cm")



########################### Generalists and specialists ##################
## Abundance (y axis) and occupancy (x axis) plot
## similar to Using network analysis to explore co-occurrence patterns in soil microbial communities

################## ASV  #################################
asv_stat_sample <- ggplot() +
    geom_point(asv_metadata,
               mapping=aes(x=n_samples, y=reads_srs_mean, color=classification)) +
    geom_errorbar(asv_metadata,
                  mapping=aes(x=n_samples,
                              y=reads_srs_mean,
                              ymin=reads_srs_mean-reads_srs_sd,
                              ymax=reads_srs_mean+reads_srs_sd,
                              alpha=0.5, colour=classification))+
    scale_x_continuous(breaks=seq(0,150,10), name="Number of samples") +
#    scale_y_continuous(trans='log10', name = "ASVs",
#                     breaks=trans_breaks('log10', function(x) 10^x),
#                     labels=trans_format('log10', math_format(10^.x))) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size=13),
          axis.title.x=element_text(face="bold", size=13),
          axis.title.y=element_text(face="bold", size=13),
          legend.position = c(0.88, 0.8))

ggsave("figures/taxonomy_asv_generalists.png",
       plot=asv_stat_sample,
       device="png",
       height = 20,
       width = 23,
       units="cm")

asv_stat_sample_cla <- ggplot(data= asv_metadata, mapping=aes(x=n_samples, y=reads_srs_mean)) +
    geom_point(asv_metadata,
               mapping=aes(x=n_samples, y=reads_srs_mean, color=classification)) +
    geom_errorbar(asv_metadata,
                  mapping=aes(x=n_samples,
                              y=reads_srs_mean,
                              ymin=reads_srs_mean-reads_srs_sd,
                              ymax=reads_srs_mean+reads_srs_sd,
                              alpha=0.5, colour=classification))+
    scale_x_continuous(breaks=seq(0,150,10), name="Number of samples") +
#    scale_y_continuous(trans='log10', name = "ASVs",
#                     breaks=trans_breaks('log10', function(x) 10^x),
#                     labels=trans_format('log10', math_format(10^.x))) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size=13),
          axis.title.x=element_text(face="bold", size=13),
          axis.title.y=element_text(face="bold", size=13),
          legend.position = "bottom") + 
    facet_wrap(vars(classification))

ggsave("figures/taxonomy_asv_generalists_cla.png",
       plot=asv_stat_sample_cla,
       device="png",
       height = 20,
       width = 23,
       units="cm")

############################## genera ####################

genera_stat_sample <- ggplot() +
    geom_point(genera_phyla_stats,
               mapping=aes(x=proportion_sample,
                           y=reads_srs_mean,
                           color=Phylum, size=average_relative)) +
#    geom_errorbar(genera_phyla_stats,
#                  mapping=aes(x=n_samples,
#                              y=reads_srs_mean,
#                              ymin=reads_srs_mean-reads_srs_sd,
#                              ymax=reads_srs_mean+reads_srs_sd,
#                              alpha=0.5, colour=Phylum))+
#    scale_x_continuous(breaks=seq(0,150,10), name="Number of samples") +
#    scale_y_continuous(trans='log10', name = "ASVs",
#                     breaks=trans_breaks('log10', function(x) 10^x),
#                     labels=trans_format('log10', math_format(10^.x))) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size=13),
          axis.title.x=element_text(face="bold", size=13),
          axis.title.y=element_text(face="bold", size=13),
          legend.position = "bottom")

ggsave("figures/taxonomy_genera_generalists.png",
       plot=genera_stat_sample,
       device="png",
       height = 20,
       width = 35,
       units="cm")

genera_stat_sample_f <- genera_stat_sample + facet_wrap(vars(Phylum))

ggsave("figures/taxonomy_genera_generalists_facet.png",
       plot=genera_stat_sample_f,
       device="png",
       height = 50,
       width = 80,
       units="cm")

########################## specialists ##############################
phyla_specialists <- phyla_stats %>%
    filter(proportion_sample<0.25) 

phyla_specialists_samples <- phyla_genera_samples_summary %>%
    filter(Phylum %in% phyla_specialists$Phylum)

base_crete_map <- ggplot() +
        geom_sf(crete_shp, mapping=aes()) +
        geom_raster(dem_crete_df, mapping=aes(x=x, y=y, fill=dem_crete))+
        scale_fill_gradientn(guide = guide_colourbar(barwidth = 0.5, barheight = 3.5,
                                      title="elevation",
                                      direction = "vertical",
                                      title.vjust = 0.8),
                            colours = c("snow3","#f0e442","#d55e00","#cc79a7"),
                            breaks = c(100, 800, 1500, 2400),
                            labels = c(100, 800, 1500, 2400))+
        coord_sf(crs="wgs84") +
        theme_bw()+
        theme(axis.title=element_blank(),
              axis.text=element_text(colour="black"),
              legend.text=element_text(size=8),
              legend.title = element_text(size=8),
              legend.position = "bottom",
              legend.box.background = element_blank())

for (phy in phyla_stats$Phylum){

    phyla_samples_loc <- phyla_genera_samples_summary %>%
        filter(Phylum==phy) %>%
        left_join(locations_spatial, by=c("ENA_RUN"="ENA_RUN"))
    
    crete_phylum <- base_crete_map +
        geom_point(phyla_samples_loc,
                mapping=aes(x=longitude, y=latitude, color=relative_srs, size=relative_srs),
                alpha=0.8,
                show.legend=T) +
        geom_label(data = phyla_samples_loc,
                   mapping=aes(x = longitude, y = latitude, label = ENA_RUN),
                   size = 1.8,
                   nudge_x = 0.07,
                   nudge_y=0.07,
                   label.padding = unit(0.1, "lines"))+
        scale_color_gradient(low = "skyblue", high = "goldenrod3")+
        ggtitle(phy) 
    
    
    ggsave(paste0("figures/map_phyla_",phy,".png"),
           plot=crete_phylum,
           height = 10,
           width = 20,
           dpi = 300,
           units="cm",
           device="png")
}

######################### ASVs distribution vs samples #####################

asv_sample_dist_t <- asv_metadata %>%
    group_by(n_samples) %>%
    summarise(n_asv=n()) %>% 
    mutate(classification="total")

asv_sample_dist_c <- asv_metadata %>%
    group_by(n_samples, classification) %>%
    summarise(n_asv=n(), .groups="keep")

asv_sample_dist <- rbind(asv_sample_dist_t, asv_sample_dist_c)


### do the generalist and specialist
### mean abundance and n sites

asv_sample_dist_plot <- ggplot() +
    geom_point(asv_sample_dist,
               mapping=aes(x=n_asv, y=n_samples, color=classification)) +
    scale_y_continuous(breaks=seq(0,150,10), name="Number of samples") +
    scale_x_continuous(trans='log10', name = "ASVs",
                     breaks=trans_breaks('log10', function(x) 10^x),
                     labels=trans_format('log10', math_format(10^.x))) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size=13),
          axis.title.x=element_text(face="bold", size=13),
          axis.title.y=element_text(face="bold", size=13),
          legend.position = c(0.88, 0.8))

ggsave("figures/taxonomy_asv_n_samples.png",
       plot=asv_sample_dist_plot,
       device="png",
       height = 20,
       width = 23,
       units="cm")


##################### Sample Diversity boxplots #########################
# Metadata to long format for diversity indices
metadata_diversity <- metadata %>%
    pivot_longer(cols=c(asvs,
                        shannon,
                        S.obs,
                        S.chao1,
                        se.chao1,
                        S.ACE,
                        se.ACE,
                        taxa),
                 names_to="diversity",
                 values_to="value") %>%
    as.data.frame()
# metadata abiotic values
metadata_abiotic <- metadata %>%
    pivot_longer(cols=c(latitude,
                        longitude,
                        elevation,
                        total_nitrogen,
                        water_content,
                        total_organic_carbon,
                        sample_volume_or_weight_for_DNA_extraction,
                        DNA_concentration),
                 names_to="abiotic_metadata",
                 values_to="value") %>%
    as.data.frame()
# metadata bioclim
metadata_bioclim <- metadata %>%
    pivot_longer(cols = grep("bio_*", colnames(.)),
                 names_to="bioclim_metadata",
                 values_to="value") %>%
    as.data.frame()

#################### Sample diversity gradients ######################
## similar to Structure and function of the global topsoil microbiome but
## without the statistic test.

# Numerical variables to plot against diversity indices

vars <- c("latitude",
          "longitude",
          "elevation",
          "total_nitrogen",
          "water_content",
          "total_organic_carbon",
          "sample_volume_or_weight_for_DNA_extraction",
          "DNA_concentration",
          "route")

for (var in vars){

    # bioclimatic variables boxplots
    gradient_scatterplot(metadata_bioclim, var, "value", "bioclim_metadata")
}

# add the bioclim variables to the variables
bioclim <- grep("bio.*", colnames(metadata_diversity),value=T)
vars <- c(vars, bioclim)
for (var in vars){
    
    gradient_scatterplot(metadata_diversity, var, "value", "diversity")
}

# abiotic pairs scatterplots
abiotic <- c("total_nitrogen",
             "water_content",
             "total_organic_carbon",
             "sample_volume_or_weight_for_DNA_extraction",
             "DNA_concentration")

abiotic_comb <- as.data.frame(t(combn(abiotic,2)))

gradient_scatterplot(metadata_diversity,"total_nitrogen", "total_organic_carbon", "elevation_bin") 
gradient_scatterplot(metadata_diversity,"total_nitrogen", "total_organic_carbon", "LABEL1") 
gradient_scatterplot(metadata_diversity,"total_nitrogen", "total_organic_carbon", "LABEL2") 

# Categorical variables to plot against diversity indices
cats <- c("vegetation_zone",
          "LABEL1",
          "LABEL2",
          "LABEL3",
          "elevation_bin",
          "location",
          "protection_status",
          "geology_na")

metadata_diversity <- metadata_diversity %>% filter(!is.na(geology_na))
for (cat in cats){
    # diversity indices boxplots
    diversity_boxplot(metadata_diversity, cat, "value", "diversity")
    # abiotic metadata boxplots
    diversity_boxplot(metadata_abiotic, cat, "value", "abiotic_metadata") 
    # bioclimatic variables boxplots
    diversity_boxplot(metadata_bioclim, cat, "value", "bioclim_metadata")
}

######################### Ordination ############################
######### plots sites ###########


nmds_isd_sites <- read_delim("results/nmds_isd_sites.tsv", delim="\t")
umap_isd_sites <- read_delim("results/umap_samples_2.tsv", delim="\t")
pcoa_isd_sites <- read_delim("results/ordination_pcoa_bray_sites.tsv", delim="\t")
#umap_isd_sites_k1 <- read_delim("results/umap_samples_1.tsv", delim="\t")
#colnames(umap_isd_sites_k1) <- c("id", "UCIE")

ordination_sites <- nmds_isd_sites %>%
    left_join(umap_isd_sites, by=c("ENA_RUN"="id")) %>%
    left_join(pcoa_isd_sites_ucie) %>%
    left_join(pcoa_isd_sites) %>% 
    left_join(metadata)
#    left_join(umap_isd_sites_k1 ,by=c("ENA_RUN"="id"))
ordination_sites$elevation_bin <- factor(ordination_sites$elevation_bin,
                        levels=unique(ordination_sites$elevation_bin)[order(sort(unique(ordination_sites$elevation_bin)))])

for (i in cats){

    if (is.character(ordination_sites[[i]]) & i!="ENA_RUN"){
        print(i)
        ordination_sites_plot(ordination_sites, i,"NMDS1","NMDS2", "nmds","UCIE")
        ordination_sites_plot(ordination_sites, i,"UMAP1","UMAP2", "umap","UCIE")
        ordination_sites_plot(ordination_sites, i,"Axis.1","Axis.2", "pcoa","UCIE")

    }else{
        next
    }
}
# ordination of sites and their 2 locations. Are there differences? 
ordination_sites_plot(ordination_sites,"location","NMDS1","NMDS2","nmds_site_loc","sites")
### 3d for inspection of the dimention reduction
# 
#library(plotly)
#plot_ly(x=umap_isd_sites_k3$UMAP1, y=umap_isd_sites_k3$UMAP2, z=umap_isd_sites_k3$UMAP3, type="scatter3d", mode="markers")

########################## Ordination and Boxplots ###########################

boxplot_single(ordination_sites, "UMAP1", "LABEL2", "elevation_bin")
boxplot_single(ordination_sites, "elevation_bin", "UMAP2", "elevation_bin")
ordination_sites_plot(ordination_sites, "elevation_bin","UMAP1","UMAP2", "umap","elevation_bin")

######### plots genera ###########

nmds_isd_taxa <- read_delim("results/nmds_isd_taxa.tsv", delim="\t")
nmds_isd_taxa_ucie <- read_delim("results/nmds_isd_taxa_ucie.tsv", delim="\t")
umap_isd_genera <- read_delim("results/umap_genera_2.tsv", delim="\t")
nmds_isd_taxa <- nmds_isd_taxa %>% left_join(nmds_isd_taxa_ucie) %>% left_join(umap_isd_genera, by=c("scientificName"="id"))

nmds_genera_plot <- ggplot() +
    geom_point(data=nmds_isd_taxa,
               mapping=aes(x=NMDS1, y=NMDS2,color=UCIE),show.legend = F) +
    scale_color_manual(values=nmds_isd_taxa$UCIE)+
    coord_equal() +
    theme_bw()

ggsave("figures/ordination_nmds_genera_plot.png",
       plot=nmds_genera_plot,
       device="png",
       height = 20,
       width = 23,
       units="cm")

nmds_genera_plot_f <- nmds_genera_plot + facet_wrap(vars(Phylum),
                                                    scales="fixed")

ggsave("figures/ordination_nmds_genera_plot_f.png",
       plot=nmds_genera_plot_f,
       device="png",
       height = 50,
       width = 53,
       units="cm")
#### umap
umap_genera_plot <- ggplot() +
    geom_point(data=nmds_isd_taxa,
               mapping=aes(x=UMAP1, y=UMAP2,color=UCIE),show.legend = F) +
    scale_color_manual(values=nmds_isd_taxa$UCIE)+
    coord_equal() +
    theme_bw()

ggsave("figures/ordination_umap_genera_plot.png",
       plot=umap_genera_plot,
       device="png",
       height = 20,
       width = 23,
       units="cm")

umap_genera_plot_f <- umap_genera_plot + facet_wrap(vars(Phylum),
                                                    scales="fixed")

ggsave("figures/ordination_umap_genera_plot_f.png",
       plot=umap_genera_plot_f,
       device="png",
       height = 50,
       width = 53,
       units="cm")



############################## Functional profiles ############################
### function
clr <- function(x, na.rm = FALSE) {
    log1p(x = x/(exp(x = sum(log1p(x = x[x > 0]), na.rm = TRUE)/length(x = x))))
}


scale2 <- function(x, na.rm = FALSE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)
### load data
faprotax_samples <- read_delim("results/faprotax_functional_table.tsv", delim="\t")
faprotax_genera <- read_delim("results/faprotax_genera_functional_table.tsv", delim="\t")

# filter empty
faprotax_genera <- faprotax_genera[rowSums(faprotax_genera[,-1])!=0,]
##### bar plot

################################# Heatmap ###############################

#faprotax_community_matrix <- crete_biodiversity %>%
#    filter(classification %in% c("Genus","Species"), !is.na(srs_abundance)) %>%
#    pivot_wider(id_cols=c(asv_id,taxonomy), names_from=ENA_RUN, values_from=srs_abundance, values_fill=0)

#write_delim(faprotax_community_matrix,"results/faprotax_community_matrix.tsv",delim="\t")

faprotax_genera_w <- faprotax_genera %>%
    filter(!(group %in% c("aerobic_chemoheterotrophy", "chemoheterotrophy"))) %>%
    as.data.frame() %>% 
    column_to_rownames("group") %>% 
    mutate_all(clr) %>%
    mutate_all(sqrt) 

#drows = vegdist(phyla_samples_w, method="bray")
dcols = vegdist(t(faprotax_genera[,-1]), method="bray")

png("figures/functions_faprotax_genera.png",
    res=300,
    width=60,
    height=40,
    unit="cm")

pheatmap(faprotax_genera_w,
         clustering_distance_cols = dcols,
         color=colorRampPalette(c("white",
                                  "skyblue",
                                  "cornflowerblue",
                                  "darkolivegreen4",
                                  "darkgoldenrod1",
                                  "palevioletred3",
                                  "darkorchid1"))(50))

dev.off()

# heatmap with geology order
my_sample_col <- metadata %>%
    filter(ENA_RUN %in% colnames(faprotax_genera_w)) %>%
    dplyr::select(geology_na,ENA_RUN) %>%
    arrange(geology_na) %>%
    column_to_rownames("ENA_RUN") %>%
    as.data.frame()


png("figures/functions_faprotax_genera_geology.png",
    res=300,
    width=60,
    height=40,
    unit="cm")

pheatmap(faprotax_genera_w[,rownames(my_sample_col)],
         annotation_col = my_sample_col, 
         cluster_cols = FALSE,
        # clustering_distance_cols = dcols,
         color=colorRampPalette(c("white",
                                  "skyblue",
                                  "cornflowerblue",
                                  "darkolivegreen4",
                                  "darkgoldenrod1",
                                  "palevioletred3",
                                  "darkorchid1"))(50))

dev.off()

# heatmap with label2 order
my_sample_col <- metadata %>%
    filter(ENA_RUN %in% colnames(faprotax_genera_w)) %>%
    dplyr::select(LABEL2,ENA_RUN) %>%
    arrange(LABEL2) %>%
    column_to_rownames("ENA_RUN") %>%
    as.data.frame()

png("figures/functions_faprotax_genera_LABEL2.png",
    res=300,
    width=70,
    height=40,
    unit="cm")

pheatmap(faprotax_genera_w[,rownames(my_sample_col)],
         annotation_col = my_sample_col, 
         cluster_cols = FALSE,
         color=colorRampPalette(c("white",
                                  "skyblue",
                                  "cornflowerblue",
                                  "darkolivegreen4",
                                  "darkgoldenrod1",
                                  "palevioletred3",
                                  "darkorchid1"))(50))

dev.off()
###### 

faprotax_genera_l <- faprotax_genera %>%
    pivot_longer(-group,names_to="ENA_RUN", values_to="value" ) %>%
    filter(value!=0) %>%
    mutate(value_clr = clr(value))

plant_pathogen <- faprotax_genera_l %>%
    filter(group=="plant_pathogen") %>%
    arrange(desc(value))

human_pathogen <- faprotax_genera_l %>%
    filter(group=="human_pathogens_all") %>%
    arrange(desc(value))

faprotax_bar <- ggplot() + 
    geom_col(faprotax_genera_l, mapping=aes(x=ENA_RUN, y=value ,fill=group)) +
    theme_bw()+
    theme(axis.text.x = element_text(face="bold",
                                     size = 10,
                                     angle = 90,
                                     vjust = 1,
                                     hjust=1)) +
    theme(legend.position="top")


ggsave("figures/faprotax_genera_bar.png",
       plot=faprotax_bar,
       device="png",
       height = 50,
       width = 53,
       units="cm")


faprotax_functions <- unique(faprotax_genera_l$group)

for (fun in faprotax_functions){

    faprotax_genera_loc <- faprotax_genera_l %>%
        filter(group==fun) %>%
        left_join(locations_spatial, by=c("ENA_RUN"="ENA_RUN"))
    
    quantile_4 <- quantile(faprotax_genera_loc$value,na.rm = T,probs = c(0.75)) 
    
    faprotax_genera_l_4<- faprotax_genera_loc %>% filter(value > quantile_4)
    
    
    crete_function <- ggplot() +
        geom_sf(crete_shp, mapping=aes()) +
        geom_raster(dem_crete_df, mapping=aes(x=x, y=y, fill=dem_crete))+
        geom_point(faprotax_genera_l_4,
                mapping=aes(x=longitude, y=latitude, color=value, size=value),
                alpha=0.8,
                show.legend=T) +
        geom_label(data = faprotax_genera_l_4,
                   mapping=aes(x = longitude, y = latitude, label = ENA_RUN),
                   size = 1.8,
                   nudge_x = 0.07,
                   nudge_y=0.07,
                   label.padding = unit(0.1, "lines"))+
        scale_fill_gradientn(guide = guide_colourbar(barwidth = 0.5, barheight = 3.5,
                                      title="elevation",
                                      direction = "vertical",
                                      title.vjust = 0.8),
                            colours = c("snow3","#f0e442","#d55e00","#cc79a7"),
                            breaks = c(100, 800, 1500, 2400),
                            labels = c(100, 800, 1500, 2400))+
        scale_color_gradient(low = "skyblue", high = "goldenrod3")+
        ggtitle(fun) +
        coord_sf(crs="wgs84") +
        theme_bw()+
        theme(axis.title=element_blank(),
              axis.text=element_text(colour="black"),
              legend.text=element_text(size=8),
              legend.title = element_text(size=8),
              legend.position = "bottom",
              legend.box.background = element_blank())
    
    
    ggsave(paste0("figures/map_faprotax_",fun,".png"),
           plot=crete_function,
           height = 10,
           width = 20,
           dpi = 300,
           units="cm",
           device="png")
}

#phyla_samples_w <- phyla_samples_summary %>%
#    pivot_wider(id_cols=ENA_RUN,
#                names_from=Phylum,
#                values_from=reads_srs_sum,
#                values_fill=0) %>%
#    as.data.frame() %>%
#    column_to_rownames("ENA_RUN")


######################## Community ################################



print("all figures were generated")
