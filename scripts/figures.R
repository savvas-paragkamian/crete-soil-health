#!/usr/bin/Rscript
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
library(tidyverse)
library(ggnewscale)
library(ggpubr)
library(sf)
library(jpeg)
library(raster)
library(scales)

################################## Load data ##################################
## biodiversity
crete_biodiversity <- read_delim("results/crete_biodiversity_asv.tsv",delim="\t")
asv_metadata <- read_delim("results/asv_metadata.tsv",delim="\t")

tax_tab <- readRDS("results/tax_tab.RDS")

# Metadata

metadata <- read_delim("results/sample_metadata.tsv", delim="\t")

## spatial
locations_spatial <- read_delim("spatial_data/ISD_sites_coordinates.tsv", delim="\t",col_names=T) %>%
    st_as_sf(coords=c("longitude", "latitude"),
             remove=F,
             crs="WGS84")
crete_shp <- sf::st_read("spatial_data/crete/crete.shp")
crete_peaks <- read_delim("spatial_data/crete_mountain_peaks.csv", delim=";", col_names=T) %>%
    st_as_sf(coords=c("X", "Y"),
             remove=F,
             crs="WGS84")

clc_crete_shp <- st_read("spatial_data/clc_crete_shp/clc_crete_shp.shp")
natura_crete <- sf::st_read("spatial_data/natura2000/natura2000_crete.shp")
wdpa_crete <- sf::st_read("spatial_data/wdpa_crete/wdpa_crete.shp")
natura_crete_land <- st_intersection(natura_crete, crete_shp)
natura_crete_land_sci <- natura_crete_land %>% filter(SITETYPE=="B")
# raster DEM hangling
dem_crete <- raster("spatial_data/dem_crete/dem_crete.tif")
dem_crete_pixel <- as(dem_crete, "SpatialPixelsDataFrame")
dem_crete_df <- as.data.frame(dem_crete_pixel) %>% filter(dem_crete>0)

##################################### MAPS #####################################
# Colorblind palette
palette.colors(palette = "Okabe-Ito")
# Crete figures
cols=c("chocolate1","cornflowerblue","darkgoldenrod1", "darkolivegreen4", "darkorchid1", "goldenrod3", "palevioletred3", "peachpuff4", "turquoise","skyblue")

crete_base <- ggplot() +
    geom_sf(crete_shp, mapping=aes()) +
    geom_raster(dem_crete_df, mapping=aes(x=x, y=y, fill=dem_crete))+
    geom_point(locations_spatial,
            mapping=aes(x=longitude, y=latitude, color=route),
            size=1,
            alpha=0.8,
            show.legend=T) +
#    geom_sf(crete_peaks,
#            mapping=aes(),
#            colour=na,
#            size=1,
#            alpha=1,
#            show.legend=f) +
    geom_label(data = crete_peaks,
               mapping=aes(x = X, y = Y, label = name),
               size = 1.8,
               nudge_x = 0.05,
               nudge_y=0.05,
               label.padding = unit(0.1, "lines"))+
    scale_fill_gradientn(guide = guide_colourbar(barwidth = 0.5, barheight = 3.5,
                                  title="elevation",
                                  direction = "vertical",
                                  title.vjust = 0.8),
                        colours = c("snow3","#f0e442","#d55e00","#cc79a7"),
                        breaks = c(100, 800, 1500, 2400),
                        labels = c(100, 800, 1500, 2400))+
    scale_color_manual(values=cols, guide=guide_legend(nrow=4))+
    coord_sf(crs="wgs84") +
    theme_bw()+
    theme(axis.title=element_blank(),
          axis.text=element_text(colour="black"),
          legend.text=element_text(size=8),
          legend.title = element_text(size=8),
          legend.position = "bottom",
          legend.box.background = element_blank())


ggsave("figures/fig1a.tiff",
       plot=crete_base,
       height = 10,
       width = 20,
       dpi = 300,
       units="cm",
       device="tiff")

ggsave("figures/fig1a.png",
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

ggsave("figures/Fig1b.tiff", 
       plot=crete_corine, 
       height = 10, 
       width = 20,
       dpi = 300, 
       units="cm",
       device="tiff")

ggsave("figures/Fig1b.png", 
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

ggsave("figures/Fig1.tiff", 
       plot=fig1, 
       height = 30, 
       width = 30,
       dpi = 300, 
       units="cm",
       device="tiff")

ggsave("figures/Fig1.png", 
       plot=fig1, 
       height = 30, 
       width = 30,
       dpi = 300, 
       units="cm",
       device="png")

ggsave("figures/Fig1.pdf", 
       plot=fig1, 
       height = 30, 
       width = 30,
       dpi = 300, 
       units="cm",
       device="pdf")

ggsave("figures/Fig1-small.png", 
       plot=fig1, 
       height = 30, 
       width = 30,
       dpi = 300, 
       units="cm",
       device="png")

######################### Biodiversity statistics ###############################

## ASVs

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
          legend.position = c(0.88, 0.1))

ggsave("figures/fig_asv_n_samples.png",
       plot=asv_sample_dist_plot,
       device="png",
       height = 20,
       width = 23,
       units="cm")

################################# Taxonomy #####################################
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
