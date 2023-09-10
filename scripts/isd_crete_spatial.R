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
library(ggnewscale)
library(ggpubr)
library(jpeg)
library(raster)
library(scales)


# Metadata

metadata_long <- read_delim("ena_metadata/ena_isd_2016_attributes.tsv", delim="\t") %>%
    mutate(VALUE=gsub("\\r(?!\\n)","", VALUE, perl=T))
# metadata to wide format

metadata_wide <- metadata_long %>% 
    dplyr::select(-c(UNITS)) %>%
    mutate(TAG=gsub(" ","_", TAG, perl=T)) %>%
    pivot_wider(names_from=TAG, 
                values_from=VALUE)

metadata_wide$ENA_RUN <- metadata_wide$`ENA-RUN`
metadata_wide$total_nitrogen <- as.numeric(metadata_wide$total_nitrogen)
metadata_wide$water_content <- as.numeric(metadata_wide$water_content)
metadata_wide$total_organic_carbon <- as.numeric(metadata_wide$total_organic_carbon)
metadata_wide$sample_volume_or_weight_for_DNA_extraction <- as.numeric(metadata_wide$sample_volume_or_weight_for_DNA_extraction)
metadata_wide$DNA_concentration <- as.numeric(metadata_wide$DNA_concentration)
metadata_wide$latitude <- as.numeric(metadata_wide$`geographic_location_(latitude)`)
metadata_wide$longitude <- as.numeric(metadata_wide$`geographic_location_(longitude)`)
metadata_wide$elevation <- as.numeric(metadata_wide$`geographic_location_(elevation)`)
metadata_wide$amount_or_size_of_sample_collected <- as.numeric(metadata_wide$amount_or_size_of_sample_collected)

metadata_spatial <- metadata_wide %>%
    dplyr::select(ENA_RUN, source_material_identifiers, total_nitrogen, water_content,total_organic_carbon,sample_volume_or_weight_for_DNA_extraction,DNA_concentration,latitude, longitude,elevation, amount_or_size_of_sample_collected, vegetation_zone) %>%
    arrange(ENA_RUN) %>%
    mutate(route = sub("isd_(.*.)_site.*" ,"\\1" , source_material_identifiers)) %>%
    st_as_sf(coords=c("longitude", "latitude"),
             remove=F,
             crs="WGS84")

# spatial data
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

# Enrichment

metadata_spatial <- st_join(metadata_spatial, clc_crete_shp, left=T) #%>%
metadata_spatial <- st_join(metadata_spatial, natura_crete, left=T)
metadata_spatial$dem <- raster::extract(dem_crete, metadata_spatial, cellnumbers=F)

metadata <- metadata_spatial %>% st_drop_geometry()

write_delim(metadata,"results/metadata_spatial.tsv", delim="\t")


# Colorblind palette
palette.colors(palette = "Okabe-Ito")
# Crete figures
cols=c("chocolate1","cornflowerblue","darkgoldenrod1", "darkolivegreen4", "darkorchid1", "goldenrod3", "palevioletred3", "peachpuff4", "turquoise","skyblue")

crete_base <- ggplot() +
    geom_sf(crete_shp, mapping=aes()) +
    geom_raster(dem_crete_df, mapping=aes(x=x, y=y, fill=dem_crete))+
    geom_point(metadata_spatial,
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

