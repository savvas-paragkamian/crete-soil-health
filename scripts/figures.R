#!/usr/bin/Rscript

# load packages and functions
library(tidyverse)
library(ggnewscale)
library(ggpubr)
library(sf)
library(jpeg)
library(raster)

# load data
locations_spatial <- sf::st_read("../results/locations_spatial/locations_spatial.shp")
crete_shp <- sf::st_read("../data/crete/crete.shp")
crete_peaks <- read_delim("../data/crete_mountain_peaks.csv", delim=";", col_names=T) %>%
    st_as_sf(coords=c("X", "Y"),
             remove=F,
             crs="WGS84")
endemic_species <- read_delim("../results/endemic_species_assessment.tsv", delim="\t")
clc_crete_shp <- st_read("../data/clc_crete_shp/clc_crete_shp.shp")
natura_crete <- sf::st_read("../data/natura2000/natura2000_crete.shp")
wdpa_crete <- sf::st_read("../data/wdpa_crete/wdpa_crete.shp")
natura_crete_land <- st_intersection(natura_crete, crete_shp)

# raster DEM hangling
dem_crete <- raster("../data/dem_crete/dem_crete.tif")
dem_crete_pixel <- as(dem_crete, "SpatialPixelsDataFrame")
dem_crete_df <- as.data.frame(dem_crete_pixel) %>% filter(dem_crete>0)


# split the SPA SCI

natura_crete_land_sci <- natura_crete_land %>% filter(SITETYPE=="B")
## Hotspots and threatspots
endemic_hotspots <- st_read("../results/endemic_hotspots/endemic_hotspots.shp")
threatspots <- st_read("../results/threatspots/threatspots.shp")
threatspots_lt <- threatspots %>%
    filter(pc_thrt>= quantile(pc_thrt,0.90))

locations_inland <- st_join(locations_shp, crete_shp, left=F)


# Colorblind palette
palette.colors(palette = "Okabe-Ito")
# Crete figures


fig1a <- readJPEG("../figures/Fig1a.jpg")

g_fig1a <- ggplot() +
    background_image(fig1a)+
    theme(plot.margin = margin(t=0.5, l=0.7, r=0.7, b=0.5, unit = "cm"))

fig1b <- readJPEG("../figures/Fig1b.jpg")

g_fig1b <- ggplot() +
    background_image(fig1b)+
    theme(plot.margin = margin(t=0.7, l=0.7, r=0.5, b=0.5, unit = "cm"))

fig1ab <- readJPEG("../figures/fig1ab.jpg")

g_fig1ab <- ggplot() +
    background_image(fig1ab)+
    theme(plot.margin = margin(t=0.5, l=0.5, r=0.5, b=0, unit = "cm"))
## Fig1c


crete_base <- ggplot() +
    geom_sf(crete_shp, mapping=aes()) +
    geom_raster(dem_crete_df, mapping=aes(x=x, y=y, fill=dem_crete))+
    scale_fill_gradientn(guide = "colourbar",colours = c("snow3","#f0e442","#d55e00","#cc79a7"),
                        breaks = c(100, 800, 1500, 2400),
                        labels = c(100, 800, 1500, 2400))+
    geom_sf(natura_crete_land_sci,
            mapping=aes(colour="natura2000 sac"),
            linewidth=0.4,
            alpha=1,
            fill=na,
            show.legend=t) +
    scale_colour_manual(values = c("natura2000 sac" = "#56b4e9"),
                        guide = guide_legend(override.aes = list(linetype="solid",shape = na)),
                        name="")+
    new_scale_color()+
    geom_point(locations_source,
            mapping=aes(x=logd, y=latd, color=source, shape=source),
            size=1.8,
            alpha=0.8,
            show.legend=t) +
#    geom_sf(crete_peaks,
#            mapping=aes(),
#            colour=na,
#            size=1,
#            alpha=1,
#            show.legend=f) +
    geom_label(data = crete_peaks,
               mapping=aes(x = x, y = y, label = name),
               size = 1.5,
               nudge_x = 0.05,
               nudge_y=0.05,
               label.padding = unit(0.1, "lines"))+
    scale_colour_manual(values = c("nhmc" = "#009e73",
                                   "bibliography" = "#999999"),
                        guide = guide_legend(override.aes = list(size= c(3,3),linetype = c("blank", "blank"))),
                        name = "sampling") +
    scale_shape_manual(values = c("nhmc" = 17,
                                   "bibliography" = 4),
                        name = "sampling") +
    guides(fill = guide_colourbar(ticks = false,
                                  label = true,
                                  title="elevation",
                                  title.vjust = 0.8),
           colour = guide_legend())+
    coord_sf(crs="wgs84") +
    theme_bw()+
    theme(axis.title=element_blank(),
          axis.text=element_text(colour="black"),
          legend.title = element_text(size=8),
          legend.position = "bottom",
          legend.box.background = element_blank())


ggsave("../figures/fig1c.tiff",
       plot=crete_base,
       height = 10,
       width = 20,
       dpi = 600,
       units="cm",
       device="tiff")

ggsave("../figures/fig1c.png",
       plot=crete_base,
       height = 10,
       width = 20,
       dpi = 600,
       units="cm",
       device="png")

