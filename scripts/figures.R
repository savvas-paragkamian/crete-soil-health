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

################################# Taxonomy #####################################
## Phyla distribution, average relative abundance and ubiquity
## Biogeography of soil bacteria and archaea across France

phyla_samples_summary <- read_delim("results/phyla_samples_summary.tsv",delim="\t")

############# Representative phyla of Cretan soils ########################

phyla_stats <- read_delim("results/phyla_stats.tsv",delim="\t")

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

ggsave("figures/representative_phyla_rel.png",
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

ggsave("figures/representative_phyla_samples.png",
       plot=representative_phyla_samples,
       device="png",
       height = 20,
       width = 23,
       units="cm")

################## Phyla and genera  #########################

genera_phyla_stats <- read_delim("results/genera_phyla_stats.tsv",delim="\t")

phyla_genera_bar <- ggplot(genera_phyla_stats, mapping=aes(x=Phylum, y=average_relative,fill=Genus)) +
    geom_col(position="stack", stat="identity") +
    coord_flip()+
    theme_bw()+
    theme(legend.position="none")


ggsave("figures/genera_phyla_stats.png",
       plot=phyla_genera_bar,
       device="png",
       height = 20,
       width = 23,
       units="cm")

#################### testing ##############
phyla_samples_m <- phyla_samples_summary %>%
    dplyr::select(ENA_RUN,relative_srs,Phylum) %>%
    pivot_wider(names_from=Phylum,
                values_from=relative_srs,
                values_fill=0) %>%
    as.data.frame()

rownames(phyla_samples_m) <- phyla_samples_m[,1]
phyla_samples_m <- phyla_samples_m[,-1]

library(vegan)
dist_long <- function(x,method){
    method <- method
    df <- as.data.frame(as.matrix(x)) %>%
    rownames_to_column() %>%
    pivot_longer(-rowname,
                 values_to=method,
                 names_to="colname")

    return(df)
}
bray <- vegdist(phyla_samples_m,
                method="bray")
plot(hclust(bray))

bray_phy <- vegdist(t(phyla_samples_m),method="bray")
plot(hclust(bray_phy))

bray_l <- dist_long(bray, "bray")
biodiversity_srs_t <- phyla_samples_m

z <- betadiver(biodiversity_srs_t, "z")
mod <- with(metadata, betadisper(z, LABEL1))
#sac <- specaccum(biodiversity_srs_t)

# Ordination
nmds <- vegan::metaMDS(biodiversity_srs_t,
                       k=2,
                       distance = "bray",
                       trymax=100)
metadata <- metadata %>%
    filter(ENA_RUN %in% rownames(phyla_samples_m))

stressplot(nmds)
ordiplot(nmds,display="sites", cex=1.25)
ordisurf(nmds,metadata$dem,main="",col="forestgreen")
ordihull(nmds,display="sites",label=T,  groups=metadata$LABEL1, cex=1.25)
ordihull(nmds,display="sites",label=T,  groups=metadata$elevation, cex=1.25)


## table for stats
phyla_dist_samples_d <- phyla_dist_samples %>% 
    as.data.frame()

rownames(phyla_dist_samples_d) <- phyla_dist_samples_d$Phylum
phyla_dist_samples_d <- phyla_dist_samples_d[,-1]

plot(hclust(dist(phyla_dist_samples_d),"average" ))

## create bar plots for each sample at family level, class level, Phylum etc
## 

########################### ASVs generalists and specialists ##################
## Abundance (y axis) and occupancy (x axis) plot
## similar to Using network analysis to explore co-occurrence patterns in soil microbial communities

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

ggsave("figures/fig_asv_generalists.png",
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

ggsave("figures/fig_asv_generalists_cla.png",
       plot=asv_stat_sample_cla,
       device="png",
       height = 20,
       width = 23,
       units="cm")
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

ggsave("figures/fig_asv_n_samples.png",
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
## function
diversity_boxplot <- function(dataset, x_axis, y_axis, grouping_var){
    plotname <- paste0("figures/",
                       grouping_var,
                       "_",
                       x_axis,
                       "_boxplot.png")
    x_lab <- x_axis
    y_lab <- y_axis
    dataset$x_axis <- dataset[,x_axis]
    dataset$y_axis <- dataset[,y_axis]
    dataset$grouping_var <- dataset[,grouping_var]

    box_diversity <- ggplot(dataset, mapping=aes(x=x_axis, y=y_axis)) +
        geom_boxplot()+
        geom_jitter(width = 0.2)+
        xlab(x_lab)+
        ylab(y_lab) +
        theme_bw()+
        theme(axis.text.x = element_text(face="bold",
                                         size = 10,
                                         angle = 45,
                                         vjust = 1,
                                         hjust=1)) +
        facet_wrap(vars(diversity), scales = "free")
    
    ggsave(plotname, 
           plot=box_diversity, 
           device="png", 
           height = 45, 
           width = 30, 
           units="cm")
}

# Numerical variables to plot against diversity indices
vars <- c( "latitude",
          "longitude",
          "elevation",
          "total_nitrogen",
          "water_content",
          "total_organic_carbon",
          "sample_volume_or_weight_for_DNA_extraction",
          "DNA_concentration",
          "route")

for (var in vars){
    
    gradient_scatterplot(metadata_diversity, var, "value", "diversity")
}
#################### Sample diversity gradients ######################
## similar to Structure and function of the global topsoil microbiome but
## without the statistic test.

gradient_scatterplot <- function(dataset, x_axis, y_axis, grouping_var){
    
    ## the dataset must be a dataframe, not a tibble, the colnames
    ## must characters

    ## keep the character names of column names to pass to plot
    ##
    plotname <- paste0("figures/",
                       grouping_var,
                       "_",
                       x_axis,
                       "_gradient.png")
    x_lab <- x_axis
    y_lab <- y_axis
    dataset$x_axis <- dataset[,x_axis]
    dataset$y_axis <- dataset[,y_axis]
    dataset$grouping_var <- dataset[,grouping_var]

    gradient <- ggplot(dataset, mapping=aes(x=x_axis, y=y_axis)) +
        geom_point()+
        xlab(x_lab)+
        ylab(y_lab) +
        theme_bw() +
        theme(panel.grid = element_blank(),
              axis.text = element_text(size=13),
              axis.title.x=element_text(face="bold", size=13),
              axis.title.y=element_text(face="bold", size=13),
              legend.position = c(0.88, 0.8)) +
        facet_wrap(vars(diversity), scales = "free")
   
    ggsave(plotname,
           plot=gradient,
           device="png",
           height = 20,
           width = 23,
           units="cm")
}

# Categorical variables to plot against diversity indices
cats <- c("vegetation_zone", "LABEL1","LABEL2","LABEL3","elevation_bin")

for (cat in cats){
    diversity_boxplot(metadata_diversity, cat, "value", "diversity")
}


######################## Community ################################




