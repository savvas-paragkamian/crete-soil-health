#!/usr/bin/env Rscript

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
library(terra)
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
wdpa_crete_wildlife <- wdpa_crete %>% filter(DESIG_ENG=="Wildlife Refugee") %>%
    mutate(DESIG_ENG = gsub("Wildlife Refugee", "Wildlife Refuge", DESIG_ENG))

natura_crete_land <- st_intersection(natura_crete, crete_shp)
natura_crete_land_sci <- natura_crete_land %>% filter(SITETYPE=="B")
# raster DEM hangling
dem_crete <- raster("spatial_data/dem_crete/dem_crete.tif")
dem_crete_pixel <- as(dem_crete, "SpatialPixelsDataFrame")
dem_crete_df <- as.data.frame(dem_crete_pixel) %>% filter(dem_crete>0)


## Global Aridity Index and Potential Evapotranspiration Database
library(rSDM)

aridity_crete <- raster("spatial_data/crete_aridity_index.tif")
# some points have 0 value
aridity_crete[aridity_crete[] == 0 ] = NA
metadata_aridity <- metadata_spatial %>% dplyr::select(ENA_RUN) 
locations_s <- locs2sp(do.call(rbind, st_geometry(metadata_aridity)))
check.coords <- points2nearestcell(locations_s, aridity_crete)

# merge the points move to closest valued cell
df_corrected_coords <- as.data.frame(check.coords)
sf_aridity <- cbind(st_drop_geometry(metadata_aridity),df_corrected_coords) %>% 
    st_as_sf(coords=c("longitude", "latitude"),
             remove=T,
             crs="WGS84")
# assign the variable to the initial file. The order of the rows is kept the same
metadata_aridity$aridity <- raster::extract(aridity_crete,sf_aridity, cellnumbers=F)  

# tranform the values to the original form as instructed by the manual. 
metadata_aridity$aridity <- metadata_aridity$aridity*0.0001
metadata_aridity$aridity_class <- cut(metadata_aridity$aridity,
                                      breaks=c(0,0.03,0.2,0.5, 0.65,0.9),
                                      labels=c("Hyper Arid", "Arid", "Semi-Arid", "Dry sub-humid", "Humid"))

## Desertification risk
esa3rdp <- rast("spatial_data/ESDAC_CATENA_Desertification2018_RasterFiles/esa3rdp/w001001.adf")
esa3rdp_wgs <- terra::project(esa3rdp, crs(metadata_aridity))

esa3rdp_crete <- crop(esa3rdp_wgs, crete_shp)
esa3rdp_attr <- data.frame(ESA=unique(esa3rdp_crete),
                           ESA_id=seq(1,nrow(unique(esa3rdp_crete)),1))
write_delim(esa3rdp_attr, "spatial_data/crete_desertification_risk/esa3rdp_crete.tsv", delim="\t")
terra::writeRaster(esa3rdp_crete, "spatial_data/crete_desertification_risk/esa3rdp_crete.tif",overwrite=TRUE)
esa3rdp_crete_r <- raster(esa3rdp_crete)

check.coords <- points2nearestcell(locations_s, esa3rdp_crete_r)
# merge the points move to closest valued cell
df_corrected_coords <- as.data.frame(check.coords)
sf_desertification <- cbind(st_drop_geometry(metadata_aridity),df_corrected_coords) %>% 
    st_as_sf(coords=c("longitude", "latitude"),
             remove=T,
             crs="WGS84")
# assign the variable to the initial file. The order of the rows is kept the same
metadata_aridity$ESA_id <- raster::extract(esa3rdp_crete_r,sf_desertification, cellnumbers=F)  

metadata_aridity <- metadata_aridity |> left_join(esa3rdp_attr, by=c("ESA_id"="ESA_id"))
metadata_spatial <- metadata_spatial %>% left_join(st_drop_geometry(metadata_aridity))

## world clim data enrichment
raster_map <- function(raster_tmp,world_clim_variable,metadata_spatial, crete_shp){
    raster_pixel <- as(raster_tmp, "SpatialPixelsDataFrame")
    raster_df <- as.data.frame(raster_pixel)

    fill_v <- colnames(raster_df)

    g_base <- ggplot() +
        geom_sf(crete_shp, mapping=aes()) +
        coord_sf(crs="WGS84") +
        theme_bw()
    
    g_dem <- g_base +
        geom_raster(raster_df, mapping=aes(x=x, y=y, fill=.data[[fill_v[1]]]))+
        geom_sf(metadata_spatial, mapping=aes(),color="firebrick", size=0.1, alpha=0.7)
    
    ggsave(paste0("figures/map_",world_clim_variable,"_crete.png",sep=""),
           plot=g_dem, 
           height = 20, 
           width = 30,
           dpi = 300, 
           units="cm",
           device="png")
}

world_clim_dir <- "spatial_data/world_clim_crete/"
world_clim_files <- list.files(world_clim_dir)
# filter the metadata 
metadata_world_clim <- metadata_spatial %>% dplyr::select(ENA_RUN)
for (f in world_clim_files){
    
    # read raster
    raster_path <- paste0(world_clim_dir,f,sep="")
    raster_tmp <- raster(raster_path)
    world_clim_variable <- gsub(".*(bio_\\d{1,2}).*","\\1", f)
    # plot the raster
    raster_map(raster_tmp, world_clim_variable,metadata_spatial, crete_shp)
    print(world_clim_variable)
    # change the location to points outside the raster
    ### 21 points were outside the raster so I use the points2nearestcell 
    ### function of the rSDM package to assign them new location
    locations_s <- locs2sp(do.call(rbind, st_geometry(metadata_world_clim)))
    check.coords <- points2nearestcell(locations_s, raster_tmp)
    # merge the data
    df_corrected_coords <- as.data.frame(check.coords)
    sf_world_clim <- cbind(st_drop_geometry(metadata_world_clim),df_corrected_coords) %>% 
        st_as_sf(coords=c("longitude", "latitude"),
                 remove=T,
                 crs="WGS84")
    # assign the variable to the initial file. The order of the rows is kept the same
    metadata_world_clim[,world_clim_variable] <- raster::extract(raster_tmp, sf_world_clim, cellnumbers=F)  
}

################################### geology ################################
### the original data have EPSF:2001 and encoding = ISO 8859-7. I changed them to 
### EPSG:4326 and encoding UTF-8
crete_geology <- st_read("spatial_data/crete_geology/crete_geology.shp") 
crete_geology <- st_make_valid(crete_geology) # there is an error for an overlapping loop
# 2 points ERR3697698 and ERR3697699 don't overlap with the polygons.
# For polygons there is a function, st_nearest_feature that finds the nearest polygon
# then with cbind join the odjects
nearest_geology = st_nearest_feature(metadata_spatial,crete_geology)
nearest_geology_join = cbind(metadata_spatial, st_drop_geometry(crete_geology)[nearest_geology,]) %>% 
    dplyr::select(ENA_RUN, geology_fo,geology_na,geology_id) %>% st_drop_geometry()

################################## Enrichment ################################
metadata_spatial <- metadata_spatial %>% left_join(nearest_geology_join,join_by(ENA_RUN)) 
metadata_spatial <- metadata_spatial %>% left_join(st_drop_geometry(metadata_world_clim))

metadata_spatial <- st_join(metadata_spatial, clc_crete_shp, left=T) #%>%
metadata_spatial <- st_join(metadata_spatial, natura_crete_land_sci, left=T)
metadata_spatial <- st_join(metadata_spatial, wdpa_crete_wildlife, left=T)
metadata_spatial$dem <- raster::extract(dem_crete, metadata_spatial, cellnumbers=F)

metadata <- metadata_spatial %>% 
    st_drop_geometry() %>% 
    dplyr::select(-c(PER.y, PER.x, MS, INSPIRE_ID, RELEASE_DA,Remark)) %>%
    mutate(elevation_bin=cut(elevation, 
                             breaks=seq.int(from=0, to=2500, by=200),
                             dig.lab = 5 )) %>%
    mutate(protection_status = ifelse(is.na(SITETYPE) & is.na(DESIG_ENG), "none",
                                      ifelse(SITETYPE=="B" & is.na(DESIG_ENG),"Natura2000",
                                             ifelse(is.na(SITETYPE) & DESIG_ENG=="Wildlife Refuge", "Wildlife Refuge", "both"))))

write_delim(metadata,"results/metadata_spatial.tsv", delim="\t")

