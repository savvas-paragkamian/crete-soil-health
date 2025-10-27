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
library(rSDM) # to find nearest cells of rasters


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

######### crete
crete_shp <- sf::st_read("spatial_data/crete/crete.shp")
crete_peaks <- read_delim("spatial_data/crete_mountain_peaks.csv", delim=";", col_names=T) %>%
    st_as_sf(coords=c("X", "Y"),
             remove=F,
             crs="WGS84")

########### crete bbox

# expand the bbox 
crete_bbox_b <- st_bbox(crete_shp) + c(-1, -1, 1, 1) 
crete_bbox_b_ext <- st_sf(st_as_sfc(crete_bbox_b))

######### other spatial data
clc_crete_shp <- st_read("spatial_data/clc_crete_shp/clc_crete_shp.shp")
natura_crete <- sf::st_read("spatial_data/natura2000/natura2000_crete.shp")
wdpa_crete <- sf::st_read("spatial_data/wdpa_crete/wdpa_crete.shp")
wdpa_crete_wildlife <- wdpa_crete %>% filter(DESIG_ENG=="Wildlife Refugee") %>%
    mutate(DESIG_ENG = gsub("Wildlife Refugee", "Wildlife Refuge", DESIG_ENG))

natura_crete_land <- st_intersection(natura_crete, crete_shp)
natura_crete_land_sci <- natura_crete_land %>% filter(SITETYPE=="B")

##################################################################
# raster DEM hangling
##################################################################

dem_crete <- rast("spatial_data/dem_crete/dem_crete.tif")

### dem
dem <- extract(dem_crete, metadata_spatial, cellnumbers=F)[, -1, drop = FALSE]

metadata_spatial <- cbind(metadata_spatial,dem)

##################################################################
## Global Aridity Index and Potential Evapotranspiration Database
##################################################################

aridity_crete <- rast("spatial_data/crete_aridity_index.tif")
# some points have 0 value
aridity_crete[aridity_crete[] == 0 ] = NA
metadata_aridity <- metadata_spatial %>% dplyr::select(ENA_RUN) 
locations_s <- locs2sf(do.call(rbind, st_geometry(metadata_aridity)))

# merge the points move to closest valued cell
check.coords <- points2nearestcell(locations_s, aridity_crete, layer=1)


sf_aridity <- cbind(st_drop_geometry(metadata_aridity),check.coords) |>
    st_as_sf(
             remove=T,
             crs="WGS84")

# assign the variable to the initial file. The order of the rows is kept the same
aridity_values <- terra::extract(aridity_crete,sf_aridity)[, -1, drop = FALSE]  
colnames(aridity_values) <- "aridity"

metadata_aridity <- cbind(metadata_aridity,aridity_values)

# tranform the values to the original form as instructed by the manual. 
metadata_aridity$aridity <- metadata_aridity$aridity*0.0001

metadata_aridity$aridity_class <- cut(metadata_aridity$aridity,
                                      breaks=c(0,0.03,0.2,0.5, 0.65,0.9),
                                      labels=c("Hyper Arid", "Arid", "Semi-Arid", "Dry sub-humid", "Humid"))

##################################################################
# --------------------Desertification risk --------------------- #
##################################################################

esa3rdp <- rast("spatial_data/ESDAC_CATENA_Desertification2018_RasterFiles/esa3rdp/w001001.adf")
esa3rdp_wgs <- terra::project(esa3rdp, crs(metadata_aridity))

esa3rdp_crete <- crop(esa3rdp_wgs, crete_bbox_b_ext)

esa3rdp_attr <- data.frame(ESA=unique(esa3rdp_crete),
                           ESA_id=seq(1,nrow(unique(esa3rdp_crete)),1))
write_delim(esa3rdp_attr, "spatial_data/crete_desertification_risk/esa3rdp_crete.tsv", delim="\t")
terra::writeRaster(esa3rdp_crete, "spatial_data/crete_desertification_risk/esa3rdp_crete.tif",overwrite=TRUE)

#esa3rdp_crete_r <- rast(esa3rdp_crete)

check.coords <- points2nearestcell(locations_s, esa3rdp_crete, distance=2000)

sf_desertification <- cbind(st_drop_geometry(metadata_aridity),check.coords) %>% 
    st_as_sf(
             remove=T,
             crs="WGS84")

# assign the variable to the initial file. The order of the rows is kept the same
desertification_values <- terra::extract(esa3rdp_crete,sf_desertification)[, -1, drop = FALSE]  
#colnames(desertification_values) <- "ESA_12CL"

# bind with the latest metadata file
metadata_aridity <- cbind(metadata_aridity,desertification_values)

##################################################################
# ---------------Harmonised world soil database v2-------------- #
##################################################################
hwsd2 <- rast("spatial_data/hwsd2_crete/hwsd2_crete.tif")

# hswd metadata
# with trimws the leading spaces are removed for the values.
HWSD2_wrb4 <- read_delim("spatial_data/hwsd2_crete/HWSD2_D_WRB4.tsv", delim="\t") |>
    mutate(VALUE=trimws(VALUE)) |>
    distinct(VALUE, CODE) 

HWSD2_SMU <- read_delim("spatial_data/hwsd2_crete/HWSD2_SMU.tsv", delim="\t") |>
    distinct(HWSD2_SMU_ID, WRB4) |>
    left_join(HWSD2_wrb4, by=c("WRB4"="CODE"))


# all points fall within a raster cell
check.coords <- points2nearestcell(locations_s, hwsd2)

# assign the variable to the initial file. The order of the rows is kept the same
HWSD2_SMU_ID <- extract(hwsd2,metadata_aridity)[, -1,drop=FALSE]

metadata_aridity <- cbind(metadata_aridity,HWSD2_SMU_ID) 

######################## soil erosion G2 #########################

# Seasonal monitoring of soil erosion at regional scale:
# An application of the G2 model in Crete focusing on agricultural land uses
g2_model_all <- terra::rast("spatial_data/1-s2.0-S0303243413001116-mmc1.kmz")

# the values are rgb values
#  has 4 layers is because it’s an RGBA image
# red, green, blue, alpha

r <- g2_model_all

# Drop alpha where transparent
r[[1:3]][r[[4]] == 0] <- NA

# Keep only RGB
r <- r[[1:3]]
# Convert RGB (0–255) to 0–1 scale for ggplot
r <- r / 255

# Convert to data frame for ggplot
r_df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
colnames(r_df) <- c("x", "y", "R", "G", "B")
r_df$col <- rgb(r_df$R, r_df$G, r_df$B)

# make a df with the categories and the colors based on the publication
t_ha_tbl <- data.frame(
                           t_ha_range = c("0 - 0.5", "0.5 - 1", "1 - 2", "2 - 5", "5 - 10", "10 - 20", ">20"),
                           col  = c("#38A800", "#6FC400", "#B0E000", "#FFFF00", "#FFAA00", "#FF5500", "#FF0000"))

t_ha_tbl$code <- row_number(t_ha_tbl)

# final raster df
g2_model <- r_df |>
    left_join(t_ha_tbl)

# make a SpatVector from the points (same CRS as r)
v <- vect(g2_model, geom = c("x","y"), crs = crs(r))

# make a 1-layer template raster on the same grid as r
tmpl <- r[[1]]            # copy extent/resolution/CRS
values(tmpl) <- NA        # clear values

# rasterize the factor field to the grid (creates a categorical raster)
t_ha_r <- rasterize(v, tmpl, field = "code", fun = "first")

# sanity check that the extent and resolution is the same
all.equal(ext(r), ext(t_ha_r))
res(r) == res(t_ha_r)

# transform raster to categories 
cat_t_ha_r <- as.factor(t_ha_r)

levels(cat_t_ha_r) <- data.frame(value=t_ha_tbl$code, label=t_ha_tbl$t_ha_range)
names(cat_t_ha_r) <- "crete_soil_erosion_g2"

# write raster
writeRaster(cat_t_ha_r,"spatial_data/crete_soil_erosion_g2.tif",overwrite=T)

## check points in raster
## there are 64 points. keeping points only below 1km distance from cell
check.coords <- points2nearestcell(locations_s, cat_t_ha_r, distance=2000)

sf_g2 <- cbind(st_drop_geometry(metadata_aridity),check.coords) %>% 
    st_as_sf(
             remove=T,
             crs="WGS84")

# assign the variable to the initial file. The order of the rows is kept the same
g2_values <- terra::extract(cat_t_ha_r,sf_g2)[, -1, drop = FALSE]  

colnames(g2_values) <- "erosion_g2"

# bind with the latest metadata file
metadata_aridity <- cbind(metadata_aridity,g2_values)

########################## world Clim #############################

metadata_spatial <- metadata_spatial %>% left_join(st_drop_geometry(metadata_aridity))

## world clim data enrichment
#raster_map <- function(raster_tmp, world_clim_variable, metadata_spatial, crete_shp){
#    raster_pixel <- as(raster_tmp, "SpatialPixelsDataFrame")
#    raster_df <- as.data.frame(raster_pixel)
#
#    fill_v <- colnames(raster_df)
#
#    g_base <- ggplot() +
#        geom_sf(crete_shp, mapping=aes()) +
#        coord_sf(crs="WGS84") +
#        theme_bw()
#    
#    g_dem <- g_base +
#        geom_raster(raster_df, mapping=aes(x=x, y=y, fill=.data[[fill_v[1]]]))+
#        geom_sf(metadata_spatial, mapping=aes(),color="firebrick", size=0.1, alpha=0.7)
#    
#    ggsave(paste0("figures/map_",world_clim_variable,"_crete.png",sep=""),
#           plot=g_dem, 
#           height = 20, 
#           width = 30,
#           dpi = 300, 
#           units="cm",
#           device="png")
#}

world_clim_dir <- "spatial_data/world_clim_crete/"
world_clim_files <- list.files(world_clim_dir)

# filter the metadata 
metadata_world_clim <- metadata_spatial %>% dplyr::select(ENA_RUN)

for (f in world_clim_files){
    # read raster
    raster_path <- paste0(world_clim_dir,f,sep="")
    raster_tmp <- rast(raster_path)
    world_clim_variable <- gsub(".*(bio_\\d{1,2}).*","\\1", f)
    # plot the raster
    #raster_map(raster_tmp, world_clim_variable,metadata_spatial, crete_shp)
    print(world_clim_variable)
    # change the location to points outside the raster
    ### 21 points were outside the raster so I use the points2nearestcell 
    ### function of the rSDM package to assign them new location
    check.coords <- points2nearestcell(locations_s, raster_tmp)
    # merge the data
    #df_corrected_coords <- as.data.frame(check.coords)
    sf_world_clim <- cbind(st_drop_geometry(metadata_world_clim),check.coords) %>% 
        st_as_sf(
                 remove=T,
                 crs="WGS84")
    # assign the variable to the initial file. The order of the rows is kept the same
    rast_ext <- terra::extract(raster_tmp,sf_world_clim)[, -1, drop = FALSE]  
    colnames(rast_ext) <- world_clim_variable
    metadata_world_clim <- cbind(metadata_world_clim,rast_ext)
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

print("Script finished.")

