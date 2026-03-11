### Title: World Wide Grids
### Author: JT Miller
### Date: 06/09/2025

### Purpose: Take BoCP's point data, NatureServe's 50x50km P/A grid cell data, and Eastern Canada's 10x10km P/A grid cell data and combine it into a P/A map for the whole world.

# load libraries
library(terra)
library(data.table)
library(dplyr)
library(ggplot2)
library(sf)

## Set up Array tasks for SLURM scheduler
start_num <- as.numeric(Sys.getenv("START_NUM"))
task_id <- as.numeric(start_num)
part <- paste0("part", task_id)

## pull out species for this task 
species_files <- list.files("/blue/guralnick/millerjared/BoCP/data/processed/fully-flagged-data/")
species_file_style <- species_files[task_id]
species_rds_style <- gsub(".csv", ".rds", species_file_style)
species <- gsub("-", " ", species_file_style)
species <- gsub(".csv", "", species)

## Mollweide the earth  & Create flat raster for the NS/CAN data #################################################################################
### Read in the shapefile of the world's botanical regions as a vect (native crs is WGS84)
bot_regions <- terra::vect("./data/raw/level3-wgsrpd/level3.shp")
### project to mollweide 
moll_proj <- "+proj=moll +datum=WGS84 +units=m +no_defs"
bot_regions_moll <- project(bot_regions, moll_proj)
### Extract extent and calc dimensions
moll_ext <- ext(bot_regions_moll)
### Set res to 50,000 m = 50km & 25,000 m = 25km
res_50_km <- 50000 # meters
res_25_km <- 25000 # meters
### Compute the num of rows and cols 
ncol_50 <- ceiling((xmax(moll_ext) - xmin(moll_ext)) / res_50_km)
nrow_50 <- ceiling((ymax(moll_ext) - ymin(moll_ext)) / res_50_km)
ncol_25 <- ceiling((xmax(moll_ext) - xmin(moll_ext)) / res_25_km)
nrow_25 <- ceiling((ymax(moll_ext) - ymin(moll_ext)) / res_25_km)
### Create an empty raster with same extent 
moll_raster_50 <- terra::rast(ncols = ncol_50, nrows = nrow_50, extent = moll_ext, crs = moll_proj)
moll_raster_25 <- terra::rast(ncols = ncol_25, nrows = nrow_25, extent = moll_ext, crs = moll_proj)
### Create a dummy value for visual (not run besides diagnostic check)
# values(moll_raster_50) <- 1
# values(moll_raster_25) <- 1
### plot to check and see how we did (not run besides diagnostic check)
plot(moll_raster_50)
lines(bot_regions_moll, col = "darkred")
# plot(moll_raster_25)
# lines(bot_regions_moll, col = "darkblue")
# res(moll_raster_50)
# res(moll_raster_25)
### do the same methods for our NS/CAN data
## Set up nature serve data
crs_102008 <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 
               +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
# create a 50x50km grid that matches the extent of the NS data, extent was retrieved manually by reviewing the same file mentioned above in QGIS
x_min <- -5105000.0000000000000000
y_min <- -2905000.0000000000000000
x_max<- 3045000.0001220712438226
y_max <- 4645000.0001220703125000

# create an empty raster template for NatureServe Data
ns_raster_template <- terra::rast(
  xmin = x_min, 
  xmax = x_max, 
  ymin = y_min, 
  ymax = y_max, 
  resolution = res_50_km, 
  crs = crs_102008
)

# set up canada data (10x10km)
## load canada shp
canada_grid <- read_sf("./data/raw/grid-data/na-rare-plant-grid-data/NA spatial phylogenetics/NAm_cell_template_10km/NAm_cell_template_10km.shp")
can_raster_template <- terra::rast(ext(canada_grid), 
                                   resolution = 10000, # 10kmx10km
                                   crs = st_crs(canada_grid)$wkt)
##########################################################################################################

### Extract Pt Occurrences for each species ##############################################################
if(any(species_file_style %in% species_files) == TRUE){
  point_occs <- fread(paste0("/blue/guralnick/millerjared/BoCP/data/processed/fully-flagged-data/", species_file_style))
} else {
  point_occs <- NULL # label NULL to avoid subsequent steps
}

## Data cleaning for point occs 
if(!is.null(point_occs)){ # if there are pts, clean
  point_occs <- point_occs %>% 
    filter(taxonomicExactMatch == TRUE) %>% 
    filter(wgs84Datum == TRUE) %>% 
    filter(coordinateIssue == FALSE) %>% 
    filter(trueCoordsWithheld == FALSE) %>% 
    filter(wcvpRangeStatus == "native") %>% # include both datasets native&introduced and native
    filter(validRecord == TRUE) %>% 
    filter(equalLatLon == FALSE) %>% 
    filter(zeroCoords == FALSE) %>% 
    filter(capitalCoord == FALSE) %>% 
    filter(centroidCoord == FALSE) %>% 
    filter(inOceanCoord == FALSE) %>% 
    filter(inGBIFHeadquarters == FALSE) %>% 
    filter(inInstitutionBounds == FALSE)
  if(nrow(point_occs) > 0){
    # Deduplicate data where applicable 
    point_occs <- point_occs %>%
      filter(!amongAggDuplicate | is.na(AggDuplicateGroupID)) %>% # Keep non-duplicates
      bind_rows( # Add in the lowest-ranked duplicates
        point_occs %>%
          filter(amongAggDuplicate) %>%
          group_by(AggDuplicateGroupID) %>%
          filter(AggDuplicateRank == min(AggDuplicateRank, na.rm = TRUE)) %>%
          slice(1) %>%  # if tied for lowest rank, just take the first
          ungroup()
      )
    
    point_occs <- point_occs %>%
      filter(!specimenDuplicate | is.na(specimenDuplicateGroupID)) %>% # Keep non-duplicates
      bind_rows( # Add in the lowest-ranked duplicates
        point_occs %>%
          filter(specimenDuplicate) %>%
          group_by(specimenDuplicateGroupID) %>%
          filter(specimenDuplicateRank == min(specimenDuplicateRank, na.rm = TRUE)) %>%
          slice(1) %>%  # if tied for lowest rank, just take the first
          ungroup()
      )
  }
  # additionally, depending on the resolution we want to filter out coordinate uncertainty that we can live with
  if(nrow(point_occs) > 0){
    point_occs_50 <- point_occs  %>% 
      filter(coordinateUncertaintyInMeters <= 25000 | is.na(coordinateUncertaintyInMeters)) # we're just going with half as a rule
    point_occs_25 <- point_occs %>% 
      filter(coordinateUncertaintyInMeters <= 12500 | is.na(coordinateUncertaintyInMeters))
  } else{
    point_occs_50 <- NULL
    point_occs_25 <- NULL
  }
  # finally assign to sf object 
  # if(!is.null(point_occs_50)){
  #   point_occs_50 <- st_as_sf(point_occs_50, coords = c("roundedLongitude", "roundedLatitude"), crs = 4326)
  # } else {
  #   point_occs_50 <- NULL
  # }
  # if(!is.null(point_occs_25)){
  #   point_occs_25 <- st_as_sf(point_occs_25, coords = c("roundedLongitude", "roundedLatitude"), crs = 4326)
  # } else {
  #   point_occs_25 <- NULL
  # }
}
# Log the number of outliers for further reference
if(!is.null(point_occs)){
  hold_df <- data.frame(
    species = species,
    n_outliers = nrow(filter(point_occs, distOutlier == TRUE))
  )
  fwrite(hold_df, file = "/home/millerjared/blue_guralnick/millerjared/BoCP/outputs/sp-occs-outliers.csv", 
         append = TRUE, col.names = !file.exists("/home/millerjared/blue_guralnick/millerjared/BoCP/outputs/sp-occs-outliers.csv"))
  
}

##########################################################################################################

### Convert the pts data to raster P/A format ############################################################
# assign the spatial data to a spatvector using terra
if(!is.null(point_occs_50)){
pts_wgs_50 <- terra::vect(point_occs_50, geom = c("roundedLongitude", "roundedLatitude"), crs = "EPSG:4326")
} else{
  pts_wgs_50 <- NULL
}
if(!is.null(point_occs_25)){
  pts_wgs_25 <- terra::vect(point_occs_25, geom = c("roundedLongitude", "roundedLatitude"), crs = "EPSG:4326")
} else{
  pts_wgs_25 <- NULL
}
# reproject the pt data to Mollweide 
if(!is.null(pts_wgs_50)){
pts_moll_50 <- project(pts_wgs_50, crs(moll_raster_50))
} else{
  pts_moll_50 <- NULL
}
if(!is.null(pts_wgs_25)){
  pts_moll_25 <- project(pts_wgs_25, crs(moll_raster_25))
} else{
  pts_moll_25 <- NULL
}
# rasterize the pts (P = 1, A = 0)
if(!is.null(pts_moll_50)){
  p_raster_50 <- rasterize(pts_moll_50, moll_raster_50, field = 1, fun = "max", touches = TRUE, background = 0)
} else{
  p_raster_50 <- moll_raster_50 # overlay empty raster
  values(p_raster_50) <- 0 # assign only absences 
}
if(!is.null(pts_moll_25)){
  p_raster_25 <- rasterize(pts_moll_25, moll_raster_25, field = 1, fun = "max", touches = TRUE,  background = 0)
} else{
  p_raster_25 <- moll_raster_25 # overlay empty raster
  values(p_raster_25) <- 0 # assign only absences 
}

### Check plotting, DNE each time, just a diagnostic ###############################################################
r_df <- as.data.frame(p_raster_25, xy = TRUE)
colnames(r_df)[3] <- "presence"
library(sf)
# Convert points to sf for ggplot
if(!is.null(pts_moll_25)){
  pts_sf <- sf::st_as_sf(pts_moll_25)
} else {
  pts_sf <- NULL
}

# Plot with ggplot2
bbox <- st_bbox(pts_sf)
ggplot() +
  geom_raster(data = r_df, aes(x = x, y = y, fill = factor(presence))) +
  scale_fill_manual(values = c("0" = "lightgrey", "1" = "darkgreen"), name = "Presence") +
  geom_sf(data = pts_sf, color = "red", size = 1, alpha = 0.6) +
  coord_sf(xlim = c(bbox["xmin"] - 10000, bbox["xmax"] + 10000), # add slight buffers so its viewable
                         ylim = c(bbox["ymin"] - 10000, bbox["ymax"] + 10000), crs = crs(p_raster_25))
  theme_minimal() +
  ggtitle("Presence/Absence Raster with Occurrence Points")
######################################################################################################################################

#### Convert the Nature Serve 50x50km centroid data to Raster P/A 50x50km and 25x25km format #########################################
# load centroid data from NatureServe for the species at hand 
ns_sp_files <- list.files("/blue/guralnick/millerjared/BoCP/data/raw/grid-data/na-rare-plant-grid-data/NA spatial phylogenetics/ns-grids-species/")
  if(any(species_rds_style %in% ns_sp_files == TRUE)){
    ns_cen_data <- readRDS(paste0("/blue/guralnick/millerjared/BoCP/data/raw/grid-data/na-rare-plant-grid-data/NA spatial phylogenetics/ns-grids-species/", species_rds_style))
    ns_cen_data <- ns_cen_data %>% rename(species = alignedParentName)
  } else{
    ns_cen_data <- NULL
  }
# load the projection string for the NatureServe data
crs_102008 <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs" # extracted from QGIS
# 1. Load and assign CRS to the vector data
if (!is.null(ns_cen_data)) {
  ns_cen_vec <- vect(ns_cen_data)
  crs(ns_cen_vec) <- crs_aea
} else {
  ns_cen_vec <- NULL
}

# 2. Reproject vector to Mollweide **before** rasterizing
if (!is.null(ns_cen_vec)) {
  ns_cen_vec_moll <- project(ns_cen_vec, crs_moll)
} else {
  ns_cen_vec_moll <- NULL
}

# 3. Define a Mollweide raster template with 50x50 km resolution
# Define extent manually or from your full study area
alb_reproj_rast <- project(ns_raster_template, crs_moll)
ext <- ext(moll_raster_50)
res_km <- 50000 # 50 km in meters
moll_template <- rast(ext = ext, resolution = res_km, crs = moll_proj)

# 4. Rasterize in Mollweide space (binary presence-absence)
# Field = 1 for presence, and use "max" to avoid fractional values
if (!is.null(ns_cen_vec_moll)) {
  ns_p_moll_raster_50 <- rasterize(ns_cen_vec_moll, moll_template, field = 1, fun = "max", background = 0)
} else {
  ns_p_moll_raster_50 <- NULL
}
# First, check out the intial rasterization of the centroid data
r_df <- as.data.frame(ns_p_raster_50, xy = TRUE)
colnames(r_df)[3] <- "presence"
# Convert points to sf for ggplot
if(!is.null(ns_cen_vec)){
  pts_sf <- sf::st_as_sf(ns_cen_vec, geom = "geometry", crs = crs_102008)
  pts_sf <- st_transform(pts_sf, crs = st_crs(ns_raster_template))
} else {
  pts_sf <- NULL
}

# Plot with ggplot2
bbox <- st_bbox(pts_sf)
ggplot() +
  geom_raster(data = r_df, aes(x = x, y = y, fill = factor(presence))) +
  scale_fill_manual(values = c("0" = "lightgrey", "1" = "darkgreen"), name = "Presence") +
  geom_sf(data = pts_sf, color = "red", size = 1, alpha = 0.6) +
  coord_sf(xlim = c(bbox["xmin"] - 1000000, bbox["xmax"] + 1000000), # add slight buffers so its viewable
           ylim = c(bbox["ymin"] - 1000000, bbox["ymax"] + 1000000))
theme_minimal() +
  ggtitle("Presence/Absence Raster with Occurrence Points")

## Second plot the reprojection to ensure things are working as intended (Tal on round 0.5, resample)
r_df <- as.data.frame(ns_p_moll_raster_50, xy = TRUE)
colnames(r_df)[3] <- "presence"
#r_df <- r_df %>% mutate(presence = round(presence, 0))
# Convert points to sf for ggplot
if(!is.null(ns_cen_vec)){
  pts_sf <- sf::st_as_sf(ns_cen_vec, geom = "geometry", crs = crs_102008)
  pts_sf <- st_transform(pts_sf, crs = st_crs(moll_raster_50))
} else {
  pts_sf <- NULL
}

# Plot with ggplot2
bbox <- st_bbox(pts_sf)
ggplot() +
  geom_raster(data = r_df, aes(x = x, y = y, fill = factor(presence))) +
  scale_fill_manual(values = c("0" = "lightgrey", "1" = "darkgreen"), name = "Presence") +
  geom_sf(data = pts_sf, color = "red", size = 1, alpha = 0.6) +
  coord_sf(xlim = c(bbox["xmin"] - 1000000, bbox["xmax"] + 1000000), # add slight buffers so its viewable
           ylim = c(bbox["ymin"] - 1000000, bbox["ymax"] + 1000000))
theme_minimal() +
  ggtitle("Presence/Absence Raster with Reproj into Mollweide Occurrence Points")
res(ns_p_moll_raster_50) # 50x50

### Now, we need to resample. 
r_50_to_50 <- resample(ns_p_moll_raster_50, moll_raster_50, method = "max")
r_50_to_25 <- resample(ns_p_moll_raster_50, moll_raster_25, method = "max")

### Further ggplotting
r_df <- as.data.frame(r_50_to_25, xy = TRUE)
colnames(r_df)[3] <- "presence"
# Convert points to sf for ggplot
if(!is.null(ns_cen_vec)){
  pts_sf <- sf::st_as_sf(ns_cen_vec, geom = "geometry", crs = crs_102008)
  pts_sf <- st_transform(pts_sf, crs = st_crs(r_50_to_25))
} else {
  pts_sf <- NULL
}
bbox <- st_bbox(pts_sf)
bot_regions_moll_sf <- st_as_sf(bot_regions_moll)
bot_regions_moll_sf <- st_transform(bot_regions_moll_sf, crs = moll_proj)
ggplot() +
  
  geom_raster(data = r_df, aes(x = x, y = y, fill = factor(presence))) +
  geom_sf(bot_regions_moll_sf, mapping = aes(), fill = NA, color = "black", size = 0.5) +
  scale_fill_manual(values = c("0" = "lightgrey", "1" = "darkgreen"), name = "Presence")  +
  geom_sf(data = pts_sf, color = "red", size = 1, alpha = 0.6) +
  coord_sf(xlim = c(bbox["xmin"] - 1000000, bbox["xmax"] + 1000000), # add slight buffers so its viewable
           ylim = c(bbox["ymin"] - 1000000, bbox["ymax"] + 1000000))
ggsave("/home/millerjared/test2.png", height = 10, width = 10)
plot(moll_raster_50)
lines(bot_regions_moll, col = "darkred") 
polys(test, col = "black")




### other scratch 
r_df <- as.data.frame(test2, xy = TRUE)
colnames(r_df)[3] <- "presence"
# Convert points to sf for ggplot
if(!is.null(ns_cen_vec)){
  pts_sf <- sf::st_as_sf(ns_cen_buff, geom = "geometry", crs = crs_102008)
  pts_sf <- st_transform(pts_sf, crs = st_crs(ns_p_moll_raster_50))
} else {
  pts_sf <- NULL
}
bbox <- st_bbox(pts_sf)
bot_regions_moll_sf <- st_as_sf(bot_regions_moll)
bot_regions_alb_sf <- st_transform(bot_regions_moll_sf, crs = crs_102008)
ggplot() +
  
  geom_raster(data = r_df, aes(x = x, y = y, fill = factor(presence))) +
  geom_sf(bot_regions_moll_sf, mapping = aes(), fill = NA, color = "black", size = 0.5) +
  scale_fill_manual(values = c("0" = "lightgrey", "1" = "darkgreen"), name = "Presence")  +
  geom_sf(data = pts_sf, color = "red", size = 1, alpha = 0.6) +
   coord_sf(xlim = c(bbox["xmin"] - 1000000, bbox["xmax"] + 1000000), # add slight buffers so its viewable
            ylim = c(bbox["ymin"] - 1000000, bbox["ymax"] + 1000000))

###
# Step 1: Get extent of buffered polygons
ext_moll <- ext(moll_raster_50)

# Step 2: Expand extent to nearest 50 km grid
res_m <- 50000

xmin <- floor(ext_moll[1] / res_m) * res_m
xmax <- ceiling(ext_moll[2] / res_m) * res_m
ymin <- floor(ext_moll[3] / res_m) * res_m
ymax <- ceiling(ext_moll[4] / res_m) * res_m

# Step 3: Create a snapped extent
snapped_ext <- ext(xmin, xmax, ymin, ymax)

# Step 4: Generate aligned raster template
moll_raster_template_ns_50 <- rast(ext = snapped_ext, resolution = res_m, crs = moll_proj)

ns_p_moll_raster_50 <- rasterize(ns_cen_buff_moll, moll_raster_template_ns_50, 
                                 field = 1, fun = "max", background = 0)



### Leaflet Check
library(leaflet)

# Reproject original centroids and buffers to WGS84
ns_cen_vec_ll <- project(ns_cen_vec, "EPSG:4326")
ns_cen_buff_ll <- project(ns_cen_buff, "EPSG:4326")
# Convert SpatVectors to sf objects for leaflet
ns_cen_vec_sf <- sf::st_as_sf(ns_cen_vec_ll)
ns_cen_buff_sf <- sf::st_as_sf(ns_cen_buff_ll)


ns_p_moll_raster_50_ll <- project(ns_p_moll_raster_50, "epsg:4326", method = "near")

leaflet() %>%
  addTiles(group = "Base Map") %>%
  addPolygons(data = ns_cen_buff_sf,
              color = "blue", weight = 1, fillOpacity = 0.2, group = "Buffered Cells") %>%
  addCircles(data = ns_cen_vec_sf,
             radius = 100, color = "black", group = "Centroids") %>%
  addRasterImage(ns_p_moll_raster_50_ll, colors = c("transparent", "red"), 
                 opacity = 0.5, project = FALSE, group = "Presence Raster") %>%
  addLayersControl(
    overlayGroups = c("Base Map", "Buffered Cells", "Centroids", "Presence Raster"),
    options = layersControlOptions(collapsed = FALSE)
  )

## 
# 1. Extract original extent
moll_ext <- ext(bot_regions_moll)

# 2. Snap xmin and ymin to nearest 50km grid cell origin
snap_origin <- function(x, resolution) floor(x / resolution) * resolution

xmin_snapped <- snap_origin(xmin(moll_ext), res_50_km)
ymin_snapped <- snap_origin(ymin(moll_ext), res_50_km)
xmax_snapped <- ceiling(xmax(moll_ext) / res_50_km) * res_50_km
ymax_snapped <- ceiling(ymax(moll_ext) / res_50_km) * res_50_km

# 3. Rebuild the extent from snapped values
snapped_ext <- ext(xmin_snapped, xmax_snapped, ymin_snapped, ymax_snapped)

# 4. Build the template raster (now aligned)
moll_raster_50 <- rast(ext = snapped_ext, resolution = res_50_km, crs = moll_proj)

# Convert raster to polygons for visual debugging
grid_polygons <- as.polygons(moll_raster_50)
grid_polygons_sf <- sf::st_as_sf(grid_polygons)
grid_polygons_sf <- st_transform(grid_polygons_sf, crs = 4326)

leaflet() %>%
  addTiles() %>%
  addPolygons(data = grid_polygons_sf, fill = FALSE, weight = 1, color = "gray") %>%
  addPolygons(data = ns_cen_buff_sf, fillOpacity = 0.3, color = "blue")


###
### format as spatvector, assign projection (original natureserve one)
if (!is.null(ns_cen_data)) {
  ns_cen_vec <- vect(ns_cen_data)
  crs(ns_cen_vec) <- crs_102008
} else {
  ns_cen_vec <- NULL
}
### reproj spatvector into Mollweide 
if (!is.null(ns_cen_vec)) {
  ns_cen_vec_moll <- project(ns_cen_vec, moll_proj)
} else {
  ns_cen_vec_moll <- NULL
}

### Buffer each centroid by 25 km Half width (to preserve the original cell level presence)
if(!is.null(ns_cen_vec_moll)){
  ns_cen_vec_moll_buff <- buffer(ns_cen_vec_moll, width = 25000) # 25,000 meters = 25km 
} else{
  ns_cen_vec_moll_buff <- NULL
}

### generate mollweide raster template after reproj
ext_moll <- ext(moll_raster_50)
res50_km <- 50000 # 50km
moll_raster_template_ns_50 <- rast(ext = ext_moll, resolution = res50_km, crs = moll_proj)
### rasterize nature serve data to the 50x50km Mollweide template, note that we use fun = max so that any fractional values are converted to 1 rather than dropped (due to clipping in reproj)
if (!is.null(ns_cen_vec_moll_buff)) {
  ns_p_moll_raster_50 <- rasterize(ns_cen_vec_moll_buff, moll_raster_template_ns_50, field = 1, fun = "max", background = 0)
} else {
  ns_p_moll_raster_50 <- NULL
}

### Resample grid to account for overlap
r_50_to_50 <- resample(ns_p_moll_raster_50, moll_raster_50, method = "max")

### Final Check Plot Diagnostics (DNR)

r_df <- as.data.frame(ns_p_moll_raster_50, xy = TRUE)
colnames(r_df)[3] <- "presence"
# Convert points to sf for ggplot
if(!is.null(ns_cen_vec)){
  pts_sf <- sf::st_as_sf(ns_cen_vec, geom = "geometry", crs = crs_102008)
  pts_sf <- st_transform(pts_sf, crs = st_crs(ns_p_moll_raster_50))
} else {
  pts_sf <- NULL
}
bbox <- st_bbox(pts_sf)
bot_regions_moll_sf <- st_as_sf(bot_regions_moll)
bot_regions_moll_sf <- st_transform(bot_regions_moll_sf, crs = moll_proj)
ggplot() +
  
  geom_raster(data = r_df, aes(x = x, y = y, fill = factor(presence))) +
  geom_sf(bot_regions_moll_sf, mapping = aes(), fill = NA, color = "black", size = 0.5) +
  scale_fill_manual(values = c("0" = "lightgrey", "1" = "darkgreen"), name = "Presence")  +
  geom_sf(data = pts_sf, color = "red", size = 1, alpha = 0.6) +
  coord_sf(xlim = c(bbox["xmin"] - 1000000, bbox["xmax"] + 1000000), # add slight buffers so its viewable
           ylim = c(bbox["ymin"] - 1000000, bbox["ymax"] + 1000000))
ggsave("/home/millerjared/test2.png", height = 10, width = 10)

