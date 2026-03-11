### Title: Merge Occ and Grid data
### Author: JT Miller
### Date: 04/18/2025

## Load Libraries 
library(tidyverse)
library(data.table)
library(sf)
library(terra)

## Set up Array tasks for SLURM scheduler
start_num <- as.numeric(Sys.getenv("START_NUM"))
task_id <- as.numeric(start_num)
part <- paste0("part", task_id)

## pull out species for this task 
species_files <- list.files("/blue/guralnick/millerjared/BoCP/data/processed/fully-flagged-data/")
### only use this code to finish stragglers that errored out
# finished <- list.files("/blue/guralnick/millerjared/BoCP/data/processed/grid-project/raster-25x25km/")
# finished <- gsub("-25km.tif", "", finished)
# species_files_grab <- gsub(".csv", "", species_files)
# species_files <- setdiff(species_files_grab, finished)
# species_files <- paste0(species_files, ".csv")
############################################################
species_file_style <- species_files[task_id]
species_rds_style <- gsub(".csv", ".rds", species_file_style)
species <- gsub("-", " ", species_file_style)
species <- gsub(".csv", "", species)

## Load in shapefile of world botanical regions to allow for occ data full world extent down the line
bot_regions <- read_sf("./data/raw/level3-wgsrpd/level3.shp")

## Read in relevant data
# point occurrences from BoCP 
if(any(species_file_style %in% species_files) == TRUE){
  point_occs <- fread(paste0("/blue/guralnick/millerjared/BoCP/data/processed/fully-flagged-data/", species_file_style))
} else {
  point_occs <- NULL
}
# centroid data from NatureServe
ns_sp_files <- list.files("/blue/guralnick/millerjared/BoCP/data/raw/grid-data/na-rare-plant-grid-data/NA spatial phylogenetics/ns-grids-species/")
if(any(species_rds_style %in% ns_sp_files == TRUE)){
  ns_cen_data <- readRDS(paste0("/blue/guralnick/millerjared/BoCP/data/raw/grid-data/na-rare-plant-grid-data/NA spatial phylogenetics/ns-grids-species/", species_rds_style))
  ns_cen_data <- ns_cen_data %>% rename(species = alignedParentName)
} else{
  ns_cen_data <- NULL
}
# centroid data from Canada
can_sp_files <- list.files("/blue/guralnick/millerjared/BoCP/data/raw/grid-data/na-rare-plant-grid-data/NA spatial phylogenetics/can-grids-species/")
if(any(species_rds_style %in% can_sp_files == TRUE)){
  can_cen_data <- readRDS(paste0("/blue/guralnick/millerjared/BoCP/data/raw/grid-data/na-rare-plant-grid-data/NA spatial phylogenetics/can-grids-species/", species_rds_style))
  can_cen_data <- can_cen_data %>% rename(species = alignedParentName)
} else{
  can_cen_data <- NULL
}
## Data cleaning for point occs 
if(!is.null(point_occs)){
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
  if(!is.null(point_occs_50)){
    point_occs_50 <- st_as_sf(point_occs_50, coords = c("roundedLongitude", "roundedLatitude"), crs = 4326)
  } else {
    point_occs_50 <- NULL
  }
  if(!is.null(point_occs_25)){
    point_occs_25 <- st_as_sf(point_occs_25, coords = c("roundedLongitude", "roundedLatitude"), crs = 4326)
  } else {
    point_occs_25 <- NULL
  }
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
### Set up Spatial Grids 

## nature serve; assign CRS manually; build 50x50km res to the extent of the original shapefile
crs_102008 <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 
               +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs" # extracted from QGIS
# create a 50x50km grid that matches the extent of the NS data, extent was retrieved manually by reviewing the same file mentioned above in QGIS
x_min <- -5105000.0000000000000000
y_min <- -2905000.0000000000000000
x_max<- 3045000.0001220712438226
y_max <- 4645000.0001220703125000
ns_bbox <- st_bbox(c(xmin = x_min, ymin = y_min, xmax = x_max, ymax = y_max), crs = crs_102008)
ns_grid <- st_make_grid(ns_bbox, cellsize = c(50000, 50000), what = "polygons", square = TRUE)
ns_grid_sf <- st_sf(id = 1:length(ns_grid), geometry = ns_grid)
# reproj the underlying centroid data to the AEA proj that came with the files
if(!is.null(ns_cen_data)){
  ns_cen_data <- st_set_crs(ns_cen_data, crs_102008)
}
## canada; load the 10x10km res grid
canada_grid <- read_sf("/blue/guralnick/millerjared/BoCP/data/raw/grid-data/na-rare-plant-grid-data/NA spatial phylogenetics/NAm_cell_template_10km/NAm_cell_template_10km.shp")
canada_grid <- canada_grid %>% 
  select(-Id, -GRID_ID) %>% 
  mutate(id = 1:n())

## For point data, we need to do some reproj on rasters so it'll come later 

### Write a function for converting whichever layer we're working in into a raster of P/A
sf_to_PA_raster <- function(sf.data, # the SF data being fed
                            res, # the resolution of the data in meters, if point data use the desired resolution instead
                            underlying.grid, # the underlying grid used to build the data (if applicable)
                            outliers = TRUE, # whether we should include outliers
                            map.outputs = FALSE # should an ggplot be generated for the points to presence comparisons?
){
  if(!is.null(sf.data)){
    if(outliers == FALSE){
      sf.data <- sf.data %>% filter(distOutlier == FALSE)
    }
  }
  if(is.null(sf.data)){
    message("sf.data is NULL, creating empty matrix sf object")
    sf.data <- sf::st_sf(
      species = character(0), 
      geometry = sf::st_sfc(crs = sf::st_crs(underlying.grid))
    )
  }
  joined <- st_join(sf.data, underlying.grid, join = st_intersects)
  presence <- unique(joined$id)
  sp_grid <- underlying.grid %>% 
    dplyr::mutate(presence = ifelse(id %in% presence, 1,0)) %>% 
    dplyr::mutate(presence = as.factor(presence))
  r_template <- rast(ext(sp_grid), 
                     resolution = res, 
                     crs = st_crs(sp_grid)$wkt)
  r_presence <- rasterize(vect(sp_grid), 
                          r_template, 
                          field = "presence", 
                          fun = "max")
  
  if(map.outputs == TRUE){
    r_pts <- as.points(r_presence)
    r_pts_sf <- st_as_sf(r_pts)
    bbox <- st_bbox(filter(r_pts_sf, presence == 1))
    sp <- unique(sf.data$species)
    p <- ggplot() +
      geom_sf(data = sf.data, color = "red", size = 1) +
      geom_sf(data = r_pts_sf, aes(color = as.factor(presence)), size = 1, shape = 5) +
      geom_sf(data = underlying.grid, fill = NA, color = "grey") +
      scale_color_manual(values = c("white", "darkgreen")) +
      coord_sf(xlim = c(bbox["xmin"] - 100000, bbox["xmax"] + 100000), # add slight buffers so its viewable
               ylim = c(bbox["ymin"] - 100000, bbox["ymax"] + 100000)) +
      ggtitle(paste0("Presence/Absence raster vs. points for ", sp )) +
      theme_minimal()
    print(p)
  }
  return(r_presence)
  
}

## Use P/A raster function on our three data sets
r_50 <- sf_to_PA_raster(sf.data = ns_cen_data, 
                        res = 50000, 
                        underlying.grid = ns_grid_sf, 
                        outliers = TRUE, # I built this slightly wrong, there are no outliers to filter in the centroided data so act like this is OFF (otherwise it'll cause errors)
                        map.outputs = FALSE)
r_10 <- sf_to_PA_raster(sf.data = can_cen_data,
                        res = 10000, 
                        underlying.grid = canada_grid, 
                        outliers = TRUE, 
                        map.outputs = FALSE)

# Build global P/A matrix for plant occ data originating from North America

## Set up a bounding box extent for typical decimal degrees lat lon of the world (wgs84 datum)
world_ext_latlon <- st_bbox(c(xmin = -180, xmax = 180, ymin = -90, ymax = 90), crs = 4326)

## Load mollweide equal area projection
moll_proj <- "+proj=moll +datum=WGS84 +units=m +no_defs"

## reproj the world latlon bbox 
world_ext_moll <- st_transform(world_ext_latlon, crs = moll_proj)

## build 50x50km grid cells over this extent (the extent of earth)
world_50_grid <- st_make_grid(world_ext_moll, 
                              cellsize = 50000, 
                              square = TRUE)
## build 25x25km grid cells over this extent (the extent of earth)
world_25_grid <- st_make_grid(world_ext_moll, 
                              cellsize = 25000, 
                              square = TRUE)
## convert grids to sf objects, add grid cell ids 
world_50_grid_sf <- st_sf(geometry = world_50_grid) %>% mutate(id = 1:n())
world_25_grid_sf <- st_sf(geometry = world_50_grid) %>% mutate(id = 1:n())

# Use the following code to check if area per cell was correctly made...ignore for running as array
# grid_test <- world_50_grid_sf %>% 
#     mutate( area_m2 = st_area(geometry),
#             area_km2 = as.numeric(area_m2) / 1e6,
#             # Compute xdist (width) and ydist (height)
#             bbox = map(geometry, st_bbox),
#             xdist = map_dbl(bbox, ~ .x["xmax"] - .x["xmin"]),
#             ydist = map_dbl(bbox, ~ .x["ymax"] - .x["ymin"])
# 
#     )

# Transform Occ Pt data to the moll proj that our world grid is made as 
if(!is.null(point_occs_50)){
  point_occs_50 <- st_transform(point_occs_50, crs = st_crs(world_50_grid_sf))
} else{
  point_occs_50 <- NULL
}

if(!is.null(point_occs_25)){
  point_occs_25 <- st_transform(point_occs_25, crs = st_crs(world_25_grid_sf))
} else{
  point_occs_25 <- NULL
}

# Run P/A raster conversion on the pt data
r_point_50 <- sf_to_PA_raster(sf.data = point_occs_50,
                              res = 50000, 
                              underlying.grid = world_50_grid_sf, 
                              outliers = FALSE, 
                              map.outputs = FALSE)

r_point_25 <- sf_to_PA_raster(sf.data = point_occs_25, 
                              res = 25000, 
                              underlying.grid = world_25_grid_sf, 
                              outliers = FALSE,
                              map.outputs = FALSE)

# Make templates of rasters at 50x50 and 25x25 km resolutions to resample to
pt_50_ext <- ext(r_point_50) # extract extent for creation of 50x50km template (will be the same as below)
pt_25_ext <- ext(r_point_25) # extract extent for creation of 25x25km template
template_50 <- rast(ext = pt_50_ext, resolution = 50000, crs = crs(world_50_grid_sf))
template_25 <- rast(ext = pt_25_ext, resolution = 25000, crs = crs(world_25_grid_sf))

# Adjust projections to match that of the mollwiede equal area cells
r_10_moll <- project(r_10, r_point_50)
r_50_moll <- project(r_50, r_point_50)

## Resample
r_50_to_50 <- resample(r_50_moll, template_50, method = "max")
r_50_to_25 <- resample(r_50_moll, template_25, method = "max")
r_10_to_50 <- resample(r_10_moll, template_50, method = "max")
r_10_to_25 <- resample(r_10_moll, template_25, method = "max")
r_point_50_to_50 <- resample(r_point_50, template_50, method = "max")
r_point_25_to_25 <- resample(r_point_25, template_25, method = "max")

## Combine these data 
combined_50 <- app(c(r_50_to_50, r_10_to_50, r_point_50_to_50), fun = max, na.rm = TRUE)
names(combined_50) <- "presence"
combined_25 <- app(c(r_50_to_25, r_10_to_25, r_point_25_to_25), fun = max, na.rm = TRUE)
names(combined_25) <- "presence"

## Write out as a raster file the species
writeRaster(combined_50, paste0("/blue/guralnick/millerjared/BoCP/data/processed/grid-project/raster-50x50km/", gsub(" ", "-", species), "-50km.tif"), overwrite = TRUE)
writeRaster(combined_25, paste0("/blue/guralnick/millerjared/BoCP/data/processed/grid-project/raster-25x25km/", gsub(" ", "-", species), "-25km.tif"), overwrite = TRUE)

## Create an overall ggplot for these data (diagnostic DNR)
r_points_c <- as.points(combined_50)
r_points_sf_c <- st_as_sf(r_points_c)
#bbox_c <- st_bbox(filter(r_points_sf_c, presence == 1))
#contributing_occ_pts <- st_transform(point_occs_50, crs = crs(world_50_grid_sf))
contributing_can_centroids <- st_transform(can_cen_data, crs = crs(world_50_grid_sf))
bbox_c <- st_bbox(contributing_can_centroids)
#contributing_ns_centroids <- st_transform(ns_cen_data, crs = crs(world_50_grid_sf))
ggplot() +
  geom_sf(data = r_points_sf_c, aes(color = as.factor(presence)), size = 1) +
  geom_sf(data = world_50_grid_sf, fill = NA, color = "grey") +
  #geom_sf(data = contributing_occ_pts, color = "red", size = 1) +
  geom_sf(data = contributing_can_centroids, color = "steelblue", size = 1.2) +
  #geom_sf(data = contributing_ns_centroids, color = "goldenrod", size = 1.4) +
  scale_color_manual(values = c("white", "darkgreen")) +
  coord_sf(xlim = c(bbox_c["xmin"] - 100000, bbox_c["xmax"] + 100000), # add slight buffers so its viewable
           ylim = c(bbox_c["ymin"] - 100000, bbox_c["ymax"] + 100000)) +
  theme_minimal()

