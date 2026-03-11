### Title: World Wide Grids
### Author: JT Miller
### Date: 06/09/2025

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
### only use this code to finish stragglers that errored out
# finished <- list.files("/blue/guralnick/millerjared/BoCP/data/processed/grid-project/worldwide/raster-50x50km/")
# finished <- gsub("-50km.tif", "", finished)
# species_files_grab <- gsub(".csv", "", species_files)
# unfinished_species_files <- setdiff(species_files_grab, finished)
# saveRDS(unfinished_species_files, "/blue/guralnick/millerjared/BoCP/data/processed/grid-project/worldwide/native-grids-unfinished-names-tmp.rds")
species_files <- readRDS("/blue/guralnick/millerjared/BoCP/data/processed/grid-project/worldwide/native-grids-unfinished-names-tmp.rds")
species_files <- paste0(species_files, ".csv")
species_file_style <- species_files[task_id]
species_rds_style <- gsub(".csv", ".rds", species_file_style)
species <- gsub("-", " ", species_file_style)
species <- gsub(".csv", "", species)

## Setup Part 1: Mollweide the Earth
### Read in the shapefile of the world's botanical regions as a vect (original crs is WGS84)
bot_regions <- terra::vect("/blue/guralnick/millerjared/BoCP/data/raw/level3-wgsrpd/level3.shp")
### project to mollweide 
moll_proj <- "+proj=moll +datum=WGS84 +units=m +no_defs"
bot_regions_moll <- project(bot_regions, moll_proj)
### Extract extent and calc dimensions
moll_ext <- ext(bot_regions_moll)
### Set res to 50,000 m = 50km & 25,000 m = 25km
res_50_km <- 50000 # meters
res_25_km <- 25000 # meters
res_10_km <- 10000 # meters
### Compute the num of rows and cols 
ncol_50 <- ceiling((xmax(moll_ext) - xmin(moll_ext)) / res_50_km)
nrow_50 <- ceiling((ymax(moll_ext) - ymin(moll_ext)) / res_50_km)
ncol_25 <- ceiling((xmax(moll_ext) - xmin(moll_ext)) / res_25_km)
nrow_25 <- ceiling((ymax(moll_ext) - ymin(moll_ext)) / res_25_km)
ncol_10 <- ceiling((xmax(moll_ext) - xmin(moll_ext)) / res_10_km)
nrow_10 <- ceiling((ymax(moll_ext) - ymin(moll_ext)) / res_10_km)
### Generate the template rasters for 50x50km and 25x25km mollweide rasters of earth
moll_raster_50 <- terra::rast(ncols = ncol_50, nrows = nrow_50, extent = moll_ext, crs = moll_proj)
moll_raster_25 <- terra::rast(ncols = ncol_25, nrows = nrow_25, extent = moll_ext, crs = moll_proj)
moll_raster_10 <- terra::rast(ncols = ncol_10, nrows = nrow_10, extent = moll_ext, crs = moll_proj)
## Setup Part 2: Read in extent and info for Nature Serve and Canada Centroid data 
### Albers Equal Area Original Projection for NatureServe and Can Data
crs_102008 <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
### Extent for NatureServe data
x_min <- -5105000.0000000000000000
y_min <- -2905000.0000000000000000
x_max<- 3045000.0001220712438226
y_max <- 4645000.0001220703125000
### Create an empty raster matching the proj, extent, and res for our NatureServe data
ns_raster_template <- terra::rast(xmin = x_min, xmax = x_max, ymin = y_min, ymax = y_max, resolution = res_50_km, crs = crs_102008
)
### read Can shapefile in (extract extent directly from file)
canada_grid <- read_sf("/blue/guralnick/millerjared/BoCP/data/raw/grid-data/na-rare-plant-grid-data/NA spatial phylogenetics/NAm_cell_template_10km/NAm_cell_template_10km.shp")
can_raster_template <- terra::rast(ext(canada_grid), 
                                   resolution = 10000, # 10kmx10km
                                   crs = st_crs(canada_grid)$wkt)

## Conversion of Data Type 1: Convert Point Occurrence Data into Raster P/A cells on the World Wide Mollweide 50x50 and 25x25 Rasters
### read in the data 
if(any(species_file_style %in% species_files) == TRUE){
  point_occs <- fread(paste0("/blue/guralnick/millerjared/BoCP/data/processed/fully-flagged-data/", species_file_style))
} else {
  point_occs <- NULL # label NULL to avoid subsequent steps
}
### Data cleaning for point occs 
if(!is.null(point_occs)){ # if there are pts, clean
  point_occs <- point_occs %>% 
    filter(taxonomicExactMatch == TRUE) %>% 
    filter(wgs84Datum == TRUE) %>% 
    filter(coordinateIssue == FALSE) %>% 
    filter(trueCoordsWithheld == FALSE) %>% 
    filter(wcvpRangeStatus == "native") %>% 
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
  if(is.null(point_occs_50) || nrow(point_occs_50) == 0){
    point_occs_50 <- NULL
  }
  if(nrow(point_occs_25) == 0 || is.null(point_occs_25)){
    point_occs_25 <- NULL
  }
}
### assign the spatial data to a spatvector using terra
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
### Reproject the pt data to Mollweide 
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
### Rasterize these pts into a P/A format where 1 = Presence, 0 = Absence
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

## Conversion of Data Type 2: Convert Nature Serve 50x50 centroid data to Mollweide Raster 50x50, and 25x25
### Load the NatureServe Centroid data
ns_sp_files <- list.files("/blue/guralnick/millerjared/BoCP/data/raw/grid-data/na-rare-plant-grid-data/NA spatial phylogenetics/ns-grids-species/")
if(any(species_rds_style %in% ns_sp_files == TRUE)){
  ns_cen_data <- readRDS(paste0("/blue/guralnick/millerjared/BoCP/data/raw/grid-data/na-rare-plant-grid-data/NA spatial phylogenetics/ns-grids-species/", species_rds_style))
  ns_cen_data <- ns_cen_data %>% rename(species = alignedParentName) %>%  mutate(presence = 1)

# define source and target projections
source_crs <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
target_crs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

ns_cen_vec <- vect(ns_cen_data)
crs(ns_cen_vec) <- source_crs
ns_cen_vec_moll <- project(ns_cen_vec, target_crs)
ns_cen_buf <- terra::buffer(ns_cen_vec_moll, 25000) # buffer the centroids by 25km so that they would fit the equivalent of a 50km grid cell
r_ns_50 <- terra::rasterize(ns_cen_buf, moll_raster_50, field = "presence", fun = "max", background = 0, touches = TRUE)

# Now resample 
r_50_to_50 <- resample(r_ns_50, moll_raster_50, method = "max")
r_50_to_25 <- resample(r_ns_50, moll_raster_25, method = "max")

} else{ # if there isnt data present, use template to denote its an empty map
  ns_cen_data <- NULL
  r_50_to_50 <- moll_raster_50
  values(r_50_to_50) <- 0 # assign only absences
  r_50_to_25 <- moll_raster_25
  values(r_50_to_25) <- 0 # assign only absences
}
## Conversion of Data Type 3: Convert Canada 10x10 centroid data to Mollweide Raster 50x50, and 25x25
can_sp_files <- list.files("/blue/guralnick/millerjared/BoCP/data/raw/grid-data/na-rare-plant-grid-data/NA spatial phylogenetics/can-grids-species/")
if(any(species_rds_style %in% can_sp_files == TRUE)){
  can_cen_data <- readRDS(paste0("/blue/guralnick/millerjared/BoCP/data/raw/grid-data/na-rare-plant-grid-data/NA spatial phylogenetics/can-grids-species/", species_rds_style))
  can_cen_data <- can_cen_data %>% rename(species = alignedParentName) %>% mutate(presence = 1)

# define source and target projections
source_crs <- crs(canada_grid)
target_crs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
can_cen_vec <- vect(can_cen_data)
crs(can_cen_vec) <- source_crs
can_cen_vec_moll <- project(can_cen_vec, target_crs)
can_cen_buf <- terra::buffer(can_cen_vec_moll, 5000) # buffer the centroids by 25km so that they would fit the equivalent of a 50km grid cell
r_can_10 <- terra::rasterize(can_cen_buf, moll_raster_10, field = "presence", fun = "max", background = 0, touches = TRUE)

# now resample 
r_10_to_50 <- resample(r_can_10, moll_raster_50, method = "max")
r_10_to_25 <- resample(r_can_10, moll_raster_25, method = "max")

} else{
  can_cen_data <- NULL
  r_10_to_50 <- moll_raster_50
  values(r_10_to_50) <- 0 # assign only absences
  r_10_to_25 <- moll_raster_25
  values(r_10_to_25) <- 0 # assign only absences
}
## Combine rasters and finish 
combined_50 <- app(c(r_50_to_50, r_10_to_50, p_raster_50), fun = max, na.rm = TRUE)
names(combined_50) <- "presence"

combined_25 <- app(c(r_50_to_25, r_10_to_25, p_raster_25), fun = max, na.rm = TRUE)
names(combined_25) <- "presence"

## Write out as a raster file the species
writeRaster(combined_50, paste0("/blue/guralnick/millerjared/BoCP/data/processed/grid-project/worldwide/raster-50x50km/", gsub(" ", "-", species), "-50km.tif"), overwrite = TRUE)
writeRaster(combined_25, paste0("/blue/guralnick/millerjared/BoCP/data/processed/grid-project/worldwide/raster-25x25km/", gsub(" ", "-", species), "-25km.tif"), overwrite = TRUE)



# vis (DNR)
# r_df <- as.data.frame(combined_25,  xy = TRUE)
# colnames(r_df)[3] <- "presence"
# # Convert points to sf for ggplot
#   ns_pts_sf <- sf::st_as_sf(ns_cen_buf, geom = "geometry", crs = target_crs)
#   #can_pts_sf <- sf::st_as_sf(can_cen_buf, geom = "geometry", crs = target_crs)
#   pts_sf <- sf::st_as_sf(pts_wgs_50, geom = "geometry", crs = target_crs)
#  # pts_sf <- st_transform(pts_sf, crs = st_crs(moll_raster_50))
# 
# bbox <- st_bbox(ext(ns_pts_sf))
# bot_regions_moll_sf <- st_as_sf(bot_regions_moll)
# bot_regions_moll_sf <- st_transform(bot_regions_moll_sf, crs = moll_proj)
# ggplot() +
#   
#   geom_raster(data = r_df, aes(x = x, y = y, fill = factor(presence))) +
#   geom_sf(bot_regions_moll_sf, mapping = aes(), fill = NA, color = "black", size = 0.5) +
#   scale_fill_manual(values = c("0" = "lightgrey", "1" = "darkgreen"), name = "Presence")  +
#   geom_sf(data = pts_sf, color = "red", alpha = 0.6) +
#   geom_sf(data = ns_pts_sf, color = "steelblue", alpha = 0.6) +
#   #geom_sf(data = can_pts_sf, color = "purple",  alpha = 0.6) +
#   coord_sf(xlim = c(bbox["xmin"] - 1000000, bbox["xmax"] + 1000000), # add slight buffers so its viewable
#            ylim = c(bbox["ymin"] - 1000000, bbox["ymax"] + 1000000))
# ggsave("/home/millerjared/test2.png", height = 10, width = 10)

