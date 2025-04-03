### Title: Spatial Indexing 
### Project: BoCP
### Author: JT Miller
### Date: 04/01/2025

# Set up array logic
start_num <- as.numeric(Sys.getenv("START_NUM"))
task_id <- as.numeric(start_num)
part <- paste0("part", task_id)

# some basic code to find out the total number of raster (5km^2) cells each species will occupy
library(data.table)
library(terra)
library(tidyverse)

# grab names 
list_found_names <- list.files("/blue/guralnick/millerjared/BoCP/data/processed/fully-flagged-data/")
list_found_names <- gsub("-", " ", list_found_names)
list_found_names <- gsub(".csv", "", list_found_names)

# extract data
names_to_retrieve <- list_found_names[task_id]
accepted_name <- names_to_retrieve[1]
accepted_name_file_style <- gsub(" ", "-", accepted_name)

# read in the species table 
print(paste("retrieving data for", accepted_name))

test_df <- fread(paste0("/blue/guralnick/millerjared/BoCP/data/processed/fully-flagged-data/", accepted_name_file_style, ".csv"))

test_df_filtered <- test_df %>% 
  filter(taxonomicExactMatch == TRUE) %>% 
  filter(wgs84Datum == TRUE) %>% 
  filter(coordinateIssue == FALSE) %>% 
  filter(trueCoordsWithheld == FALSE) %>% 
  filter(wcvpRangeStatus == "native" | wcvpRangeStatus == "introduced") %>% 
  filter(validRecord == TRUE) %>% 
  filter(equalLatLon == FALSE) %>% 
  filter(zeroCoords == FALSE) %>% 
  filter(capitalCoord == FALSE) %>% 
  filter(centroidCoord == FALSE) %>% 
  filter(inOceanCoord == FALSE) %>% 
  filter(inGBIFHeadquarters == FALSE) %>% 
  filter(inInstitutionBounds == FALSE)


test_df_filtered <- test_df_filtered %>%
  filter(!amongAggDuplicate | is.na(AggDuplicateGroupID)) %>% # Keep non-duplicates
  bind_rows( # Add in the lowest-ranked duplicates
    test_df_filtered %>%
      filter(amongAggDuplicate) %>%
      group_by(AggDuplicateGroupID) %>%
      filter(AggDuplicateRank == min(AggDuplicateRank, na.rm = TRUE)) %>%
      slice(1) %>%  # if tied for lowest rank, just take the first
      ungroup()
  )

test_df_filtered <- test_df_filtered %>%
  filter(!specimenDuplicate | is.na(specimenDuplicateGroupID)) %>% # Keep non-duplicates
  bind_rows( # Add in the lowest-ranked duplicates
    test_df_filtered %>%
      filter(specimenDuplicate) %>%
      group_by(specimenDuplicateGroupID) %>%
      filter(specimenDuplicateRank == min(specimenDuplicateRank, na.rm = TRUE)) %>%
      slice(1) %>%  # if tied for lowest rank, just take the first
      ungroup()
  )

if(nrow(test_df_filtered) > 0){

# load input raster 
# bio_clim_raster <- terra::rast("/blue/soltis/share/FL_Scrub/share/06_rasters/BioClim/wc2.1_30s_bio_1.tif")
# bio_clim_raster_coarse <- terra::aggregate(bio_clim_raster, fact = 3, fun = max)
# writeRaster(bio_clim_raster_coarse, "/blue/guralnick/millerjared/BoCP/data/processed/wc2.1_5km_bio_1.tiff")
bio_clim_raster_coarse <- terra::rast("/blue/guralnick/millerjared/BoCP/data/processed/wc2.1_5km_bio_1.tiff")
# convert spatial points to a spatVector 
pts <- terra::vect(as.data.frame(test_df_filtered), geom = c("roundedLongitude", "roundedLatitude"), crs = "EPSG:4326")
pts <- terra::project(pts, crs(bio_clim_raster_coarse))
# rasterize using unique index per raster cell
test_df_filtered$cell_id <- terra::cellFromXY(bio_clim_raster_coarse, terra::geom(pts)[, c("x", "y")])
length(unique(test_df_filtered$cell_id))
temp_df <- data.frame(
  species = unique(test_df_filtered$species), 
  num_5km_cells = length(unique(test_df_filtered$cell_id))
) 
fwrite(temp_df, file = "/blue/guralnick/millerjared/BoCP/outputs/5km_cells_occupied.csv", append = TRUE, col.names = !file.exists("/blue/guralnick/millerjared/BoCP/outputs/5km_cells_occupied.csv"))
} else {
  temp_df <- data.frame(
    species = accepted_name, 
    num_5km_cells = 0 # nothing after data filters...
  )
  fwrite(temp_df, file = "/blue/guralnick/millerjared/BoCP/outputs/5km_cells_occupied.csv", append = TRUE, col.names = !file.exists("/blue/guralnick/millerjared/BoCP/outputs/5km_cells_occupied.csv"))
  
}