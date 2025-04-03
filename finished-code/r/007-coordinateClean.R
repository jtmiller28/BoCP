# 007 Coordinate Clean
# Author: JT Miller 
# Date: 03-12-2024
# Project: BoCP 

# An additional coordinate clean step using the package Coordinate Cleaner
## Addressing:
### capitals: tests a radius around adm-9 capitals
### centroids: tests a radius around country centroids
### equal: tests for equal absolute lat and lon
### institutions: tests a radius around known biodiversity institutions from institutions
### seas: tests if coordinates fall into the ocean
### zeros: tests for plain zeros, equal latitude and longitude, and a radius around the point 0/0
### Outliers

# Load Libraries 
library(data.table)
library(tidyverse)
library(CoordinateCleaner)

# Set up array logic
start_num <- as.numeric(Sys.getenv("START_NUM"))
task_id <- as.numeric(start_num)
part <- paste0("part", task_id)

# read in temp file containing names that are missing (see name-array-setup.R)
unfinished_names <- readRDS("/blue/guralnick/millerjared/BoCP/data/processed/unfinished_flag_names_temp.rds")
# extract data
names_to_retrieve <- unfinished_names[[task_id]]
accepted_name <- names_to_retrieve[1]
accepted_name_file_style <- gsub(" ", "-", accepted_name)

# read in the species table 
print(paste("retrieving data for", accepted_name))
occ_data <- fread(paste0("/blue/guralnick/millerjared/BoCP/data/processed/dupe-flagged/", accepted_name_file_style, ".csv"))

#  remove data without coords
occ_data_w_coords <- occ_data %>% 
  filter(!is.na(decimalLatitude)) %>% 
  filter(!is.na(decimalLongitude)) %>%  # both must be clean
  filter(coordinateIssue == FALSE) %>% 
  filter(wgs84Datum == TRUE) # coordinate cleaner requires WGS84 records for its checks.

occ_data_wout_coords <- occ_data %>% 
  filter(is.na(decimalLatitude) | is.na(decimalLongitude) | wgs84Datum == FALSE | coordinateIssue == TRUE) # either or

# use coord cleaner

if(nrow(occ_data_w_coords) > 0){
  print(paste("spatially valid data surpassed 0 threshold, proceed with coord cleaning for", accepted_name))
occ_data_w_coords <- CoordinateCleaner::clean_coordinates(occ_data_w_coords, 
                                                   lon = "decimalLongitude", 
                                                   lat = "decimalLatitude", 
                                                   species = "parsedName", 
                                                   tests = c("capitals", "centroids", "equal", "gbif", "institutions", "seas",
                                                             "zeros"))
# rename for my sanity
occ_data_w_coords <- occ_data_w_coords %>% 
  rename(validRecord = `.val`,
         equalLatLon = `.equ`,
         zeroCoords = `.zer`,
         capitalCoord = `.cap`,
         centroidCoord = `.cen`,
         inOceanCoord = `.sea`,
         inGBIFHeadquarters = `.gbf`, 
         inInstitutionBounds = `.inst`) %>% 
  select(-`.summary`) %>% 
  # swap way the flags are presented for simplicity 
  mutate(equalLatLon = !equalLatLon, 
         zeroCoords = !zeroCoords, 
         capitalCoord = !capitalCoord,
         centroidCoord = !centroidCoord,
         inOceanCoord = !inOceanCoord, 
         inGBIFHeadquarters = !inGBIFHeadquarters, 
         inInstitutionBounds = !inInstitutionBounds)
  
# Perform outlier detection 
# convert to data.frame to avoid unknown issues
occ_data_w_coords <- as.data.frame(occ_data_w_coords)
dist_flag <- CoordinateCleaner::cc_outl(occ_data_w_coords, lon = "roundedLongitude", lat = "roundedLatitude", species = "parsedName",
                                      method = "distance", mltpl = 5, tdi = 250, value = "flagged", sampling_thresh = 0, 
                                      verbose = TRUE, min_occs = 7)
# Invert flag values
dist_flag <- !dist_flag

# append to data 
occ_data_w_coords$distOutlier = dist_flag
} else {
  print(paste("spatially valid data is equal to zero, skip coordinate cleaning", accepted_name))
  occ_data_w_coords <- occ_data_w_coords %>% 
    mutate(validRecord = NA,
           equalLatLon = NA, 
           zeroCoords = NA, 
           capitalCoord = NA,
           centroidCoord = NA,
           inOceanCoord = NA, 
           inGBIFHeadquarters = NA, 
           inInstitutionBounds = NA,
           distOutlier = NA)
}
# prep binding the dataset back together via NA values
occ_data_wout_coords <- occ_data_wout_coords %>% 
  mutate(validRecord = NA,
         equalLatLon = NA, 
         zeroCoords = NA, 
         capitalCoord = NA,
         centroidCoord = NA,
         inOceanCoord = NA, 
         inGBIFHeadquarters = NA, 
         inInstitutionBounds = NA,
         distOutlier = NA)

print(paste("bind dataset for", accepted_name))
occ_data <- rbind(occ_data_w_coords, occ_data_wout_coords)
print(paste("write coordinate flagged dataset for", accepted_name))
fwrite(occ_data, paste0("/blue/guralnick/millerjared/BoCP/data/processed/fully-flagged-data/", accepted_name_file_style, ".csv"))

