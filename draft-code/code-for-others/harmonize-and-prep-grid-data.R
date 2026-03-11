### Title: Harmonize & Prep Grid data
### Author: JT Miller
### Date: 04/30/2025

## Load Packages
library(sf)
library(data.table)
library(tidyverse)

## Take Nature Serve centroid data, Harmonize taxonomy with BoCP's, Create species files to be read in piece by piece for optimization
# load ns and harmonized name table
ns_cen_data <- read_sf("/blue/guralnick/millerjared/BoCP/data/raw/grid-data/na-rare-plant-grid-data/NA spatial phylogenetics/tns_points.shp")
ns_names_aligned <- fread("/blue/guralnick/millerjared/BoCP/data/processed/nature-serve-names-aligned.csv")
# harmonize taxonomy via merge, rename original name field
ns_cen_harmonized <- merge(ns_cen_data, ns_names_aligned, by.x = "SNAME", by.y = "user_supplied_name")
ns_cen_harmonized <- ns_cen_harmonized %>% 
  rename(original_ns_name = SNAME)
ns_sp_vector <- unique(ns_cen_harmonized$alignedParentName)
ns_sp_filestyle <- gsub(" ", "-", ns_sp_vector)
# create species files
for(i in 1:length(ns_sp_vector)){
  ns_sp <- ns_cen_harmonized %>% 
    filter(alignedParentName == ns_sp_vector[i])
  write_rds(ns_sp, paste0("/blue/guralnick/millerjared/BoCP/data/raw/grid-data/na-rare-plant-grid-data/NA spatial phylogenetics/ns-grids-species/", ns_sp_filestyle[i], ".rds" ))
}

## Take Canada Centroid data, harmonize taxonomy with BoCP's, Create species files to read in piece by piece for optimization
# load can centroid data and harmonized taxonomy 
can1_cen_data <- read_sf("/blue/guralnick/millerjared/BoCP/data/raw/grid-data/na-rare-plant-grid-data/NA spatial phylogenetics/canada1.shp")
can2_cen_data <- read_sf("/blue/guralnick/millerjared/BoCP/data/raw/grid-data/na-rare-plant-grid-data/NA spatial phylogenetics/canada2.shp")
can_names_aligned <- fread("/blue/guralnick/millerjared/BoCP/data/processed/canada-names-aligned.csv")
# standardize fields for both datasets
can1_cen_data <- can1_cen_data %>% 
  select(NNAME) %>% 
  rename(SNAME = NNAME)
can2_cen_data <- can2_cen_data %>% 
  select(GNAME) %>% 
  rename(SNAME = GNAME)
# combine datasets
canada_centroids <- rbind(can1_cen_data, can2_cen_data)
# harmonize taxonomy 
can_cen_harmonized <- merge(canada_centroids, can_names_aligned, by.x = "SNAME", by.y = "user_supplied_name")
can_cen_harmonized <- can_cen_harmonized %>% 
  rename(original_ns_name = SNAME)
can_sp_vector <- unique(can_cen_harmonized$alignedParentName)
can_sp_filestyle <- gsub(" ", "-", can_sp_vector)
# use loop to write out files 
for(i in 1:length(can_sp_vector)){
  can_sp <- can_cen_harmonized %>% 
    filter(alignedParentName == can_sp_vector[i])
  write_rds(can_sp, paste0("/blue/guralnick/millerjared/BoCP/data/raw/grid-data/na-rare-plant-grid-data/NA spatial phylogenetics/can-grids-species/", can_sp_filestyle[i], ".rds" ))
}
