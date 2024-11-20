# Full spatial clean prototype 
# Author: JT Miller 
# Date: 08-03-2024
# Project: BoCP 

## Load libraries 
library(data.table)
library(dplyr)
library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(sp)

## Load in alignment to create a vector for pulling relevant data
# Read in the name alignment
name_alignment <- fread("./data/processed/finalized_name_alignment_wcvp.csv")
# Create an accepted name vector & filestyle version.
accepted_name_v <- unique(name_alignment$acceptedNameParent)
accepted_name_filestyle_v <- gsub(" ", "-", accepted_name_v)


## Compute the world grid, commented out as this should only be ran once to make the rds grid object, else call. 
# world <- ne_countries(scale = "medium", returnclass = "sf")
# world_mollweide <- st_transform(world, crs="+proj=moll +datum=WGS84")
# 
# word_m_grid <- st_make_grid(
#   x = world_mollweide, 
#   cellsize = 5000 # 5km 
# )
# 
# # save the grid as an RDS file 
# saveRDS(word_m_grid, file = "./data/processed/world_m_grid.rds")

world_m_grid <- readRDS("./data/processed/world_m_grid.rds")

world_m_grid_sf <- st_sf(grid_id = 1:length(world_m_grid),
                         geometry = world_m_grid)

for(i in 1:length(accepted_name_v)){
  if(file.exists(paste0("./data/processed/soft-spatial-clean-data/", accepted_name_filestyle_v[i], ".csv"))){
    occur_data <- fread(paste0("./data/processed/soft-spatial-clean-data/", accepted_name_filestyle_v[i], ".csv"))
  } else{
    print(paste("No occur data for", accepted_name_v[[i]]))
  } # end of else no data. 
  init_num_records <- nrow(occur_data)
  # set data as spatial & transform to mollweide projection
  occur_data_sf <- st_as_sf(occur_data, coords = c("longitude", "latitude"), crs = 4326)
  occur_data_sf <- st_transform(occur_data_sf, crs = st_crs(world_m_grid))
  
  # grid the occurrence data
  occur_data_gridded <- st_join(occur_data_sf, world_m_grid_sf, join = st_intersects)
  
  # Only retain one record per grid cell, also prioritize the occurrence with the least amount of associated uncertainty. 
  occur_data_gridded <- occur_data_gridded %>%
    group_by(grid_id) %>% 
    arrange(coordinateUncertaintyInMeters) %>%
    distinct(grid_id, .keep_all = TRUE) %>% 
    ungroup() %>% 
    st_drop_geometry()
  
  # write out report and data
  full_spatial_thin_num_records <- nrow(occur_data_gridded)
  
  full_spatial_clean_info_tb <- data.frame(
    name = accepted_name_v[i],
    priorFullSpatialCleanNumRecords = init_num_records,
    postFullSpatialCleanNumRecords = full_spatial_thin_num_records
  )
  
  fwrite(occur_data_gridded, paste0("./data/processed/full-spatial-clean-data/", accepted_name_filestyle_v[i], ".csv"))
  fwrite(full_spatial_clean_info_tb, file = "./data/processed/cleaning-summaries/full-spatial-clean-info-tb.csv", append = TRUE, col.names = !file.exists("./data/processed/cleaning-summaries/full-spatial-clean-info-tb.csv"))
  print(i) # print where we're at in the loop...
  
  
} # end of for loop