# 003D-spThin
# Author: JT Miller 
# Date: 09-18-2024
# Project: BoCP 

library(sf)
library(data.table)
library(spThin)
library(tidyverse)

# Set up array 
task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
part <- paste0("part", task_id)
# set dir 
setwd("/home/millerjared/blue_guralnick/millerjared/BoCP/")
## Load in alignment to create a vector for pulling relevant data
# Read in the name alignment
name_alignment <- fread("./data/processed/finalized_name_alignment_wcvp.csv")
# Create an accepted name vector & filestyle version.
accepted_name_v <- unique(name_alignment$acceptedNameParent)
# Index due to incomplete runs (these should be removed at finished for simplicity)
accepted_name_v <- accepted_name_v[1164:32867]
accepted_name_filestyle_v <- gsub(" ", "-", accepted_name_v)

# Load in the shapefile that delimits the regions of the world. 
level3_regions <- read_sf("./data/raw/level3-wgsrpd/level3.shp")

for(i in 1:length(accepted_name_v)){
  if(file.exists(paste0("./data/processed/coord-clean-data/", accepted_name_filestyle_v[i], ".csv"))){
    occur_data <- fread(paste0("./data/processed/coord-clean-data/", accepted_name_filestyle_v[i], ".csv"))
  } else{
    print(paste("No occur data for", accepted_name_v[[i]]))
    
  } # end of else no data. 
  
  if(nrow(occur_data) < 45000){
  
  # Create a new field for name. 
  occur_data <- occur_data[, name := accepted_name_v[i]]
  num_occur_data_prior_spthin <- nrow(occur_data)
  spThin::thin(loc.data = occur_data, 
               lat.col = "latitude",
               long.col = "longitude",
               spec.col = "name",
               thin.par = 1,
               max.files = 1,
               out.base = paste0(accepted_name_filestyle_v[i]),
               reps = 100, 
               out.dir = "./data/processed/sp-thinned-data/1km/")
  spThin::thin(loc.data = occur_data, 
               lat.col = "latitude",
               long.col = "longitude",
               spec.col = "name",
               thin.par = 5,
               max.files = 1,
               out.base = paste0(accepted_name_filestyle_v[i]),
               reps = 100, 
               out.dir = "./data/processed/sp-thinned-data/5km/")
  
  
  thinned_data_1km <- fread(paste0("./data/processed/sp-thinned-data/1km/", accepted_name_filestyle_v[i], "_thin1.csv"))
  num_occur_data_1km_post_spthin <- nrow(thinned_data_1km)
  thinned_data_5km <- fread(paste0("./data/processed/sp-thinned-data/5km/", accepted_name_filestyle_v[i], "_thin1.csv"))
  num_occur_data_5km_post_spthin <- nrow(thinned_data_5km)
  occur_data_sf <- st_as_sf(occur_data, coords = c("longitude", "latitude"), crs = 4326)
  thinned1km_data_sf <- st_as_sf(thinned_data_1km, coords = c("longitude", "latitude"), crs = 4326)
  thinned5km_data_sf <- st_as_sf(thinned_data_5km, coords = c("longitude", "latitude"), crs = 4326)
  occur_data_hull <- st_convex_hull(occur_data_sf) # Build a convex hull
  occur_data_hull_b <- st_buffer(occur_data_hull, 10)
  hull_bbox <- st_bbox(occur_data_hull_b)
  
  
  
  one_km_plot <- ggplot() + 
    geom_sf(level3_regions, mapping = aes()) +
    geom_sf(occur_data_sf, mapping = aes(color = "blue"),  alpha = 0.5) +
    geom_sf(thinned1km_data_sf,mapping = aes(color = "orange"), alpha = 0.5) + 
    coord_sf(xlim = c(hull_bbox[1], hull_bbox[3]), ylim = c(hull_bbox[2], hull_bbox[4])) +
    theme_bw() +
    ggtitle(paste("1km Thinned Occurrence Data for", accepted_name_v[i])) +
    theme(plot.title=element_text(hjust = 0.5)) +
    xlab("Longitude") +
    ylab("Latitude") + 
    scale_color_manual(name = "Thinning Status", 
                       values = c("blue" = "blue",
                                  "orange" = "orange"),
                       labels = c(paste("Original Data", "n = ", nrow(occur_data_sf)), 
                                  paste("Thinned Data", "n = ", nrow(thinned1km_data_sf))))
  # Save as a jpg
  ggsave(one_km_plot, file = paste0("./outputs/sp-thin-plots/1km/", accepted_name_filestyle_v[i], "thinned.jpg"), height = 10, width = 10)
  
  five_km_plot <- ggplot() + 
    geom_sf(level3_regions, mapping = aes()) +
    geom_sf(occur_data_sf, mapping = aes(color = "blue"),  alpha = 0.5) +
    geom_sf(thinned5km_data_sf,mapping = aes(color = "orange"), alpha = 0.5) + 
    coord_sf(xlim = c(hull_bbox[1], hull_bbox[3]), ylim = c(hull_bbox[2], hull_bbox[4])) +
    theme_bw() +
    ggtitle(paste("5km Thinned Occurrence Data for", accepted_name_v[i])) +
    theme(plot.title=element_text(hjust = 0.5)) +
    xlab("Longitude") + 
    ylab("Latitude") + 
    scale_color_manual(name = "Thinning Status", 
                       values = c("blue" = "blue",
                                  "orange" = "orange"),
                       labels = c(paste("Original Data", "n = ", nrow(occur_data_sf)), 
                                  paste("Thinned Data", "n = ", nrow(thinned5km_data_sf))))
  
  # Save as jpg
  ggsave(five_km_plot, file = paste0("./outputs/sp-thin-plots/5km/", accepted_name_filestyle_v[i], "thinned.jpg"), height = 10, width = 10)
  
  # Create info tb
  # create table 
  sp_thin_info_tb <- data.frame(
    name = accepted_name_v[i],
    priorspThin = num_occur_data_prior_spthin,
    postspThin1km = num_occur_data_1km_post_spthin,
    postspThin5km = num_occur_data_5km_post_spthin)
  
  fwrite(sp_thin_info_tb, file = "./data/processed/cleaning-summaries/sp_thin_info_tb.csv", append = TRUE, col.names = !file.exists("./data/processed/cleaning-summaries/sp_thin_info_tb.csv"))
  } else{
    # Too much data with low mem allocations, report and rerun with larger mem
    sp_thin_too_big_tb <- data.frame(name = accepted_name_v[i], 
                                     tooBig = TRUE, 
                                     numOfRecords = nrow(occur_data))
    fwrite(sp_thin_too_big_tb, file = "./data/processed/cleaning-summaries/spThin_too_big.csv", append = TRUE, col.names = !file.exists("./data/processed/cleaning-summaries/spThin_too_big.csv"))
  }
  print(i)
}