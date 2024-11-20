# 003C-soft-spatial-clean-prototype
# Author: JT Miller 
# Date: 08-03-2024
# Project: BoCP 

# Load Libraries 
library(dplyr)
library(data.table)
library(geosphere)

# set dir 
setwd("/home/millerjared/blue_guralnick/millerjared/BoCP/")

# Read in the name alignment
name_alignment <- fread("./data/processed/finalized_name_alignment_wcvp.csv")
# Create an accepted name vector & filestyle version.
accepted_name_v <- unique(name_alignment$acceptedNameParent)
accepted_name_filestyle_v <- gsub(" ", "-", accepted_name_v)

# Create a function that filters out data that is 1km distance to eachother by decimal degrees when longitude is held constant, prioritize records that have the least amount of coordinateUncertainty. 
filter_latitude <- function(group){
  group <- group %>% arrange(coordinateUncertaintyInMeters) # arrange acsending
  lat_diff <- abs(diff(group$latitude)) # compute latitudnal diff 
  keep_indices <- c(TRUE, lat_diff > 0.01) # keep records where lat diff is greater to 0.01 ~ 1km 
  filtered_group <- group[keep_indices, ] # subset by our keep indices 
  return(filtered_group)
}

# Create a function that filters out data that is 1km distance to eachother by decimal degrees when latitude is held constant, prioritize records that have the least amount of coordinateUncertainty. 
filter_longitude <- function(group){
  group <- group %>% arrange(coordinateUncertaintyInMeters) # arrange acsending
  lon_diff <- abs(diff(group$longitude)) # compute latitudnal diff 
  keep_indices <- c(TRUE, lon_diff > 0.01) # keep records where lat diff is greater to 0.01 ~ 1km 
  filtered_group <- group[keep_indices, ] # subset by our keep indices 
  return(filtered_group)
}

# Load in occur data and soft clean for records that are too close to one another to be useful. 
for(i in 1:length(accepted_name_v)){
  if(file.exists(paste0("./data/processed/taxa-cleaned-data/", accepted_name_filestyle_v[i], ".csv"))){
    occur_data <- fread(paste0("./data/processed/taxa-cleaned-data/", accepted_name_filestyle_v[i], ".csv"))
  } else{
    print(paste("No occur data for", accepted_name_v[[i]]))
  }   
  init_num_records <- nrow(occur_data)
  
  occur_data <- occur_data %>% 
    group_by(longitude) %>% 
    group_modify(~ filter_latitude(.x)) %>% 
    ungroup()
  
  occur_data <- occur_data %>% 
    group_by(latitude) %>% 
    group_modify(~ filter_longitude(.x)) %>% 
    ungroup()
  
  soft_spatial_clean_num_records <- nrow(occur_data)
  
  soft_spatial_clean_info_tb <- data.frame(
    name = accepted_name_v[i],
    priorSoftSpatialCleanNumRecords = init_num_records,
    postSoftSpatialCleanNumRecords = soft_spatial_clean_num_records
  )
  
  fwrite(occur_data, paste0("./data/processed/soft-spatial-clean-data/", accepted_name_filestyle_v[i], ".csv"))
  fwrite(soft_spatial_clean_info_tb, file = "./data/processed/cleaning-summaries/soft_spatial_clean_info_tb.csv", append = TRUE, col.names = !file.exists("./data/processed/cleaning-summaries/soft_spatial_clean_info_tb.csv"))
  print(i) # print where we're at in the loop...
  
  
}

