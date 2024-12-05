### Flagging Outliers  
# Author: JT Miller 
# Date: 10-31-2024
# Project: BoCP

# Load packages
library(data.table)
library(sf)
library(dplyr)
library(ggplot2)
library(CoordinateCleaner)
library(spThin)
library(gridExtra)

# Set up array logic
start_num <- as.numeric(Sys.getenv("START_NUM"))
task_id <- as.numeric(start_num)
part <- paste0("part", task_id)


# use name alignment to read in data
na_plants_alignment <- fread("./data/processed/finalized_name_alignment_wcvp.csv")

# build a vector of names to index by, then index according to array tasks
accepted_name_v <- unique(na_plants_alignment$acceptedNameParent)
accepted_name <- accepted_name_v[task_id] # index by task ID
accepted_name_filestyle <- gsub(" ", "-", accepted_name)

# Check for input file
file_path <- file.path(paste0("./data/processed/regioned-data/native/", accepted_name_filestyle, ".csv"))
if (!file.exists(file_path)) {
  message("No occurrence data for ", accepted_name)
  return(NULL)
} else{
  print(paste("Found", accepted_name))
  # load the data
  occur_data <- fread(paste0("./data/processed/coord-clean-data/", accepted_name_filestyle, ".csv"))
}

occur_data <- mutate(occur_data, id = 1:n())
occur_data <- mutate(occur_data, parentTaxon = accepted_name)
# Perform outlier detection
id_dist <- CoordinateCleaner::cc_outl(occur_data, lon = "roundedLongitude", lat = "roundedLatitude", species = "parentTaxon",
                                      method = "distance", mltpl = 5, tdi = 250, value = "flagged", sampling_thresh = 0, 
                                      verbose = TRUE, min_occs = 7)

# Invert flag values
id_dist <- !id_dist 

# Append flags as a boolean value
occur_data$distOutlierFlagged <- id_dist 

fwrite(occur_data, paste0("./data/processed/outlier-flagged-distM-data/", accepted_name_filestyle, ".csv"))

