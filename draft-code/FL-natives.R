### A script for determining florida plants from our datasets. 
# Author: JT Miller 
# Date: 09-25-2024
# Project: BoCP 

# Load Libraries
library(data.table)
library(dplyr)

# set dir 
setwd("/home/millerjared/blue_guralnick/millerjared/BoCP/")

# Set up array logic
start_num <- as.numeric(Sys.getenv("START_NUM"))
task_id <- as.numeric(start_num)
part <- paste0("part", task_id)

# Read in all data 
na_plants_alignment <- fread("./data/processed/finalized_name_alignment_wcvp.csv")

# build a vector of names to index by, then index according to array tasks
accepted_name_v <- unique(na_plants_alignment$acceptedNameParent)
accepted_name <- accepted_name_v[task_id] # index by task ID
accepted_name_filestyle <- gsub(" ", "-", accepted_name)

if(file.exists(paste0("./data/processed/regioned-data/all/", accepted_name_filestyle, ".csv"))){
  # load the data
  occur_data <- fread(paste0("./data/processed/regioned-data/all/", accepted_name_filestyle, ".csv"))
  occur_data_fl <- occur_data[area == "Florida"]
  if(nrow(occur_data_fl > 0)){
    fl_region_summary <- data.frame(
      parentTaxon = accepted_name, 
      numOfRecordsInFL = nrow(subset(occur_data_fl, occur_data_fl$area == "Florida")),
      nativeOfFL = unique(occur_data_fl$nativeRangeStatus)
    )
    fwrite(fl_region_summary, "./data/processed/fl-plants-project/name-summary/fl-plant-summary.csv", append = TRUE, col.names = !file.exists("./data/processed/fl-plants-project/name-summary/fl-plant-summary.csv"))
    fwrite(occur_data_fl, paste0("./data/processed/fl-plants-project/occurrence-data/", accepted_name_filestyle, ".csv"))
  }
} else {
  print(paste("no data for", accepted_name))
}