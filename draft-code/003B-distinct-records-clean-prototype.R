# 003B-distinct-records-clean-prototype
# Author: JT Miller 
# Date: 08-02-2024
# Project: BoCP 

# load libraries 
library(dplyr)
library(data.table)

# set dir 
setwd("/home/millerjared/blue_guralnick/millerjared/BoCP/")

# Read in the name alignment
name_alignment <- fread("./data/processed/finalized_name_alignment_wcvp.csv")
# Create an accepted name vector & filestyle version.
accepted_name_v <- unique(name_alignment$acceptedNameParent)
accepted_name_filestyle_v <- gsub(" ", "-", accepted_name_v)

# Build a loading bar for simplicity
# Function to print the loading bar
print_loading_bar <- function(current, total) {
  width <- 50
  progress <- current / total
  bar <- paste(rep("=", floor(progress * width)), collapse = "")
  space <- paste(rep(" ", width - floor(progress * width)), collapse = "")
  percentage <- sprintf("%3.0f%%", progress * 100)
  cat(sprintf("\r|%s%s| %s", bar, space, percentage))
}

# init name length
total_names <- length(accepted_name_v)

# load names and take distinct records for relevant fields.
for(i in 1:length(accepted_name_v)){
  if(file.exists(paste0("./data/processed/taxa-cleaned-data/", accepted_name_filestyle_v[i], ".csv"))){
    occur_data <- fread(paste0("./data/processed/taxa-cleaned-data/", accepted_name_filestyle_v[i], ".csv"))
  } else{
    print(paste("No occur data for", accepted_name_v[[i]]))
  } 
  
  init_num_records <- nrow(occur_data) # store the intial num of records
  
  # First clean data for relevant fields: latitude, longitude, & eventDate. 
  occur_data <- occur_data %>% 
    filter(!is.na(longitude)) %>% 
    filter(!is.na(latitude)) %>% 
    filter(!latitude == 0.00) %>% 
    filter(!longitude == 0.00) %>% 
    filter(!is.na(eventDate))
  
  field_cleaned_num_records <- nrow(occur_data)
  
  # take the distinct num of records
  occur_data <- occur_data %>% 
    group_by(latitude, longitude, scientificNameParsed, eventDate) %>% 
    arrange(coordinateUncertaintyInMeters) %>% # just in case there is a duplicate, but uncertainty can be minimized before we take distinct records. 
    distinct(latitude, longitude, scientificNameParsed, eventDate, .keep_all = TRUE)
  
  distinct_num_records <- nrow(occur_data)
  
  distinct_clean_info_tb <- data.frame(
    name = accepted_name_v[i],
    priorDistinctRecords = init_num_records,
    relevantFieldCleanedRecords = field_cleaned_num_records, 
    distinctRecords = distinct_num_records
  )
  fwrite(occur_data, paste0("./data/processed/distinct-clean-data/", accepted_name_filestyle_v[i], ".csv"))
  fwrite(distinct_clean_info_tb, file = "./data/processed/cleaning-summaries/distinct_clean_info_tb.csv", append = TRUE, col.names = !file.exists("./data/processed/cleaning-summaries/distinct_clean_info_tb.csv"))
  print(i) # print where we're at in the loop...
}