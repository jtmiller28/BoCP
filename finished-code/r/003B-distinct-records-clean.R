# 003B-distinct-records-clean
# Author: JT Miller 
# Date: 08-02-2024
# Project: BoCP 

# load libraries 
library(dplyr)
library(data.table)

# set dir 
setwd("/home/millerjared/blue_guralnick/millerjared/BoCP/")
# call edited gatoRs fxns
source("finished-code/r/gatoRs-fxns-edited.R") # edited gatoRs fxns to allow for only pulling idigbio
# Read in the name alignment
name_alignment <- fread("./data/processed/finalized_name_alignment_wcvp.csv")
# Create an accepted name vector & filestyle version.
accepted_name_v <- unique(name_alignment$acceptedNameParent)
# index when necessary
accepted_name_v <- accepted_name_v[19458:32867]
accepted_name_filestyle_v <- gsub(" ", "-", accepted_name_v)

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
    filter(!latitude > 90 & !latitude < -90) %>% 
    filter(!longitude > 180 & !longitude < -180) %>% 
    mutate(roundedLatitude = round(latitude, 2)) %>% 
    mutate(roundedLongitude = round(longitude, 2))
  
  field_cleaned_num_records <- nrow(occur_data)
  
  # Deal with a slight day-month-year issue (redundant, remove)
  # occur_data <- occur_data %>% 
  #   mutate(day = ifelse(day == "", NA, day)) %>% 
  #   mutate(month = ifelse(month == "", NA, month)) %>% 
  #   mutate(year = ifelse(year == "", NA, year))

    
# Identify true duplications of the data (using gatoRs::remove_duplicates with the addition of the collector field (recordedBy) being used. Additionally, Ive added the functionality that prioritizes coordinateUncertainty being lower when choosing from probable duplicates to drop. 
occur_data <- as.data.frame(occur_data) # must only be of df type in order to avoid index issues 
occur_data <- remove_duplicates_mod(occur_data, remove.unparseable =  TRUE)
# Appears that depending on the notation of the data, the collector can vary only slightly (by naming convention, multiple names, etc) so these will be marked as probable duplicates for choice removal. 
occur_data <- occur_data %>% 
  dplyr::group_by(occurrenceID, roundedLatitude, roundedLongitude, scientificName, year, month, day) %>% 
  dplyr::mutate(probableDuplicate = n() > 1) %>% 
  dplyr::ungroup()

  # View the dataframe with the new columns
  occur_data
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