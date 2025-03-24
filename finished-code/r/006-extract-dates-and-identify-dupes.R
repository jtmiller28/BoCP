# 005 Identify Duplicates
# Author: JT Miller 
# Date: 03-08-2024
# Project: BoCP 

### Purpose: extract dates where necessary, then create a flagged field that denotes whether a record is likely a duplicate
#### Duplicates can either be a specimen duplicate or aggregate duplicate
### Specimen duplicates can represent multiple specimens made from one organism 
### Aggregator duplicates can either be duplicates made during the data storing process or multiple aggregators holding the same specimen record

# Load Libraries 
library(data.table)
library(tidyverse)
library(lubridate)

# Set up array logic
start_num <- as.numeric(Sys.getenv("START_NUM"))
task_id <- as.numeric(start_num)
part <- paste0("part", task_id)

# read the taxonomic harmonized list of names + synonyms
name_list <- readRDS("/blue/guralnick/millerjared/BoCP/data/processed/name_list.rds")

# check for names that have occurrence data
list_found_names <- list.files("/blue/guralnick/millerjared/BoCP/data/processed/flagged-species-occs/")
list_found_names <- gsub("-", " ", list_found_names)
list_found_names <- gsub(".csv", "", list_found_names)
# extract data
names_to_retrieve <- list_found_names[[task_id]]
accepted_name <- names_to_retrieve[1]
accepted_name_file_style <- gsub(" ", "-", accepted_name)
print(paste("reading in data for", accepted_name))
# read in the species occurrence table 
occ_data <- fread(paste0("/blue/guralnick/millerjared/BoCP/data/processed/flagged-species-occs/", accepted_name_file_style, ".csv"))
print(paste("Parsing date if possible", accepted_name))
if(class(occ_data$eventDate) == "character"){ # note that this will remove nuanced records. 
  occ_data$eventDate <- as.Date(occ_data$eventDate, format = "%Y-%m-%d")
} else {
  
}
occ_data <- occ_data %>% # note this doesnt solve partials, but too complex for our needs
  mutate( # if there is data present in the day,month,year fields keep it, if not try the event date
    year = ifelse(is.na(year) & !is.na(eventDate), as.numeric(format(as.Date(eventDate), "%Y")), year),
    month = ifelse(is.na(month) & !is.na(eventDate), as.numeric(format(as.Date(eventDate), "%m")), month),
    day = ifelse(is.na(day) & !is.na(eventDate), as.numeric(format(as.Date(eventDate), "%d")), day)
  ) 
# identify duplicates, these can be within aggregate, among aggregates or specimen 
occ_data <- occ_data %>% 
  # group by the unique identifier used within each aggregate, check if any are ever more than one. 
  group_by(uuid, aggregator) %>% 
  mutate(withInAggDuplicate = ifelse(!is.na(uuid) & uuid != "" & n() > 1, TRUE, FALSE)) %>% 
  ungroup() %>% # ungroup for next dupe building
  mutate(occID_lower = tolower(occurrenceID)) %>% # lower occurrenceID just to be sure
  group_by(occID_lower) %>% # group by lowered occurrenceID
  # if there are ever more than one aggregator associated with the same occurrenceID, we're going to assume these are duplicates across aggregates
  # assign a uniq group id to such records 
  mutate(AggDuplicateGroupID = ifelse(n_distinct(aggregator) > 1 & !is.na(occurrenceID), cur_group_id(), NA)) %>% 
  # if a uniq group id exists for this field, its a duplicate so return TRUE
  mutate(amongAggDuplicate = !is.na(AggDuplicateGroupID)) %>% 
  # identify specimenDuplicates via looking whether there are specimens that share occurrenceIDs, all temporal info, and are not aggregate duplicates
  mutate(specimenDuplicateGroupID = ifelse(n_distinct(year) == 1 &
                                             n_distinct(month) == 1 &
                                             n_distinct(day) == 1 &
                                             n() > 1 &
                                             !is.na(occurrenceID) &
                                             amongAggDuplicate == FALSE, 
                                           cur_group_id(), NA)) %>% 
  ungroup() %>% 
  # identify specimen duplicates via looking whether there are specimens without occurrenceIDs that contain the samespatial coords or date collection
  group_by(decimalLatitude, decimalLongitude, day, month, year) %>% 
  mutate(specimenDuplicateGroupID = ifelse(is.na(occurrenceID) &
                                             (sum(!is.na(decimalLatitude)) > 0 & n_distinct(decimalLatitude) == 1) & 
                                             (sum(!is.na(decimalLongitude)) > 0 & n_distinct(decimalLongitude) == 1) & 
                                             (sum(!is.na(day)) > 0 & n_distinct(day) == 1) & 
                                             (sum(!is.na(month)) > 0 & n_distinct(month) == 1) & 
                                             (sum(!is.na(year)) > 0 & n_distinct(year) == 1) & 
                                             n() > 1, 
                                           cur_group_id(), specimenDuplicateGroupID)) %>% 
  ungroup() %>% 
  mutate(specimenDuplicate = ifelse(!is.na(specimenDuplicateGroupID), TRUE, FALSE))

occ_data <- occ_data %>%
  group_by(AggDuplicateGroupID) %>%
  mutate(
    completeness_score = rowSums(
      cbind(
        !is.na(geodeticDatum) & geodeticDatum != "",
        !is.na(decimalLatitude) & decimalLatitude != "",
        !is.na(decimalLongitude) & decimalLongitude != "",
        !is.na(coordinateUncertaintyInMeters) & coordinateUncertaintyInMeters != "",
        !is.na(year) & year != "",
        !is.na(month) & month != "",
        !is.na(day) & day != ""
      )
    ),
    AggDuplicateRank = ifelse(amongAggDuplicate == TRUE, dense_rank(desc(completeness_score)), NA) # rank with 1 being highest priority
  ) %>%
  ungroup() %>% 
  # do the same for specimen duplicates
  group_by(specimenDuplicateGroupID) %>%
  mutate(
    completeness_score = rowSums(
      cbind(
        !is.na(geodeticDatum) & geodeticDatum != "",
        !is.na(decimalLatitude) & decimalLatitude != "",
        !is.na(decimalLongitude) & decimalLongitude != "",
        !is.na(coordinateUncertaintyInMeters) & coordinateUncertaintyInMeters != "",
        !is.na(year) & year != "",
        !is.na(month) & month != "",
        !is.na(day) & day != ""
      )
    ),
    specimenDuplicateRank = ifelse(specimenDuplicate == TRUE, dense_rank(desc(completeness_score)), NA) # rank with 1 being highest priority
  ) %>%
  ungroup() 

# Select and clean up df
occ_data <- occ_data %>% 
  select(-AggDuplicateGroupID, -specimenDuplicateGroupID, -occID_lower, -completeness_score)
fwrite(occ_data, paste0("/blue/guralnick/millerjared/BoCP/data/processed/dupe-flagged/", accepted_name_file_style, ".csv"))
