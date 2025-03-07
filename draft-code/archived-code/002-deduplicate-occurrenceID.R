# 002 Deduplicate Occurrence ID 
# Author: JT Miller 
# Date: 02-28-2024
# Project: BoCP 

## Purpose: To join the raw and processed datasets of idigbio and gbif, its required that we use the occurrenceID as a joining parameter.
## Issue is that occurrenceID is commonly duplicated, presumably to denote organismal duplicates or data entry duplicates, more concerning would be misapplications of occurrenceIDs however, as these would cause joins to manufacture data 

# Load libraries 
library(tidyverse)
library(duckdb)
library(data.table)
library(dplyr)
library(dbplyr)
## dbplyr way for R clarity 
# Load the parquet file as a table reference
con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
# register parquet file as a table 
dbExecute(con, "INSTALL parquet; LOAD parquet;")
dbExecute(con, "CREATE TABLE idigbio_raw_occ AS SELECT * FROM read_parquet('/home/millerjared/blue_guralnick/millerjared/BoCP/data/raw/parquet-occs/idigbio_raw_occ.parquet');")

# run con as the database connection, and lazy eval idigbio_raw_occ table
idigbio_raw_occ <- tbl(con, "idigbio_raw_occ")

filtered_records <- idigbio_raw_occ %>%
  left_join( # join the results of the following to the full table
    idigbio_raw_occ %>%
      group_by(occurrenceID) %>% 
      summarize(
        occ_count = n(), # count the num of times occurrenceID exists in the table
        unique_names = n_distinct(verbatimScientificName), # count the number of distinct verbatimScientificNames 
        .groups = "drop" # remove the group by for the resulting table
      ),
    by = "occurrenceID" # join key
  ) %>%
  mutate(
    dupID = case_when( # create a new field that identifies what type of duplicate we're dealing with
      occ_count == 1 ~ NA_character_, # if a occurrenceID is uniq to the tbl, then NA 
      unique_names == 1 ~ "organismalDup", # if its non-uniq, but only has one unique name assocaited with that occurrenceID then its an organismal duplicate
      TRUE ~ "misAppliedDup" # if its non-uniq, and has multiple names associated with it its misapplied 
    ),
    # do a best case count, where we score records higher based upon how many relevant fields are populated
    non_null_count = 
      (case_when(!is.na(verbatimDecimalLatitude) & verbatimDecimalLatitude != "" ~ 1, TRUE ~ 0) +
         case_when(!is.na(verbatimDecimalLongitude) & verbatimDecimalLongitude != "" ~ 1, TRUE ~ 0) +
         case_when(!is.na(verbatimYear) & verbatimYear != "" ~ 1, TRUE ~ 0) +
         case_when(!is.na(verbatimMonth) & verbatimMonth != "" ~ 1, TRUE ~ 0) +
         case_when(!is.na(verbatimDay) & verbatimDay != "" ~ 1, TRUE ~ 0) +
         case_when(!is.na(verbatimGeodeticDatum) & verbatimGeodeticDatum != "" ~ 1, TRUE ~ 0) +
         case_when(!is.na(verbatimCoordinateUncertaintyInMeters) & verbatimCoordinateUncertaintyInMeters != "" ~ 1, TRUE ~ 0) +
         case_when(!is.na(verbatimInformationWithheld) & verbatimInformationWithheld != "" ~ 1, TRUE ~ 0))
  ) %>%
  # arrange the data in order to of these priority fields grouped by occurrenceID
  arrange(occurrenceID, desc(non_null_count)) %>%
  group_by(occurrenceID, dupID) %>%
  mutate(rank = row_number()) %>% # rank them
  ungroup() %>%
  # filter these records, if not a duplicate always keep, if organismal duplicate take the best ranked record, if a misapplied duplicate always remove
  filter(
    is.na(dupID) | 
      (dupID == "organismalDup" & rank == 1) | 
      (!(dupID %in% "misAppliedDup") & rank == 1)
  )

# Collect results into an R dataframe
filtered_records <- filtered_records %>% collect()

## Now link the processed data to the correct occurrenceID, removing duplicates. 
dbExecute(con, "CREATE TABLE idigbio_processed_occ AS SELECT * FROM read_parquet('/home/millerjared/blue_guralnick/millerjared/BoCP/data/raw/parquet-occs/idigbio_processed_occ.parquet');")

# run con as the database connection, and lazy eval idigbio_processed_occ table
idigbio_processed_occ <- tbl(con, "idigbio_processed_occ")

# run filter, note that its not necessary to deal with misApplied records as these occurrence IDs have been removed (and therefore will be removed by the INNER JOIN of the last step to combine these data)
# this filter is simple, any record with non-uniq occurrenceIDs will be ranked by number of fields filled and filtered via that method. this shouldnt matter because its an organismal duplicate (same organism)

idigbio_processed_occ <- idigbio_processed_occ %>%
  mutate(
    coordinateUncertaintyInMeters = as.numeric(na_if(coordinateUncertaintyInMeters, "")), 
    eventDate = na_if(eventDate, ""), 
    basisOfRecord = na_if(basisOfRecord, "")
  )


filtered_processed_records <- idigbio_processed_occ %>% 
  mutate(coordinateUncertaintyInMeters = as.numeric(na_if(coordinateUncertaintyInMeters, ""))) %>%  # Convert empty strings to NA and ensure numeric type
  left_join(
    idigbio_processed_occ %>% 
      group_by(occurrenceID) %>% 
      summarize(occ_count = n(), # num of occurrenceID instances in the table
                .groups = "drop"), 
    by = "occurrenceID" # join key
  ) %>% 
  mutate(
    non_null_count = 
      case_when(!is.na(basisOfRecord) & basisOfRecord != "" ~ 1, TRUE ~ 0) +
      case_when(!is.na(eventDate) & eventDate != "" ~ 1, TRUE ~ 0) +
      case_when(!is.na(coordinateUncertaintyInMeters) & coordinateUncertaintyInMeters != "" ~ 1, TRUE ~ 0)
  ) %>% 
  arrange(occurrenceID, desc(non_null_count)) %>% 
  group_by(occurrenceID, occ_count) %>% 
  mutate(rank = row_number())
  
filtered_processed_records <- filtered_processed_records %>% collect()
