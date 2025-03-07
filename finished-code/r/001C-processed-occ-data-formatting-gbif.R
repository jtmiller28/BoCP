# 001B Occurrence Data Formatting (GBIF) processed data
# Author: JT Miller 
# Date: 02-18-2024
# Project: BoCP 

# Load packages 
library(tidyverse)
library(data.table)
library(arrow)
library(DBI)
library(duckdb)

# set dir
setwd("/home/millerjared/blue_guralnick/millerjared/BoCP/")

# Set up array logic
start_num <- as.numeric(Sys.getenv("START_NUM"))
task_id <- as.numeric(start_num)
part <- paste0("part", task_id)

print(task_id)


# load in the processed data snapshot 
dw_processed <- fread(paste0("/home/millerjared/blue_guralnick/millerjared/BoCP/data/raw/gbif-tracheophyte-dw/occurrence", task_id, ".txt"))

# SELECT fields of interest
dw_processed <- dw_processed %>%
  select(
    gbifID, 
    occurrenceID, 
    basisOfRecord, 
    eventDate, 
    coordinateUncertaintyInMeters, 
    decimalLatitude, 
    decimalLongitude, 
    informationWithheld, 
    year, 
    month, 
    day
    ) %>% 
  # RENAME fields 
  rename(
    processedCoordinateUncertaintyInMeters = coordinateUncertaintyInMeters, 
    processedDecimalLatitude = decimalLatitude, 
    processedDecimalLonitude = decimalLongitude, 
    processedInformationWithheld = informationWithheld, 
    processedYear = year, 
    processedMonth = month, 
    processedDay = day
  )

fwrite(dw_processed, paste0("/home/millerjared/blue_guralnick/millerjared/BoCP/data/raw/gbif-tracheophyte-dw/select-occurrence", task_id, ".txt"))