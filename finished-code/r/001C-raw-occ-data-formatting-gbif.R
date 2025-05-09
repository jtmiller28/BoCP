# 001B Occurrence Data Formatting (GBIF)
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

# load in verbatim data snapshot 
dw_raw <- fread(paste0("/home/millerjared/blue_guralnick/millerjared/BoCP/data/raw/gbif-tracheophyte-dw/verbatim", task_id, ".txt"))


# SELECT for fields on interest
dw_raw <- dw_raw %>% 
  select(
    gbifID,
    occurrenceID,
    scientificName,
    year,
    month,
    day, 
    institutionCode, 
    recordedBy, 
    locality, 
    decimalLatitude,
    decimalLongitude, 
    geodeticDatum, 
    coordinateUncertaintyInMeters, 
    informationWithheld, 
    genus, 
    specificEpithet, 
    infraspecificEpithet
  ) %>% 
  mutate(verbatimScientificName = ifelse(scientificName != "" & !is.na(scientificName), scientificName, trimws(paste(genus, specificEpithet, infraspecificEpithet)))) %>% 
  select(-genus, -specificEpithet, -infraspecificEpithet) %>% 
  rename(
    verbatimYear = year, 
    verbatimMonth = month,
    verbatimDay = day, 
    verbatimInstitutionCode = institutionCode, 
    verbatimRecordedBy = recordedBy,
    verbatimDecimalLatitude = decimalLatitude, 
    verbatimDecimalLongitude = decimalLongitude, 
    verbatimGeodeticDatum = geodeticDatum, 
    verbatimCoordinateUncertaintyInMeters = coordinateUncertaintyInMeters, 
    verbatimInformationWithheld = informationWithheld
  )

fwrite(dw_raw, paste0("/home/millerjared/blue_guralnick/millerjared/BoCP/data/raw/gbif-tracheophyte-dw/select-verbatim", task_id, ".txt"))