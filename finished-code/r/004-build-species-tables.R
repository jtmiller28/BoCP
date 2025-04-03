# 004 Build Species Tables
# Author: JT Miller 
# Date: 03-06-2024
# Project: BoCP 

# load packages
library(data.table)
library(arrow)
library(duckdb)
library(DBI)
library(dbplyr)
library(rgnparser)
library(tidyverse)

#### Set up pathing to use gnparser
my_path <- Sys.getenv("PATH") # grab our path
Sys.setenv(PATH = paste0(my_path, "/home/millerjared/gnparser"))

# Set up array logic
start_num <- as.numeric(Sys.getenv("START_NUM"))
task_id <- as.numeric(start_num)
part <- paste0("part", task_id)

# read the taxonomic harmonized list of names + synonyms
name_list <- readRDS("/blue/guralnick/millerjared/BoCP/data/processed/name_list.rds")

# check for names that have already been retrieved
list_found_names <- list.files("/blue/guralnick/millerjared/BoCP/data/processed/species-occs/")
list_found_names <- gsub("-", " ", list_found_names)
list_found_names <- gsub(".csv", "", list_found_names)
# retrieve accepted names 
accepted_name_retrieval_list <- list()
for(i in 1:length(name_list)){
accepted_names_retrieval <- name_list[[i]][1]
accepted_name_retrieval_list[[i]] <- accepted_names_retrieval
}
accepted_name_retrieval_vector <- do.call(rbind, accepted_name_retrieval_list)
# check diff between dir names and total names to call
names_not_found <- setdiff(accepted_name_retrieval_vector, list_found_names)
# extract 
# Extract elements where any name in the list matches a target name
name_list_update <- name_list[sapply(name_list, function(x) any(x %in% names_not_found))]
# use task id as an index to grab out the set of names, format for SQL LIKE query, format taxon name for file storage. 
names_to_retrieve <- name_list_update[[task_id]]
names_to_retrieve_query <- toupper(names_to_retrieve)
names_to_retrieve_query <- gsub(" ", "%", names_to_retrieve_query)
accepted_name <- names_to_retrieve[1]
accepted_name_file_style <- gsub(" ", "-", accepted_name)

# print to track
print(paste("retrieving", accepted_name, "occurrence data"))
# Connect to DuckDB
print("opening db connection")
connect <- dbConnect(duckdb::duckdb())

# Set nulls
idigbio_df <- NULL
gbif_df <- NULL
symbiota_df <- NULL

# Set up directory paths for each aggregator's parquet files
idigbio_file_path <- "/blue/guralnick/millerjared/BoCP/data/raw/parquet-occs/idigio_full_occ.parquet"
gbif_file_path <- "/blue/guralnick/millerjared/BoCP/data/raw/parquet-occs/gbif_full_occ.parquet"
symbiota_file_path <- "/blue/guralnick/millerjared/BoCP/data/raw/parquet-occs/symbiota_occ.parquet"

# Construct the SQL query that grabs the accepted name & all of its valid synonyms
query_idigbio <- paste0("SELECT * FROM '", idigbio_file_path, "' WHERE UPPER(verbatimScientificName) LIKE ", 
                        paste0("'%", names_to_retrieve_query, "%'", collapse = " OR UPPER(verbatimScientificName) LIKE "))
query_gbif <- paste0("SELECT * FROM '", gbif_file_path, "' WHERE UPPER(verbatimScientificName) LIKE ", 
                     paste0("'%", names_to_retrieve_query, "%'", collapse = " OR UPPER(verbatimScientificName) LIKE "))
query_symbiota <- paste0("SELECT * FROM '", symbiota_file_path, "' WHERE UPPER(verbatimScientificName) LIKE ", 
                         paste0("'%", names_to_retrieve_query, "%'", collapse = " OR UPPER(verbatimScientificName) LIKE "))
print("executing queries")
# Execute queries on the databases
idigbio_df <- dbGetQuery(connect, query_idigbio)
gbif_df <- dbGetQuery(connect, query_gbif)
symbiota_df <- dbGetQuery(connect, query_symbiota)

# Taxonomically Clean the data for exact matches with our names 
print("idigbio parsing name, and taxonomic check")
## idigbio
if(nrow(idigbio_df) > 0){
idigbio_df_names_parsed <- gn_parse_tidy(na.omit(idigbio_df$verbatimScientificName))
idigbio_df_names_parsed <- idigbio_df_names_parsed %>% mutate(n = 1:n()) %>% select(canonicalfull, n)
idigbio_df <- idigbio_df %>% mutate(n = 1:n())
idigbio_df <- merge(idigbio_df, idigbio_df_names_parsed, by = "n")
idigbio_df <- idigbio_df %>% mutate(taxonomicExactMatch = ifelse(toupper(canonicalfull) %in% toupper(names_to_retrieve), TRUE, FALSE)) 
}
print("gbif parsing name, and taxonomic check")
## gbif
if(nrow(gbif_df) > 0){
gbif_df_names_parsed <- gn_parse_tidy(na.omit(gbif_df$verbatimScientificName)) # note, there are extreme cases where parsing fails causing memory leak error. 
gbif_df_names_parsed <- gbif_df_names_parsed %>% mutate(n = 1:n()) %>% select(canonicalfull, n)
gbif_df <- gbif_df %>% mutate(n = 1:n())
gbif_df <- merge(gbif_df, gbif_df_names_parsed, by = "n")
gbif_df <- gbif_df %>% mutate(taxonomicExactMatch = ifelse(toupper(canonicalfull) %in% toupper(names_to_retrieve), TRUE, FALSE)) 
}
print("symbiota parsing name, and taxonomic check")
## symbiota
if(nrow(symbiota_df) > 0){
symbiota_df_names_parsed <- gn_parse_tidy(na.omit(symbiota_df$verbatimScientificName))
symbiota_df_names_parsed <- symbiota_df_names_parsed %>% mutate(n = 1:n()) %>% select(canonicalfull, n)
symbiota_df <- symbiota_df %>% mutate(n = 1:n())
symbiota_df <- merge(symbiota_df, symbiota_df_names_parsed, by = "n")
symbiota_df <- symbiota_df %>% mutate(taxonomicExactMatch = ifelse(toupper(canonicalfull) %in% toupper(names_to_retrieve), TRUE, FALSE)) 
}
print("formatting tables")
# Standardize formatting, deal with verbatim vs processed fields 
if(nrow(idigbio_df) > 0){
## idigbio
idigbio_df <- idigbio_df %>% 
  rename(parsedName = canonicalfull, 
         year = verbatimYear,
         month = verbatimMonth, 
         day = verbatimDay, 
         uuid = coreid, 
         informationWithheld = verbatimInformationWithheld,
         decimalLatitude = verbatimDecimalLatitude, 
         decimalLongitude = verbatimDecimalLongitude, 
         geodeticDatum = verbatimGeodeticDatum,
         institutionCode = verbatimInstitutionCode, 
         recordedBy = verbatimRecordedBy, 
         locality = verbatimLocality
         ) %>% 
  mutate(species = accepted_name, 
         aggregator = "idigbio") %>% 
  select(uuid,
         aggregator,
         occurrenceID,
         basisOfRecord, 
         species, 
         verbatimScientificName, 
         parsedName, 
         geodeticDatum, 
         decimalLatitude, 
         decimalLongitude,
         coordinateUncertaintyInMeters, 
         eventDate, 
         year, 
         month, 
         day, 
         informationWithheld, 
         recordedBy, 
         institutionCode, 
         locality,
         taxonomicExactMatch
         )
}
## gbif 
if(nrow(gbif_df) > 0){
gbif_df <- gbif_df %>% 
  rename(parsedName = canonicalfull, 
         uuid = gbifID
  ) %>% 
  mutate(species = accepted_name, 
         aggregator = "gbif") %>% 
  select(uuid,
         aggregator,
         occurrenceID,
         basisOfRecord, 
         species, 
         verbatimScientificName, 
         parsedName, 
         geodeticDatum, 
         decimalLatitude, 
         decimalLongitude,
         coordinateUncertaintyInMeters, 
         eventDate, 
         year, 
         month, 
         day, 
         informationWithheld, 
         recordedBy, 
         institutionCode, 
         locality,
         taxonomicExactMatch
  )
}
## symbiota 
if(nrow(symbiota_df) > 0){
symbiota_df <- symbiota_df %>% 
  rename(parsedName = canonicalfull) %>% 
  mutate(species = accepted_name, 
         aggregator = "symbiota", 
         uuid = NA) %>% 
  select(uuid,
         aggregator,
         occurrenceID,
         basisOfRecord, 
         species, 
         verbatimScientificName, 
         parsedName, 
         geodeticDatum, 
         decimalLatitude, 
         decimalLongitude,
         coordinateUncertaintyInMeters, 
         eventDate, 
         year, 
         month, 
         day, 
         informationWithheld, 
         recordedBy, 
         institutionCode, 
         locality,
         taxonomicExactMatch
  )
}
print("compile datasets")
# Create a data list 
data_list <- list(symbiota_df, idigbio_df, gbif_df)
# Filter out anything that contains nothing.
non_null_data_list <- Filter(Negate(is.null), data_list)
# Combine all non null dataframes
if(length(non_null_data_list) > 0){
  raw_occur_data <- do.call(rbind, non_null_data_list)
} else { 
  raw_occur_data <- NULL
}
if(is.null(raw_occur_data)){
  empty_df <- data.frame(
    acceptedParentName = accepted_name
  )
  fwrite(empty_df, file = "/home/millerjared/blue_guralnick/millerjared/BoCP/outputs/names-without-occ-data.csv", 
         append = TRUE, col.names = !file.exists("/home/millerjared/blue_guralnick/millerjared/BoCP/outputs/names-without-occ-data.csv"))
  
}
print("writing data")
# write out the data
fwrite(raw_occur_data, paste0("/blue/guralnick/millerjared/BoCP/data/processed/species-occs/", accepted_name_file_style, ".csv"))
print("closing connection")
# Close the connection
dbDisconnect(connect)
