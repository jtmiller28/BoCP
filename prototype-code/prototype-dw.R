# Prototype Download Script for iDigBio occurrence records
## Goals:
# 1. Pull idigbio occurrence data for taxa of interest
# 2. Filter down data to ensure taxonomy is congruent 
# 3. Create appropriate error handlers to deal with automation

# Load packages
library(gatoRs)
library(data.table)
library(dplyr)
library(rgnparser)
source("prototype-code/gatoRs-fxns-edited.R") # edited gatoRs fxns to allow for only pulling idigbio

# Read in the name alignment
name_alignment <- fread("/home/jtmiller/my_elements/jtmiller/BoCP/data/processed/finalized_name_alignment.csv")

# Structure these names in a list with acceptedName + synonym(s)
accepted_name_v <- unique(name_alignment$acceptedName)
name_list <- list() # initialize empty list

bench::bench_time({
for(i in 1:length(accepted_name_v)){
  p <- name_alignment[acceptedName == accepted_name_v[[i]]] # create a for-loop that goes through by the accepted name and grabs synonyms, storing them both as a vector in a list. 
  s <- p[!is.na(synonym)]
  x <- print(s$synonym)
  a <- print(unique(p$acceptedName))
  name_list[[i]] <- unique(c(a,x))
}
})


for(i in 1:length(name_list)){
# Proceed to download
bench::bench_time({
occurrence_dw <- gators_download_edited(name_list[[1]], call.idigbio = TRUE, call.gbif = FALSE)
# First Clean taxonomy for overfetched names
occurrence_dw <- as.data.table(occurrence_dw)
occurrence_parsed <- gn_parse_tidy(occurrence_dw$concatenatedName)  # break up names
occurrence_parsed_names_df <- merge(occurrence_dw, occurrence_parsed, by.x = 'concatenatedName', by.y = 'verbatim')
occurrence_dw <- occurrence_dw[scientificName %in% name_list[[1]]] # assure that what we are downloading for is the same

# Initiate cleaning for records without georeference (lat lon)

})
}


