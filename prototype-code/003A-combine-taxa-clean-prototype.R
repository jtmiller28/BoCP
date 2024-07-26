# 003A-combine-taxa-clean-prototype
# Author: JT Miller 
# Date: 07-24-2024
# Project: BoCP 

## This script is designed to pull all data together from gbif, idigbio, and symbiota downloads,
## and taxonomically clean it according to our name alignment

# Load libraries
library(data.table)
library(dplyr)
library(lubridate)
library(rgnparser)

# call edited gatoRs fxns
source("finished-code/r/gatoRs-fxns-edited.R")
# Set up pathing to use gnparser
my_path <- Sys.getenv("PATH") # grab our path
Sys.setenv(PATH = paste0(my_path, "/home/millerjared/gnparser"))

# Read in the name alignment
name_alignment <- fread("./data/processed/finalized_name_alignment_wcvp.csv")
# Create an accepted name vector & filestyle version.
accepted_name_v <- unique(name_alignment$acceptedNameParent)
accepted_name_filestyle_v <- gsub(" ", "-", accepted_name_v)

# Load in some specific information about basisOfRecord for symbiota
sym_basis_info <- fread("./data/symbiota-dws/symbiota/SEINet_Observations_or_Qustionable.csv")
sym_basis_info <- select(sym_basis_info, -collid, -institutionCode, -collectionCode,
                         -collectionName, -occid)
## Create a loop that first combines the data from the three download directories, then taxonomically cleans it for the name alignment 
for(i in 1:length(accepted_name_v)){
  if(file.exists(paste0("./data/idigbio-dws/raw/", accepted_name_filestyle_v[1], ".csv"))){
    idigbio_occur <- fread(paste0("./data/idigbio-dws/raw/", accepted_name_filestyle_v[1], ".csv"))
  } else{
    print(paste("no idigbio data for", accepted_name_v[[1]]))
  } # end of idigbio if-else
  if(file.exists(paste0("./data/gbif-dws/raw/", accepted_name_filestyle_v[1], ".csv"))){
    gbif_occur <- fread(paste0("./data/gbif-dws/raw/", accepted_name_filestyle_v[1], ".csv"))
  } else{
    print(paste("no gbif data for", accepted_name_v[[1]]))
  } # end of gbif if-else
  if(file.exists(paste0("./data/symbiota-dws/raw/", accepted_name_filestyle_v[1], ".csv"))){
    symbiota_occur <- fread(paste0("./data/symbiota-dws/raw/", accepted_name_filestyle_v[1], ".csv"))
  } else{
    print(paste("no symbiota data for", accepted_name_v[[1]]))
  } # end of symbiota if-else
  
  # Standardize fields between idigbio and gbif
  ## First deal with basisOfRecord, Ed provided a list of known records with non-specimen info, so join them up and all NAs will be perservedSpecimens.
  symbiota_occur <- left_join(symbiota_occur, sym_basis_info, by = "occurrenceID") # create the relational table to query in an if-else
  symbiota_occur <- symbiota_occur %>% 
    mutate(basisOfRecord = ifelse(!is.na(basisOfRecord), basisOfRecord, "PerservedSpecimen"))
  symbiota_occur <- symbiota_occur %>% 
    rename(scientificName = sciname, scientificNameAuthorship = scientificNameAuthorship, ID = occid, 
           latitude = decimalLatitude, longitude = decimalLongitude) %>% 
    select(-collectionCode, -collectionName,  -catalogNumber, -family, -taxonRemarks, -identifiedBy, 
           -dateIdentified, -identificationReferences, -identificationRemarks, -typeStatus, -recordNumber, 
           -associatedCollectors, -verbatimEventDate, -occurrenceRemarks, -associatedOccurrences, 
           -dataGeneralizations, -associatedTaxa, -dynamicProperties, -verbatimAttributes, -reproductiveCondition,
           -cultivationStatus, -establishmentMeans, -municipality, -verbatimCoordinates, -minimumElevationInMeters, 
           -maximumElevationInMeters, -verbatimElevation) %>% 
    mutate(                         # its in symbiota so this must be the case.
           aggregator = "Symbiota",
           year = lubridate::year(lubridate::ymd(eventDate)),
           month = lubridate::month(lubridate::ymd(eventDate)),
           day = lubridate::day(lubridate::ymd(eventDate)))
  symbiota_occur <- correct_class(symbiota_occur) # see gatoRs fxns edited for details on class correction. 
  idigbio_occur <- idigbio_occur %>% 
    select(-genus, -specificEpithet, -infraspecificEpithet) %>% 
    mutate(scientificNameAuthorship = NA) # this is not parsed in the parent data as of yet, but we'll use NA for organization purposes. 
  idigbio_occur <- correct_class(idigbio_occur)
  gbif_occur <- gbif_occur %>% 
    select(-genus, -specificEpithet, -infraspecificEpithet) %>% 
    mutate(scientificNameAuthorship = NA)
  gbif_occur <- correct_class(gbif_occur)
  # Combine all of the data
  raw_occur_data <- rbind(symbiota_occur, idigbio_occur, gbif_occur)
  
  ## Taxonomically clean the data, then write it to a associated dir. 
  
  
  
  # Summary (built at the end!)
  taxa_clean_info_tb <- data.frame(
    acceptedNameParent = accepted_name_v[1],
    totalAggregatorRecords = nrow(raw_occur_data)
  )
  fwrite(taxa_clean_info_tb, file = "./data/processed/cleaning-summaries/taxa_clean_info_tb.csv", append = TRUE, col.names = !file.exists("./data/processed/cleaning-summaries/taxa_clean_info_tb.csv"))
} # end of i loop
