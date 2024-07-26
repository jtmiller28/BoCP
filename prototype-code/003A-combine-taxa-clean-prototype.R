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
  if(file.exists(paste0("./data/idigbio-dws/raw/", accepted_name_filestyle_v[i], ".csv"))){
    idigbio_occur <- fread(paste0("./data/idigbio-dws/raw/", accepted_name_filestyle_v[i], ".csv"))
  } else{
    print(paste("no idigbio data for", accepted_name_v[[i]]))
  } # end of idigbio if-else
  if(file.exists(paste0("./data/gbif-dws/raw/", accepted_name_filestyle_v[i], ".csv"))){
    gbif_occur <- fread(paste0("./data/gbif-dws/raw/", accepted_name_filestyle_v[i], ".csv"))
  } else{
    print(paste("no gbif data for", accepted_name_v[[i]]))
  } # end of gbif if-else
  if(file.exists(paste0("./data/symbiota-dws/raw/", accepted_name_filestyle_v[i], ".csv"))){
    symbiota_occur <- fread(paste0("./data/symbiota-dws/raw/", accepted_name_filestyle_v[i], ".csv"))
  } else{
    print(paste("no symbiota data for", accepted_name_v[[i]]))
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
  name_relations <- name_alignment[acceptedNameParent == accepted_name_v[i]] # change to i
  raw_occur_data <- raw_occur_data %>% mutate(uuid = 1:n()) # add a unique identifier field. 
  occurrence_sn_parsed <- rgnparser::gn_parse_tidy(raw_occur_data$scientificName) # parse scientific name into relevant fields
  occurrence_sn_parsed <- occurrence_sn_parsed %>% select(canonicalfull, authorship) %>% mutate(uuid = 1:n())
  
  if(nrow(raw_occur_data) == nrow(occurrence_sn_parsed)){ # conditional to assure we're combining like things. 
    occurrence_parsed_df <- merge(raw_occur_data, occurrence_sn_parsed, by = "uuid")
    occurrence_parsed_df <- occurrence_parsed_df %>% 
      rename(scientificNameParsed = canonicalfull) %>% 
      mutate(authorshipParsed = ifelse(aggregator == "Symbiota", scientificNameAuthorship, authorship)) # conditional as symbiota already does an excellent job of parsing their authorship.
  } else {
    return(NULL) # this would be problematic, and things would not line up. 
  }
  
  accepted_name_relation <- name_relations %>% filter(acceptedName == acceptedNameParent)
  accepted_name <- unique(accepted_name_relation$acceptedName)
  possible_synonyms <- unique(name_relations$synonym)
  all_possible_names <- c(accepted_name, possible_synonyms)
  # An initial filtering of all names possible according to our taxonomic alignment. 
  occurrence_name_filtered <- occurrence_parsed_df[scientificNameParsed %in% all_possible_names]
  # A more complex filter, checking for multiple mapping boolean values as a conditional for requiring more stringent filter
  occurrence_only_acceptedName_data <- occurrence_name_filtered[scientificNameParsed == accepted_name] 
  occurrence_synonym_data <- occurrence_name_filtered[scientificNameParsed != accepted_name]
  ## If any circumstances exist where the acceptedName can multiple map, we need to strictly filter for appropriate authorship
  if(any(name_relations$acceptedNameMultMapsPossible) == TRUE){
    print("Accepted Name Multiple Mappings Possible, Proceed to Authority Matching")
    # First standardize the authorship within the occurrence data. 
    accepted_name_relation_less <- accepted_name_relation %>% filter(acceptedNameScientificNameAuthorship != "") # slight bug in name alignment seems to have made thing on occasion...
    occurrence_only_acceptedName_data <- occurrence_only_acceptedName_data %>% 
      mutate(authorshipParsed = tolower(authorshipParsed)) %>% 
      mutate(authorshipParsed = gsub(" ", "", authorshipParsed))
    acceptedAuthorship <- unique(accepted_name_relation_less$acceptedNameScientificNameAuthorship)
    filtered_acceptedName_occur_data <- occurrence_only_acceptedName_data[authorshipParsed == acceptedAuthorship]
    # And create the variable that will be used for storage conditionally...
    occurrence_acceptedName_checked <- filtered_acceptedName_occur_data
  } else{
    print("Accepted Name Does Not Have Multiple Mappings Possible, Skip Authority Matching.")
    occurrence_acceptedName_checked <- occurrence_only_acceptedName_data
  }
  
  ## Now deal with synonym matching conditional 
  viable_synonyms <- unique(name_relations$synonym) # note these encompass the acceptable name subspecific variants. 
  occur_data_synonyms <- unique(occurrence_synonym_data$scientificNameParsed) # possibly synonyms
  inner_synonyms <- intersect(viable_synonyms, occur_data_synonyms) # matching our alignment synonyms in the occurrence data.
  synonym_name_relations <- name_relations %>% filter(synonym %in% inner_synonyms) # extract relevant synonym relations
  
  synonym_df_holder <- data.table()
  for(j in 1:length(inner_synonyms)){ # move through and check each synonym case-by-case. 
    synonym_relation <- synonym_name_relations[synonym == inner_synonyms[j]]
    occur_single_synonym_data <- NULL # intialize a NULL value
    occur_single_synonym_data <- occurrence_synonym_data[scientificNameParsed == inner_synonyms[j]] # grab out the synonym
    if(unique(synonym_relation$synonymMultMapPossible) == TRUE){
      print("Synonym Name Multiple Mappings Possible, Proceed to Synonym Name Authority Matching")
      occur_single_synonym_data  <- occur_single_synonym_data %>% 
        mutate(authorshipParsed = tolower(authorshipParsed)) %>% 
        mutate(authorshipParsed = gsub(" ", "", authorshipParsed))
      synonymAuthorship <- unique(synonym_relation$synonymNameSpacelessWCVPAuthorship)
      filtered_occurrence_synonym_data <- occur_single_synonym_data[authorshipParsed == synonymAuthorship]
      # and write out the filtered data
      occurrence_synonym_checked <- filtered_occurrence_synonym_data
      synonym_df_holder <- rbind(synonym_df_holder, occurrence_synonym_checked, fill=TRUE)
    } else {
      print("Synonym Name Does Not Have Multiple Mappings Possible, No Filtering By Authority Necessary")
      synonym_df_holder <- rbind(synonym_df_holder, occur_single_synonym_data, fill=TRUE)
    }
  }

occurrence_taxonomy_cleaned <- rbind(occurrence_acceptedName_checked, synonym_df_holder)
  # Summary (built at the end!)
  taxa_clean_info_tb <- data.frame(
    acceptedNameParent = accepted_name_v[i],
    totalAggregatorRecords = nrow(raw_occur_data),
    initialMatchingClean = nrow(occurrence_name_filtered), 
    numOfAcceptedNameRecords = nrow(occurrence_only_acceptedName_data), 
    AcceptedNameMultipleMatchingClean = nrow(occurrence_acceptedName_checked), 
    numOfSynonymRecords = nrow(occurrence_synonym_data),
    SynonymNameMultipleMatchingClean = nrow(synonym_df_holder), 
    totalTaxonomyCleanedRecords = nrow(occurrence_taxonomy_cleaned)
  )
  fwrite(taxa_clean_info_tb, file = "./data/processed/cleaning-summaries/taxa_clean_info_tb.csv", append = TRUE, col.names = !file.exists("./data/processed/cleaning-summaries/taxa_clean_info_tb.csv"))
  fwrite(occurrence_taxonomy_cleaned, file = paste0("./data/processed/taxa-cleaned-data/", accepted_name_filestyle_v[i], ".csv"))
} # end of i loop
