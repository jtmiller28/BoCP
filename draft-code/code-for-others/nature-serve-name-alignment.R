### Title: Nature Serve Name Alignment
### Project: BoCP
### Date: 04/01/2025
### Author: JT Miller

### Aligning Israel's nature serve taxonomy so that those names are usable with BoCP data

# load packages 
library(data.table)
library(tidyverse)

# load nature serve taxonomy 
ns_names <- fread("/blue/guralnick/millerjared/BoCP/data/raw/natureserve_final_list(Sheet1).csv")

# from looking at the user_supplied name, it appears it has already been parsed for authorship, so skip gnparsing 
ns_names <- ns_names %>% select(user_supplied_name) # select only field of interest

# load in BoCP's taxonomic harmonized table
bocp_names <- fread("/blue/guralnick/millerjared/BoCP/data/processed/wcvp-ncbi-alignment-na.csv")

# we're going to try merging the the user supplied name upon the name field within bocp's harmonized taxons
merged_names <- merge(ns_names, bocp_names, by.x = "user_supplied_name", by.y = "name", all.x = TRUE)

merged_names <- merged_names %>% 
  select(user_supplied_name, nameStatus, alignedName, alignedParentName, ncbiPossibleIDs, source) %>% 
  distinct()

failed_names <- merged_names %>% 
  filter(is.na(nameStatus))

# filter out these fails, write them out in case Israel needs them. 
ns_names_aligned <- merged_names %>% 
  filter(!is.na(nameStatus)) # remove records that did not find taxonomic harm 

fwrite(ns_names_aligned, "/blue/guralnick/millerjared/BoCP/data/processed/nature-serve-names-aligned.csv")
fwrite(failed_names, "/blue/guralnick/millerjared/BoCP/data/processed/nature-serve-alignFail-names.csv")


## Adding Canada names 
canada1 <- fread("/blue/guralnick/millerjared/BoCP/data/raw/canada1_final_list.csv")
canada2 <- fread("/blue/guralnick/millerjared/BoCP/data/raw/canada2_final_list.csv")

canada1 <- select(canada1, user_supplied_name)
canada2 <- select(canada2, user_supplied_name)
canada_names <- rbind(canada1, canada2)
canada_names <- canada_names %>% distinct()

# align 
# we're going to try merging the the user supplied name upon the name field within bocp's harmonized taxons
merged_names <- merge(canada_names, bocp_names, by.x = "user_supplied_name", by.y = "name", all.x = TRUE)

merged_names <- merged_names %>% distinct()

failed_canada_names <- merged_names %>% 
  filter(is.na(nameStatus))

ca_aligned_names <- merged_names %>% 
  filter(!is.na(nameStatus))

fwrite(ca_aligned_names, "/blue/guralnick/millerjared/BoCP/data/processed/canada-names-aligned.csv")
fwrite(failed_canada_names, "/blue/guralnick/millerjared/BoCP/data/processed/canada-names-failedAlignment.csv")
