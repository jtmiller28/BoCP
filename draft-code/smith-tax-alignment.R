### Title: Taxon Alignment for stephen's Smith & Brown phylogeny
### Author: JT Miller
### Date: 02-24-2025

# Load packages
library(data.table)
library(rgnparser)
library(tidyverse)
library(taxadb)

#### Set up pathing to use gnparser ti access
my_path <- Sys.getenv("PATH") # grab our path
Sys.setenv(PATH = paste0(my_path, "/home/millerjared/gnparser"))

#### Load taxonomic backbone
wcvp_backbone <- fread("/home/millerjared/blue_guralnick/millerjared/BoCP/data/raw/wcvp_names.csv")
wcvp_backbone_s <- select(wcvp_backbone, # select relevant fields within the wcvp backbone for alignment
                          wcvpPlantNameID = plant_name_id, 
                          wcvpFamily = family, 
                          wcvpName = taxon_name, 
                          wcvpNameAuthors = taxon_authors, 
                          wcvpNameStatus = taxon_status,
                          wcvpAcceptedPlantNameID = accepted_plant_name_id, 
                          )
#### Load smith & brown name list to be reconciled against wcvp taxonomy
# note that this will generate a warning, presumably due to \n issues at the end of the name list, I manaully checked the head and tail and compared against the txt file's contents so we should be okay
smith_names <- readLines("/home/millerjared/blue_guralnick/millerjared/BoCP/data/raw/global_tree_search_trees_by_country_updated(1).tree_list")
smith_names_fix <- smith_names[-1] # header read issues, remove the first item
smith_names_fix <- gsub('[\'"]', '', smith_names_fix)

#### Use parsing to separate out the name & authors, then create a dataframe to align against, remove some troublesome cases that will cause issues with alignment
smith_names_parsed <- gn_parse_tidy(smith_names_fix)
nrow(smith_names_parsed) == length(smith_names_fix) # if TRUE, no need to dedup
smith_names_parsed_df <- smith_names_parsed %>% # retain only important info
  select(smithVerbatimName = verbatim,
         smithName = canonicalfull, 
         smithAuthor = authorship, 
         numNameComponents = cardinality) %>% 
  # Specific to this name list, I've identified that the cases that have relevant subspecific info will cause issues when matching. 
  mutate(smithName = ifelse(numNameComponents > 2, word(smithName, 1,2), smithName)) %>% 
  filter(!numNameComponents <= 1) # remove these names, as they will not be resolvable

#### Merge smith names with wcvp backbone
# merge our name table on the wcvp backbone using the smith verbatim names
matched_names <- merge(smith_names_parsed_df, wcvp_backbone_s, by.x = "smithName", by.y = "wcvpName", all.x = TRUE)
matched_names <- distinct(matched_names)

#### Identify multiple name mappers, check if removal can proceed. 
wcvp_backbone_multmaps <- wcvp_backbone_s %>% 
  group_by(wcvpName) %>% 
  mutate(numDifferentWCVPNamePaths = n_distinct(wcvpAcceptedPlantNameID)) %>% 
  filter(numDifferentWCVPNamePaths > 1)

matched_names <- matched_names %>% 
  mutate(nameMultMapStatus = ifelse(smithName %in% wcvp_backbone_multmaps$wcvpName, TRUE, FALSE))

# large num of occurrences, attempt recovery of some of these names using exact matching w/author 
# rm spaces & capitilization, then perform an exact match 
matched_names <- matched_names %>% 
  mutate(spacelessSmithAuthor = tolower(gsub(" ", "", smithAuthor)),
         spacelessWCVPNameAuthors = tolower(gsub(" ", "", wcvpNameAuthors)))

matched_names <- matched_names %>% 
  mutate(multMapAuthorMatch = ifelse(nameMultMapStatus == TRUE, spacelessSmithAuthor == spacelessWCVPNameAuthors, NA))

#### Remove mult mapping names unless they can be mapped directly with author
matched_names <- matched_names %>% 
  filter(nameMultMapStatus == FALSE | multMapAuthorMatch == TRUE)

#### Setup accepted name mappings
wcvp_accepted <- wcvp_backbone_s[wcvpNameStatus == "Accepted"]

#### Merge and denote aligned name
aligned_data <- matched_names %>% 
  left_join(select(wcvp_accepted, wcvpAcceptedPlantNameID, wcvpAlignedName = wcvpName, wcvpAlignedNameStatus = wcvpNameStatus)) %>% 
  distinct() # ensure distinct

# select for simplicity
aligned_data_s <- aligned_data %>% 
  select(verbatimName = smithVerbatimName,
         name = smithName, 
         nameAuthors = smithAuthor, 
         nameStatusInWCVP = wcvpNameStatus, 
         wcvpAlignedName, 
         wcvpAlignedNameStatus, 
         wcvpAcceptedPlantNameID) %>% # run distinct, we dont care about mult authors as the names with mult mappings have been removed. 
  distinct(verbatimName, name, nameStatusInWCVP, wcvpAlignedName, wcvpAlignedNameStatus, wcvpAcceptedPlantNameID)
  
# Denote what happens to unplaced taxa, either when they themselves are unplaced OR if they are a synonym/somehow else related to an unplaced name 
aligned_data_s <- aligned_data_s %>% 
  mutate(wcvpAlignedNameStatus = case_when(
    !is.na(wcvpAlignedName) ~ wcvpAlignedNameStatus, 
    is.na(wcvpAlignedName) & !is.na(nameStatusInWCVP) ~ "Unplaced", 
    is.na(wcvpAlignedName) & is.na(nameStatusInWCVP) ~ NA
  )) %>% 
  mutate(foundInWCVP = ifelse(!is.na(wcvpAlignedNameStatus), TRUE, FALSE))
 


# Take these data and merge with stephen's intial list, allowing for any taxa that failed to align to be shown 
smith_names_reconciled <- smith_names_parsed_df %>% 
  select(verbatimName = smithVerbatimName) %>% 
  left_join(aligned_data_s) %>% 
  distinct() %>% 
  group_by(verbatimName) %>% 
  mutate(n = n()) %>% 
  ungroup()

# deal with these names that have multiple statuses, this is due to authors
smith_names_mult_stats <- filter(smith_names_reconciled, n > 1)
smith_names_mult_fixed <- filter(smith_names_mult_stats, nameStatusInWCVP %in% c("Accepted", "Synonym"))
smith_names_mult_fixed <- smith_names_mult_fixed %>% 
  group_by(name) %>% 
  mutate(shouldTakeAcceptedStatus = ifelse(
    any(nameStatusInWCVP == "Accepted") & any(nameStatusInWCVP == "Synonym"),
    TRUE, 
    FALSE
  )) %>% ungroup() %>% 
  filter((shouldTakeAcceptedStatus & nameStatusInWCVP == "Accepted")| !shouldTakeAcceptedStatus) %>% 
  select(-shouldTakeAcceptedStatus)
n_distinct(smith_names_mult_stats$name) == n_distinct(smith_names_mult_fixed$name) # works!

# remove and readd based on fix
smith_names_reconciled <- smith_names_reconciled %>% 
  filter(n == 1) %>% 
  rbind(smith_names_mult_fixed) %>% 
  select(-n)

# identify cases where there are both Accepted, and Synonym

check <- setdiff(smith_names_parsed_df$smithVerbatimName, smith_names_reconciled$verbatimName)

# write out
fwrite(smith_names_reconciled, "/home/millerjared/blue_guralnick/millerjared/BoCP/outputs/global-tree-wcvp-name-alignment.csv")
