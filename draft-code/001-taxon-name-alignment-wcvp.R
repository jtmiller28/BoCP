# 001-taxon-name-alignment wcvp rewrite
# Author: JT Miller 
# Date: 06-19-2024
# Project: BoCP 

## Load Libraries
library(data.table)
library(tidyverse)
library(rgnparser)
library(taxadb)

## Load Data
raw_bocp_names <- fread("/home/jt-miller/Gurlab/BoCP/data/raw/USACANMEX_plantNamesFieldsAdded.csv",encoding = "UTF-8")
wcvp_backbone <- fread("/home/jt-miller/Gurlab/BoCP/data/raw/taxonomic_backbones/wcvp/wcvp_names.csv")

### Data Prep
## Concatenate Authorship, then parse out using gnparser
raw_bocp_names$concatenatedName <- apply(raw_bocp_names[, c("Name", "ExteriorNameAuthorship", "ExteriorNameAuthorship2", "V5")], 1, function(x) paste(x[x != ""], collapse = ", ")) # Collapse fields into one by adding a space when there is a string present
bocp_names <- raw_bocp_names[, .(TropicosNameID = NameID, concatenatedName)] # Select only fields of interest 
bocp_parsed_names <- gn_parse_tidy(bocp_names$concatenatedName)  # break up names
bocp_parsed_names_df <- merge(bocp_names, bocp_parsed_names, by.x = 'concatenatedName', by.y = 'verbatim')

## Identify special cases in the data 
copied_names <- bocp_parsed_names_df %>% # Find the cases of non-unique names 
  group_by(concatenatedName) %>%
  mutate(n = n()) %>%
  filter(n > 1)

hybrids <- bocp_parsed_names_df %>% # Grab a list of names that are hybrids
  filter(grepl("×", concatenatedName)) 

## Remove special cases (all hybrids for now)
bocp_parsed_names_df2 <- bocp_parsed_names_df %>%
  filter(!concatenatedName %in% hybrids$concatenatedName)

### Align with WCVP Taxonomy 
## Using the wcvp backbone, merge on the canonicalFull name 
matched_names <- merge(bocp_parsed_names_df2, wcvp_backbone, by.x = "canonicalfull", by.y = "taxon_name", all.x = TRUE)
matched_names <- matched_names %>% # Simplify to select fields of interest
  select(nameMatch = canonicalfull, concatenatedName, TropicosNameID, authorship, wcvp_accepted_plant_name_id = accepted_plant_name_id, taxon_rank, wcvp_taxon_authors = taxon_authors, family, genus, specificEpithet = species, taxon_status, reviewed, powo_id) 

## Identify Multiple Mapping Cases based on underlying Authorship 
wcvp_backbone_m <- wcvp_backbone %>%
  group_by(taxon_name) %>%
  mutate(num_different_paths = n_distinct(accepted_plant_name_id)) # adjusting this to match wcvp 

multiple_mappings <- wcvp_backbone_m %>%
  filter(num_different_paths > 1) # it would seem that most of these cases are actually subspecific. 

matched_names$multipleMappingsPossible <- ifelse(matched_names$nameMatch %in% multiple_mappings$taxon_name, TRUE, FALSE)

matched_name_summary <- matched_names %>%
  group_by(taxon_status, multipleMappingsPossible) %>%
  count() %>%
  rename(count = n)

# Create some simple summary reports 
num_of_instances_w_multiple_maps_possible <- matched_names %>% 
  filter(multipleMappingsPossible == TRUE) %>% 
  nrow()
num_distinct_names_w_multiple_maps_possible <- matched_names %>% 
  filter(multipleMappingsPossible == TRUE) %>% 
  distinct(nameMatch) %>%  
  nrow()


ggplot(matched_name_summary, aes(x = taxon_status, y = count, fill = multipleMappingsPossible)) +
  geom_col(width = 0.5) +
  geom_text(aes(label = count), vjust = -0.5, position = position_stack(vjust = 1.0)) +
  scale_fill_manual(values = c("TRUE" = "darkred", "FALSE" = "darkgrey"), name = "Multiple Name\nStatus Based On\n Authority") +  # Set colors for TRUE and FALSE
  theme_bw() +
  ggtitle("Taxon Status of North American Plant Name List According to the WorldFloraOnline Backbone") +
  theme(legend.title = element_text(hjust = 0.5),
        legend.text.align = 0.5)
### Appears from my spot checking that taxon_status with Illegitimate, Misapplied, and Invalid all still have accepted_plant_name_id mappings, therefore for now we will treat these classes as synonyms 

## Identify 'recoverable' names based on authorship present in our dataset 
### rm spaces, then preform an exact match 
# Remove spacing in all strings in 'column1' and 'column2'
matched_names <- matched_names[, c("spacelessOurAuthorship", "spacelessWCVPAuthorship") := lapply(.SD, function(x) gsub(" ", "", x)), .SDcols = c("authorship", "wcvp_taxon_authors")]
matched_names$authorshipMatch <- matched_names$spacelessOurAuthorship == matched_names$spacelessWCVPAuthorship

matched_names2 <- matched_names # create a copy of the variable for testing purposes 
matched_names3 <- matched_names2 %>% # Create a new field that designates whether authorships need and sufficient for resolution
  mutate(multMapAuthorshipMatch = case_when(
    multipleMappingsPossible & authorshipMatch == TRUE ~ TRUE, 
    multipleMappingsPossible & authorshipMatch != TRUE ~ FALSE)
  )

resolved_df <- matched_names3 %>% # Grab resolvable names that are multiple mapping but have a usable authorship
  filter(multipleMappingsPossible == TRUE) %>% 
  filter(multMapAuthorshipMatch == TRUE) %>% 
  group_by(nameMatch) %>% 
  mutate(numDiffTaxStatsWithSameAuthor = n_distinct(taxon_status)) %>% # num of distinct paths
  ungroup() %>% 
  filter(numDiffTaxStatsWithSameAuthor == 1) # Filter to the ones that are usable 

matched_names4 <- matched_names3 %>% 
  mutate(multMapResolutionPossible = case_when(
    nameMatch %in% resolved_df$nameMatch & multMapAuthorshipMatch == TRUE ~ TRUE,
    multipleMappingsPossible == FALSE ~ NA, # if we dont need this leave as NA
    TRUE ~ FALSE # else make it false...
  ))

# Report number of names recovered using Authorship as mapping 
mult_map_recovery_name_num <- matched_names4 %>% 
  filter(multMapResolutionPossible) %>% 
  distinct(nameMatch) %>% 
  nrow()

# Remove those redundant names
matched_names5 <- matched_names4 %>% 
  ungroup() %>% 
  filter(multMapResolutionPossible == TRUE & multipleMappingsPossible == TRUE | multipleMappingsPossible == FALSE)

# Check to see if there are any left over cases of multiple mappings being possible prior to filtering. 
check <- matched_names5 %>% 
  group_by(nameMatch) %>% 
  mutate(n = n()) %>% 
  arrange() %>% 
  filter(multipleMappingsPossible == TRUE) # Check shows that there are only 2 names that are redundant that are multiple mappings, however this is just due to the catalogues they derive their info from so safe to filter down. 

# Remove duplicates...
grab_dups <- matched_names5 %>% 
  group_by(nameMatch) %>% 
  mutate(n = n()) %>% 
  filter(multipleMappingsPossible == TRUE & n > 1)

deduped <- grab_dups %>% 
  distinct(nameMatch, .keep_all = TRUE) %>% # verify manually, looks like these cases are fine
  select(!n) # remove extra info 

# Replace dups 
matched_names6 <- matched_names5 %>% 
  filter(!nameMatch %in% grab_dups$nameMatch) %>% 
  rbind(deduped)

# Produce a summary noting alignment status
matched_names6 %>% 
  group_by(taxon_status) %>% 
  summarize(n = n())

# By definition Unplaced taxa cannot be resolved. 

# Split datasets based on whether we can get a resolution via WFO 
matched_names_wcvp <- matched_names6 %>% 
  mutate(source = "wcvp") %>% 
  filter(!is.na(taxon_status) & !taxon_status == "Unplaced")

unmatched_names <- matched_names6 %>% 
  filter(is.na(taxon_status))

# For the unmatched names, use a selection of taxonomic backbones from the taxadb package to match names
match_df <- unmatched_names %>%  # Note that for multiple mapping names, defaullt action is to remove in taxadb 
  mutate(col_id = get_ids(nameMatch, "col")) %>% 
  mutate(col_harmonizedName = get_names(col_id, "col")) %>% 
  mutate(itis_id = get_ids(nameMatch, "itis")) %>% 
  mutate(itis_harmonizedName = get_names(itis_id, "itis")) %>% 
  mutate(gbif_id = get_ids(nameMatch, "gbif")) %>% 
  mutate(gbif_harmonizedName = get_names(gbif_id, "gbif")) %>% 
  select(nameMatch, concatenatedName, TropicosNameID, authorship,col_id, col_harmonizedName, itis_id, itis_harmonizedName, gbif_id, gbif_harmonizedName, spacelessOurAuthorship) %>% 
  mutate(family = NA, taxonRank = NA , scientificNameAuthorship = NA, family = NA, genus = NA, specificEpithet = NA, wfoTaxonomicStatus = NA, references = NA, multipleMappingsPossible = NA,spacelessWFOAuthorship = NA, authorshipMatch = NA, multMapAuthorshipMatch = NA, multMapResolutionPossible = NA)
# Coalesce the names via priorty col > itis > gbif
coal_df <- match_df %>% 
  mutate(nameAligned = coalesce(col_harmonizedName, itis_harmonizedName, gbif_harmonizedName),
         source = case_when(
           !is.na(col_harmonizedName) ~ "col",
           !is.na(itis_harmonizedName) ~ "itis",
           !is.na(gbif_harmonizedName) ~ "gbif",
           TRUE ~ NA_character_
         ),
         taxonomicSourceID = case_when(
           source == "col" ~ col_id, 
           source == "itis" ~ itis_id, 
           source == "gbif" ~ gbif_id,
           TRUE ~ NA_character_
         ))
# Note how many names cant find alignment at all 
coal_name_fails <- coal_df %>% 
  filter(is.na(source))
coal_name_distinct_fails <- coal_name_fails %>% 
  distinct(nameMatch) %>%
  nrow()
# Remove records that do not find a match at all
coal_df <- coal_df %>% 
  filter(!is.na(source))
# There are dups of names, assure that there are no multiple mappings present
col_align_n <- coal_df %>% 
  filter(source == "col") %>% 
  group_by(nameMatch) %>% 
  mutate(n_maps = n_distinct(nameAligned)) 

unique(col_align_n$n_maps) == 1 # true, so we can narrow things down to only distinct 

col_df <- coal_df %>% 
  filter(source == "col") %>% 
  distinct(nameMatch, nameAligned, .keep_all = TRUE) # note that using both fields *should* be redundant 
# ITIS check 
itis_align_n <- coal_df %>% 
  filter(source == "itis") %>% 
  group_by(nameMatch) %>% 
  mutate(n_maps = n_distinct(nameAligned)) 

unique(itis_align_n$n_maps) == 1 # true, so we can narrow things down to only distinct 

itis_df <- coal_df %>% 
  filter(source == "itis") %>% 
  distinct(nameMatch, nameAligned, .keep_all = TRUE) # note that using both fields *should* be redundant 

# GBIF check 
gbif_align_n <- coal_df %>% 
  filter(source == "gbif") %>% 
  group_by(nameMatch) %>% 
  mutate(n_maps = n_distinct(nameAligned)) 

unique(gbif_align_n$n_maps) == 1 # true, so we can narrow things down to only distinct 

gbif_df <- coal_df %>% 
  filter(source == "gbif") %>% 
  distinct(nameMatch, nameAligned, .keep_all = TRUE) # note that using both fields *should* be redundant 



### Attach Synonyms: wcvp
# First find the accepted subspecific assignments
wcvp_backbone_subspecific <- wcvp_backbone %>%
  filter(taxon_rank %in% c("Variety", "Subspecies", "Form", "Subvariety", "Subform", "Convariety", "subspecioid")) %>%
  filter(taxon_status == "Accepted")
# Find what their parent taxon is 
wcvp_backbone_subspecific <- wcvp_backbone_subspecific %>%
  mutate(parentTaxon = stringr::word(taxon_name, 1, 2)) %>%
  mutate(taxon_status = "Synonym") # manually change these to Synonyms.

updated_wcvp_backbone <- merge(wcvp_backbone, wcvp_backbone_subspecific, by.x = "taxon_name", by.y = "parentTaxon")
updated_wcvp_backbone <- updated_wcvp_backbone %>%
  mutate(statusMatch = ifelse(taxon_status.x == taxon_status.y, TRUE, FALSE)) %>%
  filter(statusMatch == FALSE) %>%  # grab the cases we designed
  #mutate(accepted_plant_name_id.y = ifelse(accepted_plant_name_id.y == "", accepted_plant_name_id.x, accepted_plant_name_id.y)) %>%
  rename(accepted_plant_name_id = accepted_plant_name_id.x, parentTaxon = taxon_name, taxon_name = taxon_name.y,
         taxon_authors = taxon_authors.y, taxon_status = taxon_status.y, powo_id = powo_id.y) %>%
  select(accepted_plant_name_id, taxon_name, taxon_authors, parentTaxon, taxon_status, powo_id)

test_mult_maps <- updated_wcvp_backbone %>%
  group_by(taxon_name) %>%
  mutate(num_different_paths = n_distinct(taxon_status))

test_multiple_mappings <- test_mult_maps  %>% # tis all good!
  filter(num_different_paths > 1)

wcvp_backbone2 <- wcvp_backbone %>% 
   filter(!taxon_name %in% updated_wcvp_backbone$taxon_name)
 
wcvp_backbone2 <- wcvp_backbone2 %>%
  mutate(parentTaxon = NA) %>%
  select(accepted_plant_name_id, taxon_name, taxon_authors, parentTaxon, taxon_status, powo_id) %>%
  rbind(updated_wcvp_backbone)

# WCVP Mapping of Names 
wcvp_accepted_mapping <- wcvp_backbone2[taxon_status == "Accepted"]
wcvp_synonym_mapping <- wcvp_backbone2[taxon_status %in% c("Synonym", "Illegitimate", "Invalid", "Misapplied", "Orthographic")] # all will be regarded as equivalent for the purposes of mapping.
non_unique_accepted_names <- wcvp_accepted_mapping[, .N, by = "taxon_name"]# Identify non unique names 
# Create a simplified version of each respective dt, then attach synonyms
accepted_mapping <- wcvp_accepted_mapping[, c("accepted_plant_name_id", "taxon_name", "powo_id"), with = FALSE]
synonym_mapping <- wcvp_synonym_mapping[, c("accepted_plant_name_id", "taxon_name", "powo_id", "taxon_authors"), with = FALSE]

# Append acceptedNameUsageID with that of its parent taxon. 
matched_names_wcvp <- matched_names_wcvp %>% data.table()  %>% as.data.table() # awful black magic is required to solve internal indexing errors
wcvp_df_accepted <- matched_names_wcvp[taxon_status == "Accepted"]
wcvp_df_synonym <- matched_names_wcvp[taxon_status %in% c("Synonym", "Illegitimate", "Invalid", "Misapplied", "Orthographic")]
wcvp_relations <- merge(accepted_mapping, synonym_mapping, by.x = "accepted_plant_name_id", by.y = "accepted_plant_name_id", all.x = TRUE)
wcvp_relations <- wcvp_relations %>% 
  rename(acceptedName = taxon_name.x, synonyms = taxon_name.y,
         acceptedNamePOWO_ID = powo_id.x, synonymPOWO_ID = powo_id.y) 
# subset backbone for our list of name matched taxa
wcvp_relations_f <- wcvp_relations %>%  
  filter(acceptedName %in% wcvp_df_accepted$nameMatch | synonyms %in% wcvp_df_synonym$nameMatch) %>% 
  rename(synonym_taxon_authors = taxon_authors) %>% 
  mutate(catalogID = "wcvp") %>% 
  distinct(accepted_plant_name_id, acceptedName, synonyms, catalogID, acceptedNamePOWO_ID, synonymPOWO_ID, .keep_all = TRUE)

## other catalogues
# Load catalogues 
col <-taxa_tbl("col") %>%
  select(scientificName,taxonRank,acceptedNameUsageID,taxonomicStatus) %>%
  filter(taxonRank == "species")
col %>% show_query()
col_names <- col %>% collect() #retrieve results
col_names <- as.data.table(col_names)
col_accepted_names <- col_names[taxonomicStatus == "accepted"] 
col_synonymous_names <- col_names[taxonomicStatus == "synonym"]
col_synonymous_names_d <- col_synonymous_names %>% 
  distinct(scientificName, taxonRank, taxonomicStatus, .keep_all = TRUE) # there are cases of duplicates. 
# Check for ambiguous mapping on the acceptedNameUsageID as that would be problematic
col_non_unique_accepted_names_check <- col_accepted_names %>% 
  group_by(scientificName) %>% 
  mutate(multipleAcceptableUsageIDs = n_distinct(acceptedNameUsageID))
col_non_unique_synonym_names_check <- col_synonymous_names %>% 
  group_by(scientificName) %>% 
  mutate(multipleAcceptableUsageIDs = n_distinct(acceptedNameUsageID))
# remove ambigious mapping as there is no simple way to deal with it 
col_non_unique_accepted_names <- col_non_unique_accepted_names_check %>%
  filter(multipleAcceptableUsageIDs >= 2) # vast majority of these are just genera, but there are some cases of genus+specificEpithet having 2 valid mapping paths. This is just too much to deal with though
col_non_unique_synonym_names <- col_non_unique_synonym_names_check %>% 
  filter(multipleAcceptableUsageIDs >= 2)
# It appears that many of these non-unique mappings are simply due to non-existent IDs cluttering the mapping (presumably these a relictual and uncleaned from the db), lets remove them if they dont exist as an accepted name and remap to try again. 
col_accepted_usage_ids <- unique(col_accepted_names$acceptedNameUsageID) # collect the name usage IDs
col_synonymous_names2 <- col_synonymous_names[acceptedNameUsageID %in% col_accepted_usage_ids]
col_non_unique_synonym_names_check2 <- col_synonymous_names2 %>% 
  group_by(scientificName) %>% 
  mutate(multipleAcceptableUsageIDs = n_distinct(acceptedNameUsageID))
col_non_unique_synonym_names2 <- col_non_unique_synonym_names_check2 %>% 
  filter(multipleAcceptableUsageIDs >= 2)
# Those are true issues (likely related to authorship), we'll have to remove them from our alignment. 
col_synonymous_names_drop <- col_synonymous_names[!scientificName %in% col_non_unique_synonym_names2$scientificName]
col_accepted_names_drop <- col_accepted_names[!scientificName %in% col_non_unique_accepted_names]
# Now, filter out relations according to our data...
col_df <- col_df[, colTaxonomicStatus := ifelse(nameMatch == nameAligned, "ACCEPTED", "SYNONYM")]
col_df_accepted <- col_df[colTaxonomicStatus == "ACCEPTED"]
col_df_synonym <- col_df[colTaxonomicStatus == "SYNONYM"]
# Now create the relations 
col_relations <- merge(col_accepted_names_drop, col_synonymous_names_drop, by = "acceptedNameUsageID", all.x = TRUE)
col_relations <- col_relations %>% 
  rename(acceptedName = scientificName.x, synonyms = scientificName.y) %>% 
  mutate(catalogID = "col")
# Now filter the relations to only those that matter for the names we're resolving in our dataset
col_relations_f <- col_relations %>%  
  filter(acceptedName %in% col_df_accepted$nameAligned | synonyms %in% col_df_synonym$nameMatch) %>% 
  rename(taxonID = acceptedNameUsageID) %>% 
  select(taxonID, acceptedName, synonyms, catalogID)

### ITIS 
itis <-taxa_tbl("itis") %>%
  select(scientificName,taxonRank,acceptedNameUsageID,taxonomicStatus) %>%
  filter(taxonRank == "species")
itis %>% show_query()
itis_names <- itis %>% collect() #retrieve results
itis_names <- as.data.table(itis_names)
itis_accepted_names <- itis_names[taxonomicStatus == "accepted"] 
itis_synonymous_names <- itis_names[taxonomicStatus == "synonym"]
itis_synonymous_names_d <- itis_synonymous_names %>% 
  distinct(scientificName, taxonRank, taxonomicStatus, .keep_all = TRUE) # there are cases of duplicates. 
# Check for ambiguous mapping on the acceptedNameUsageID as that would be problematic
itis_non_unique_accepted_names_check <- itis_accepted_names %>% 
  group_by(scientificName) %>% 
  mutate(multipleAcceptableUsageIDs = n_distinct(acceptedNameUsageID))
itis_non_unique_synonym_names_check <- itis_synonymous_names %>% 
  group_by(scientificName) %>% 
  mutate(multipleAcceptableUsageIDs = n_distinct(acceptedNameUsageID))
# remove ambigious mapping as there is no simple way to deal with it 
itis_non_unique_accepted_names <- itis_non_unique_accepted_names_check %>%
  filter(multipleAcceptableUsageIDs >= 2) # vast majority of these are just genera, but there are some cases of genus+specificEpithet having 2 valid mapping paths. This is just too much to deal with though
itis_non_unique_synonym_names <- itis_non_unique_synonym_names_check %>% 
  filter(multipleAcceptableUsageIDs >= 2)
# It appears that many of these non-unique mappings are simply due to non-existent IDs cluttering the mapping (presumably these a relictual and uncleaned from the db), lets remove them if they dont exist as an accepted name and remap to try again. 
itis_accepted_usage_ids <- unique(itis_accepted_names$acceptedNameUsageID) # itislect the name usage IDs
itis_synonymous_names2 <- itis_synonymous_names[acceptedNameUsageID %in% itis_accepted_usage_ids]
itis_non_unique_synonym_names_check2 <- itis_synonymous_names2 %>% 
  group_by(scientificName) %>% 
  mutate(multipleAcceptableUsageIDs = n_distinct(acceptedNameUsageID))
itis_non_unique_synonym_names2 <- itis_non_unique_synonym_names_check2 %>% 
  filter(multipleAcceptableUsageIDs >= 2)
# Those are true issues (likely related to authorship), we'll have to remove them from our alignment. 
itis_synonymous_names_drop <- itis_synonymous_names[!scientificName %in% itis_non_unique_synonym_names2$scientificName]
itis_accepted_names_drop <- itis_accepted_names[!scientificName %in% itis_non_unique_accepted_names]
# Now, filter out relations according to our data...
itis_df <- itis_df[, itisTaxonomicStatus := ifelse(nameMatch == nameAligned, "ACCEPTED", "SYNONYM")]
itis_df_accepted <- itis_df[itisTaxonomicStatus == "ACCEPTED"]
itis_df_synonym <- itis_df[itisTaxonomicStatus == "SYNONYM"]
# Now create the relations 
itis_relations <- merge(itis_accepted_names_drop, itis_synonymous_names_drop, by = "acceptedNameUsageID", all.x = TRUE)
itis_relations <- itis_relations %>% 
  rename(acceptedName = scientificName.x, synonyms = scientificName.y) %>% 
  mutate(catalogID = "itis")
# Now filter the relations to only those that matter for the names we're resolving in our dataset
itis_relations_f <- itis_relations %>%  
  filter(acceptedName %in% itis_df_accepted$nameAligned | synonyms %in% itis_df_synonym$nameMatch) %>% 
  rename(taxonID = acceptedNameUsageID) %>% 
  select(taxonID, acceptedName, synonyms, catalogID)

### GBIF 
gbif <-taxa_tbl("gbif") %>%
  select(scientificName,taxonRank,acceptedNameUsageID,taxonomicStatus) %>%
  filter(taxonRank == "species")
gbif %>% show_query()
gbif_names <- gbif %>% collect() #retrieve results
gbif_names <- as.data.table(gbif_names)
gbif_accepted_names <- gbif_names[taxonomicStatus == "accepted"] 
gbif_synonymous_names <- gbif_names[taxonomicStatus == "synonym"]
gbif_synonymous_names_d <- gbif_synonymous_names %>% 
  distinct(scientificName, taxonRank, taxonomicStatus, .keep_all = TRUE) # there are cases of duplicates. 
# Check for ambiguous mapping on the acceptedNameUsageID as that would be problematic
gbif_non_unique_accepted_names_check <- gbif_accepted_names %>% 
  group_by(scientificName) %>% 
  mutate(multipleAcceptableUsageIDs = n_distinct(acceptedNameUsageID))
gbif_non_unique_synonym_names_check <- gbif_synonymous_names %>% 
  group_by(scientificName) %>% 
  mutate(multipleAcceptableUsageIDs = n_distinct(acceptedNameUsageID))
# remove ambigious mapping as there is no simple way to deal with it 
gbif_non_unique_accepted_names <- gbif_non_unique_accepted_names_check %>%
  filter(multipleAcceptableUsageIDs >= 2) # vast majority of these are just genera, but there are some cases of genus+specificEpithet having 2 valid mapping paths. This is just too much to deal with though
gbif_non_unique_synonym_names <- gbif_non_unique_synonym_names_check %>% 
  filter(multipleAcceptableUsageIDs >= 2)
# It appears that many of these non-unique mappings are simply due to non-existent IDs cluttering the mapping (presumably these a relictual and uncleaned from the db), lets remove them if they dont exist as an accepted name and remap to try again. 
gbif_accepted_usage_ids <- unique(gbif_accepted_names$acceptedNameUsageID) # gbiflect the name usage IDs
gbif_synonymous_names2 <- gbif_synonymous_names[acceptedNameUsageID %in% gbif_accepted_usage_ids]
gbif_non_unique_synonym_names_check2 <- gbif_synonymous_names2 %>% 
  group_by(scientificName) %>% 
  mutate(multipleAcceptableUsageIDs = n_distinct(acceptedNameUsageID))
gbif_non_unique_synonym_names2 <- gbif_non_unique_synonym_names_check2 %>% 
  filter(multipleAcceptableUsageIDs >= 2)
# Those are true issues (likely related to authorship), we'll have to remove them from our alignment. 
gbif_synonymous_names_drop <- gbif_synonymous_names[!scientificName %in% gbif_non_unique_synonym_names2$scientificName]
gbif_accepted_names_drop <- gbif_accepted_names[!scientificName %in% gbif_non_unique_accepted_names]
# Now, filter out relations according to our data...
gbif_df <- gbif_df[, gbifTaxonomicStatus := ifelse(nameMatch == nameAligned, "ACCEPTED", "SYNONYM")]
gbif_df_accepted <- gbif_df[gbifTaxonomicStatus == "ACCEPTED"]
gbif_df_synonym <- gbif_df[gbifTaxonomicStatus == "SYNONYM"]
# Now create the relations 
gbif_relations <- merge(gbif_accepted_names_drop, gbif_synonymous_names_drop, by = "acceptedNameUsageID", all.x = TRUE)
gbif_relations <- gbif_relations %>% 
  rename(acceptedName = scientificName.x, synonyms = scientificName.y) %>% 
  mutate(catalogID = "gbif")
# Now filter the relations to only those that matter for the names we're resolving in our dataset
gbif_relations_f <- gbif_relations %>%  
  filter(acceptedName %in% gbif_df_accepted$nameAligned | synonyms %in% gbif_df_synonym$nameMatch) %>% 
  rename(taxonID = acceptedNameUsageID) %>% 
  select(taxonID, acceptedName, synonyms, catalogID)

### Now compare accepted/synonymous names across catalogues to ensure exclusive status
# Make a reproducible sample to sample testcase this on...
x <- data.frame(acceptedName = c("X", "X", "Y", "Y", "Y", "Z", "V", "T"),
                synonymName = c("xxx", "xyx", "yyy", "yyy", "yml", "zzz", "vvv", "ttt"), 
                backbone = c("col", "col", "col", "col", "col", "col", "col", "col"))
y <- data.frame(acceptedName = c("A", "A", "B", "C", "D"),
                synonymName = c("aaa", "axa", "bbb", "ccc", "ddd"), 
                backbone = c("itis", "itis", "itis", "itis", "itis"))
z <- data.frame(acceptedName = c("J", "J", "Xi", "T", "A"),
                synonymName = c("jjj", "xyx", "xxx", "tit", "uuu"), 
                backbone = c("gbif", "gbif", "gbif", "gbif", "gbif")) 

xy <- bind_rows(x,y) # an example where no duplication should take place

xyz <- bind_rows(x,y,z) # an example where duplication definitely takes place

add_duplication_flags <- function(df) {
  df <- df %>%
    group_by(acceptedName) %>% 
    mutate(a_dupes = n_distinct(backbone)) %>% 
    ungroup() %>% 
    group_by(synonymName) %>% 
    mutate(s_breaks = n_distinct(backbone)) %>% 
    ungroup()
  return(df)
}

# Apply the function to your combined data frames
xy_with_flags <- add_duplication_flags(xy)
xyz_with_flags <- add_duplication_flags(xyz)

# Try something else out
xyz_with_flags <- xyz %>% 
  group_by(acceptedName) %>% 
  mutate(a_dupes = n_distinct(backbone)) %>% 
  ungroup() %>% 
  group_by(synonymName) %>% 
  mutate(s_breaks = n_distinct(backbone)) %>% 
  ungroup()


# Combine the relational datatables
relational_tb <- bind_rows(wcvp_relations_f, col_relations_f, itis_relations_f, gbif_relations_f )

relational_tb <- relational_tb %>% 
  mutate(synonyms = case_when(
    is.na(synonyms) == TRUE ~ acceptedName, 
    is.na(synonyms) == FALSE ~ synonyms
  )) # this is mostly a fix, though we dont have unique cases now...keep in mind 
add_duplication_flags <- function(df) {
  df <- df %>%
    group_by(acceptedName) %>%
    mutate(a_dupes = n_distinct(catalogID)) %>%
    ungroup() %>%
    group_by(synonyms) %>%
    mutate(s_breaks = n_distinct(catalogID)) %>%
    ungroup()
  return(df)
}

relational_tb_wflags <- add_duplication_flags(relational_tb)
priority_order <- c("wcvp", "col", "itis", "gbif")

relational_tb_filtered <- relational_tb_wflags %>% 
  mutate(priority = case_when(
    a_dupes > 1 | s_breaks > 1  ~ match(catalogID, priority_order), 
    TRUE ~ 0)) %>% 
  group_by(acceptedName) %>% 
  filter(priority == min(priority))  %>% 
  ungroup() %>% 
  group_by(synonyms) %>% 
  filter(priority == min(priority)) %>% 
  ungroup()

# Final Resolution Info:
final_df <- relational_tb_filtered

final_accepted_count <- length(unique(final_df$acceptedName))
final_synonym_count <- final_df %>% 
  mutate(copy = case_when(
    acceptedName == synonyms ~ TRUE,
    acceptedName != synonyms ~ FALSE
  )) %>% 
  filter(copy == FALSE) %>% 
  group_by(synonyms) %>% 
  mutate(n = n())

# Create a join to list authorship and multiple mapping statuses 
accepted_matched_names_wcvp <- filter(matched_names_wcvp, taxon_status == "Accepted")
accepted_half_merge <- merge(final_df, accepted_matched_names_wcvp, by.x = "acceptedNamePOWO_ID", by.y = "powo_id")
synonym_matched_names_wcvp <- filter(matched_names_wcvp, taxon_status != "Accepted")
synonym_half_merge <- merge(final_df, synonym_matched_names_wcvp, by.x = "synonymPOWO_ID", by.y = "powo_id") 
final_df2 <- rbind(accepted_half_merge, synonym_half_merge)
final_df3 <- distinct(final_df2, accepted_plant_name_id, acceptedName, acceptedNamePOWO_ID, synonyms, synonymPOWO_ID, synonym_taxon_authors, catalogID, authorship, multipleMappingsPossible, spacelessOurAuthorship, spacelessWCVPAuthorship, authorshipMatch, multMapAuthorshipMatch, multMapResolutionPossible)
# There are cases where there are repeated synonyms for the same accepted name, this is presumably due to synonyms having multiple mapping paths depending on their underlying authorship. Lets resolve this. 
#final_df2 <- merge(final_df, matched_names_wcvp, by.x = "taxonID", by.y = "taxonomicSourceID", all.x = TRUE)
#final_df2 <- select(final_df2, names(final_df), scientificNameAuthorship, multipleMappingsPossible, spacelessWFOAuthorship, authorshipMatch, multMapResolutionPossible)
final_df4 <- final_df3 %>% # rename for clarity that our multiple mapping resolution was depedent on the acceptedNames, now we need to deal with the synonyms
  rename(acceptedNameScientificNameAuthorship = authorship,
         acceptedNameMultMapsPossible = multipleMappingsPossible,
         acceptedNameSpacelessWCVPAuthorship = spacelessWCVPAuthorship, 
         acceptedNameAuthorshipMatch = authorshipMatch,
         acceptedNameMultMapResolutionPossible =  multMapResolutionPossible, 
         synonymNameWCVPAuthorship = synonym_taxon_authors)

# Assure multiple maps prior to simplifying 
multiple_s_mappings <- final_df4 %>% 
  group_by(synonyms) %>% 
  mutate(numDiffAcceptedNames = length(unique((acceptedName))))
syn_mult_maps <- filter(final_df4 , multiple_s_mappings$numDiffAcceptedNames > 1)

final_df4$synonymMultMapPossible <- ifelse(final_df4$synonyms %in% syn_mult_maps$synonyms, TRUE, FALSE)

final_synonym_count <- final_df4 %>% 
  mutate(copy = case_when(
    acceptedName == synonyms ~ TRUE,
    acceptedName != synonyms ~ FALSE
  )) %>% 
  filter(copy == FALSE) %>% 
  distinct(synonyms, synonymNameWCVPAuthorship) %>% 
  nrow()

final_multiple_map_on_synonym_count <- final_df4 %>% 
  mutate(copy = case_when(
    acceptedName == synonyms ~ TRUE,
    acceptedName != synonyms ~ FALSE
  )) %>% 
  filter(copy == FALSE & synonymMultMapPossible == TRUE) %>% 
  distinct(synonyms, synonymNameWCVPAuthorship) %>% 
  nrow()



final_df5 <- final_df4 %>% 
  mutate(synonymNameSpacelessWCVPAuthorship = gsub(" ", "", synonymNameWCVPAuthorship))

final_df6 <- final_df5 %>% 
  mutate(remove = case_when(
    catalogID == "wcvp" & is.na(acceptedNameMultMapsPossible) ~ TRUE, 
    TRUE ~ FALSE )
) %>% 
  filter(!remove == TRUE) %>% 
  select(-remove)


fwrite(final_df6, "/home/jt-miller/Gurlab/BoCP/data/processed/finalized_name_alignment_wcvp.csv")


















































