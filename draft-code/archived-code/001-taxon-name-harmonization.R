# 001-taxon-name-harmonization (matches prototype 001-taxon-name-alignment-wcvp-rw.Rs)
# Author: JT Miller 
# Date: 06-25-2024
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

## Remove special cases (all hybrids)
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

## Identify 'recoverable' names based on authorship present in our dataset 
### rm spaces & capitilization, then perform an exact match 
matched_names <- matched_names[, c("spacelessOurAuthorship", "spacelessWCVPAuthorship") := 
                                 lapply(.SD, function(x) tolower(gsub(" ", "", x))), 
                               .SDcols = c("authorship", "wcvp_taxon_authors")]
matched_names$authorshipMatch <- matched_names$spacelessOurAuthorship == matched_names$spacelessWCVPAuthorship

## Create new a new field called multMapAuthorship match that requires authorship to match when multipleMappings is Possible
matched_names <- matched_names %>% 
  mutate(multMapAuthorshipMatch = case_when(
    multipleMappingsPossible & authorshipMatch == TRUE ~ TRUE, 
    multipleMappingsPossible & authorshipMatch != TRUE ~ FALSE))

## Check to make sure that mapping is 1:1 while accounting for authorship in multiple mappings.
resolvable_names_df <- matched_names %>% 
  filter(multipleMappingsPossible == TRUE) %>% 
  filter(multMapAuthorshipMatch == TRUE) %>% 
  group_by(nameMatch) %>% 
  mutate(numDiffTaxStatsWithSameAuthor = n_distinct(taxon_status)) %>% # num of distinct paths
  ungroup() %>% 
  filter(numDiffTaxStatsWithSameAuthor == 1) # Filter to the ones that are usable 

## Create a field called multMapResolutionPossible that denotes whether in case there is multiple maps if AuthorshipMatch is possible, else other cases
matched_names <- matched_names %>% 
  mutate(multMapResolutionPossible = case_when(
    nameMatch %in% resolvable_names_df$nameMatch & multMapAuthorshipMatch == TRUE ~ TRUE,
    multipleMappingsPossible == FALSE ~ NA, # if we dont need this leave as NA
    TRUE ~ FALSE # else make it false...
  ))

## Report number of names recovered using Authorship as mapping 
matched_names %>% 
  filter(multMapResolutionPossible) %>% 
  distinct(nameMatch) %>% 
  nrow()

## Remove the multiple matching names that cannot be resolved via authorship
matched_names <- matched_names %>% 
  filter(multMapResolutionPossible == TRUE & multipleMappingsPossible == TRUE | multipleMappingsPossible == FALSE) %>% 
  mutate(authorshipInAcceptedName = ifelse(taxon_status == "Accepted", TRUE, FALSE))

## Create a summary of alignment status 
alignment_summary <- matched_names %>% 
  group_by(taxon_status) %>% 
  summarize(n = n()) 
print(alignment_summary)

### Split the dataset. Taxa that have any form of alignment outside of Unplaced and NA are ok for wcvp.
## Unplaced names will be thrown away as the creators of wcvp have issues with these names 
## NA names will be fed through the catalogs COL, ITIS, and GBIF in order to assign aligmnet 
matched_names_wcvp <- matched_names %>% 
  mutate(source = "wcvp") %>% 
  filter(!taxon_status == "Unplaced" & !is.na(taxon_status))

unmatched_names <- matched_names %>% 
  filter(is.na(taxon_status))

## For unmatched names, harmonize against COL, ITIS, and GBIF backbones
match_df <- unmatched_names %>%  
  mutate(col_id = get_ids(nameMatch, "col")) %>% 
  mutate(col_harmonizedName = get_names(col_id, "col")) %>% 
  mutate(itis_id = get_ids(nameMatch, "itis")) %>% 
  mutate(itis_harmonizedName = get_names(itis_id, "itis")) %>% 
  mutate(gbif_id = get_ids(nameMatch, "gbif")) %>% 
  mutate(gbif_harmonizedName = get_names(gbif_id, "gbif")) %>% 
  select(nameMatch, concatenatedName, TropicosNameID, authorship,col_id, 
         col_harmonizedName, itis_id, itis_harmonizedName, gbif_id, 
         gbif_harmonizedName, spacelessOurAuthorship) %>% 
  mutate(family = NA, taxonRank = NA , scientificNameAuthorship = NA, 
         family = NA, genus = NA, specificEpithet = NA, wfoTaxonomicStatus = NA, 
         references = NA, multipleMappingsPossible = NA,spacelessWFOAuthorship = NA, 
         authorshipMatch = NA, multMapAuthorshipMatch = NA, multMapResolutionPossible = NA)

### Coalesce the names via priorty col > itis > gbif
coal_df <- match_df %>% 
  mutate(nameAligned = coalesce(col_harmonizedName, itis_harmonizedName, gbif_harmonizedName),
         source = case_when(
           !is.na(col_harmonizedName) ~ "col",
           !is.na(itis_harmonizedName) ~ "itis",
           !is.na(gbif_harmonizedName) ~ "gbif",
           TRUE ~ NA_character_),
         taxonomicSourceID = case_when(
           source == "col" ~ col_id, 
           source == "itis" ~ itis_id, 
           source == "gbif" ~ gbif_id,
           TRUE ~ NA_character_)
  )

## Report successes and fails in coalescent resolution on the unmatched names 
coal_df %>% filter(is.na(source)) %>% distinct(nameMatch) %>% nrow()
coal_df %>% filter(!is.na(source)) %>% distinct(nameMatch) %>% nrow()

## Remove names that failed to find a match in coalescent resolution 
coal_df <- filter(coal_df, !is.na(source))

### Attach Synonyms to each Backbone
## wvcp
# Organize by accepted vs everything else that we're lumping into synonyms
wcvp_accepted_mapping <- wcvp_backbone[taxon_status == "Accepted"]
wcvp_synonym_mapping <- wcvp_backbone[taxon_status %in% c("Synonym", "Illegitimate", "Invalid", "Misapplied", "Orthographic")]

# Select cols of interest
accepted_mapping <- wcvp_accepted_mapping[, c("accepted_plant_name_id", "taxon_name", "powo_id", "taxon_authors"), with = FALSE]
synonym_mapping <- wcvp_synonym_mapping[, c("accepted_plant_name_id", "taxon_name", "powo_id", "taxon_authors"), with = FALSE]

# Create a relational table for wcvp relations
wcvp_relations <- merge(accepted_mapping, synonym_mapping, 
                        by.x = "accepted_plant_name_id",
                        by.y = "accepted_plant_name_id", 
                        all.x = TRUE)

# Rename fields for clarity
wcvp_relations <- wcvp_relations %>% 
  rename(acceptedName = taxon_name.x,
         synonym = taxon_name.y,
         accepted_powo_id = powo_id.x, 
         synonym_powo_id = powo_id.y, 
         accepted_taxon_authors = taxon_authors.x, 
         synonym_taxon_authors = taxon_authors.y)

## Deal with infraspecifics, we are including infraspecifics as part of the underlying taxonomy 
# Detect names with infraspecifics (more than 2 words)
wcvp_relations <- wcvp_relations %>% 
  mutate(word_count = str_count(acceptedName, "\\s+") + 1)

# Duplicate these instances into a synonym for the name (so we know they exist)
# then create a acceptedNameParent field that is just genus + specificEpithet of the acceptedName. 
# Being sure to associate the authorship if present into the synonymy. 
wcvp_relations <- wcvp_relations %>% 
  bind_rows(wcvp_relations %>% 
              filter(word_count > 2) %>% 
              distinct(acceptedName, .keep_all = TRUE) %>% 
              mutate(synonym = acceptedName, 
                     synonym_taxon_authors = accepted_taxon_authors, 
                     synonym_powo_id = accepted_powo_id)) %>% 
  mutate(acceptedNameParent = ifelse(word_count > 2, 
                                     stringr::word(acceptedName, 1, 2), 
                                     acceptedName))

## Filter name relations down to North America Plants, utilizing an OR statement in acceptedNames or Synonyms
wcvp_df_accepted <- matched_names_wcvp[taxon_status == "Accepted"]
wcvp_df_synonym <- matched_names_wcvp[taxon_status != "Accepted"] # all other cases
wcvp_relations_na <- wcvp_relations %>% 
  filter(acceptedName %in% wcvp_df_accepted$nameMatch | synonym %in% wcvp_df_synonym$nameMatch | acceptedNameParent %in% wcvp_df_accepted$nameMatch | acceptedNameParent %in% wcvp_df_synonym$nameMatch) %>% 
  select(-word_count) %>% # no longer necessary. 
  mutate(catalogID = "wcvp")

### COL 
col <- taxa_tbl("col") %>% 
  select(scientificName,taxonRank,acceptedNameUsageID,taxonomicStatus) %>%
  filter(taxonRank == "species")

# Load into memory, convert to data.table for faster operations
col_names <- col %>% collect() # load into memory
col_names <- as.data.table(col_names)

# Organize by taxonomic status
col_accepted_names <- col_names[taxonomicStatus == "accepted"]
col_synonymous_names <- col_names[taxonomicStatus == "synonym"]

# Check for ambigious mapping in COL on accepted names
col_non_uniq_accepted_names <- col_accepted_names %>% 
  group_by(scientificName) %>% 
  mutate(multipleAcceptableUsageIDs = n_distinct(acceptedNameUsageID)) %>% 
  filter(multipleAcceptableUsageIDs >= 2)

# There are non-existent acceptedNameUsageIDs cluttering in the synonyms, 
# first remove these, then check for multiple cases
col_accepted_usage_ids <- unique(col_accepted_names$acceptedNameUsageID)
col_synonymous_names <- col_synonymous_names[acceptedNameUsageID %in% col_accepted_usage_ids]
col_non_uniq_synonymous_names <- col_synonymous_names %>% 
  group_by(scientificName) %>% 
  mutate(multipleAcceptableUsageIDs = n_distinct(acceptedNameUsageID)) %>% 
  filter(multipleAcceptableUsageIDs >= 2)

# these non-unique cases are true issues likely related to authorship. We're going to choose to drop them
col_accepted_names_drop <- col_accepted_names[!scientificName %in% col_non_uniq_accepted_names$scientificName]
col_synonymous_names_drop <- col_synonymous_names[!scientificName %in% col_non_uniq_synonymous_names$scientificName]

# create relational table for COL taxonomy 
col_relations <- merge(col_accepted_names_drop, col_synonymous_names_drop, by = "acceptedNameUsageID", all.x = TRUE)
col_relations <- col_relations %>% 
  rename(acceptedName = scientificName.x,
         synonym = scientificName.y) %>% 
  mutate(catalogID = "col")

## Deal with infraspecifics in COL, also deal use regular expressions to deal with issues regarding [sensu lato]
col_relations <- col_relations %>% 
  mutate(acceptedName = gsub("\\[.*?\\]", "", acceptedName)) %>% 
  mutate(acceptedName = gsub("  ", " ", acceptedName)) %>% 
  mutate(synonym = gsub("\\[.*?\\]", "", synonym)) %>% 
  mutate(synonym = gsub("  ", " ", synonym)) %>%
  mutate(word_count = str_count(acceptedName, "\\s+") + 1) 

# Duplicate these infraspecific instances into the synonymy so we know that they exist
col_relations <- col_relations %>% 
  rename(accepted_name_taxon_rank = taxonRank.x,
         synonym_name_taxon_rank = taxonRank.y,
         accepted_taxon_status = taxonomicStatus.x,
         synonym_taxon_status = taxonomicStatus.y) %>% 
  bind_rows(col_relations %>% 
              filter(word_count > 2) %>% 
              distinct(acceptedName, .keep_all = TRUE) %>% 
              mutate(synonym = acceptedName)) %>% 
  mutate(acceptedNameParent = ifelse(word_count > 2,
                                     stringr::word(acceptedName, 1, 2),
                                     acceptedName)
  )
# Filter name relations in COL down to North American plants 
col_df <- coal_df %>% 
  filter(source == "col")
col_df <- col_df[, colTaxonomicStatus := ifelse(nameMatch == nameAligned, "ACCEPTED", "SYNONYM")]
col_df_accepted <- col_df[colTaxonomicStatus == "ACCEPTED"]
col_df_synonym <- col_df[colTaxonomicStatus == "SYNONYM"]
col_relations_na <- col_relations %>% 
  filter(acceptedName %in% col_df_accepted$nameMatch | synonym %in% col_df_synonym$nameMatch) %>% 
  select(-word_count) %>% 
  rename(taxonID = acceptedNameUsageID) %>% 
  select(taxonID, acceptedName, acceptedNameParent, synonym, catalogID)

### ITIS 
itis <- taxa_tbl("itis") %>% 
  select(scientificName,taxonRank,acceptedNameUsageID,taxonomicStatus) %>%
  filter(taxonRank == "species")

# Load into memory, convert to data.table for faster operations
itis_names <- itis %>% collect() # load into memory
itis_names <- as.data.table(itis_names)

# Organize by taxonomic status
itis_accepted_names <- itis_names[taxonomicStatus == "accepted"]
itis_synonymous_names <- itis_names[taxonomicStatus == "synonym"]

# Check for ambigious mapping in itis on accepted names
itis_non_uniq_accepted_names <- itis_accepted_names %>% 
  group_by(scientificName) %>% 
  mutate(multipleAcceptableUsageIDs = n_distinct(acceptedNameUsageID)) %>% 
  filter(multipleAcceptableUsageIDs >= 2)

# There are non-existent acceptedNameUsageIDs cluttering in the synonyms, 
# first remove these, then check for multiple cases
itis_accepted_usage_ids <- unique(itis_accepted_names$acceptedNameUsageID)
itis_synonymous_names <- itis_synonymous_names[acceptedNameUsageID %in% itis_accepted_usage_ids]
itis_non_uniq_synonymous_names <- itis_synonymous_names %>% 
  group_by(scientificName) %>% 
  mutate(multipleAcceptableUsageIDs = n_distinct(acceptedNameUsageID)) %>% 
  filter(multipleAcceptableUsageIDs >= 2)

# these non-unique cases are true issues likely related to authorship. We're going to choose to drop them
itis_accepted_names_drop <- itis_accepted_names[!scientificName %in% itis_non_uniq_accepted_names$scientificName]
itis_synonymous_names_drop <- itis_synonymous_names[!scientificName %in% itis_non_uniq_synonymous_names$scientificName]

# create relational table for itis taxonomy 
itis_relations <- merge(itis_accepted_names_drop, itis_synonymous_names_drop, by = "acceptedNameUsageID", all.x = TRUE)
itis_relations <- itis_relations %>% 
  rename(acceptedName = scientificName.x,
         synonym = scientificName.y) %>% 
  mutate(catalogID = "itis")

## Deal with infraspecifics in itis
itis_relations <- itis_relations %>% 
  mutate(word_count = str_count(acceptedName, "\\s+") + 1) 

# Duplicate these infraspecific instances into the synonymy so we know that they exist
itis_relations <- itis_relations %>% 
  rename(accepted_name_taxon_rank = taxonRank.x,
         synonym_name_taxon_rank = taxonRank.y,
         accepted_taxon_status = taxonomicStatus.x,
         synonym_taxon_status = taxonomicStatus.y) %>% 
  bind_rows(itis_relations %>% 
              filter(word_count > 2) %>% 
              distinct(acceptedName, .keep_all = TRUE) %>% 
              mutate(synonym = acceptedName)) %>% 
  mutate(acceptedNameParent = ifelse(word_count > 2,
                                     stringr::word(acceptedName, 1, 2),
                                     acceptedName)
  )
# Filter name relations in itis down to North American plants 
itis_df <- coal_df %>% 
  filter(source == "itis")
itis_df <- itis_df[, itisTaxonomicStatus := ifelse(nameMatch == nameAligned, "ACCEPTED", "SYNONYM")]
itis_df_accepted <- itis_df[itisTaxonomicStatus == "ACCEPTED"]
itis_df_synonym <- itis_df[itisTaxonomicStatus == "SYNONYM"]
itis_relations_na <- itis_relations %>% 
  filter(acceptedName %in% itis_df_accepted$nameMatch | synonym %in% itis_df_synonym$nameMatch) %>% 
  select(-word_count) %>% 
  rename(taxonID = acceptedNameUsageID) %>% 
  select(taxonID, acceptedName, acceptedNameParent, synonym, catalogID)

### GBIF 
gbif <- taxa_tbl("gbif") %>% 
  select(scientificName,taxonRank,acceptedNameUsageID,taxonomicStatus) %>%
  filter(taxonRank == "species")

# Load into memory, convert to data.table for faster operations
gbif_names <- gbif %>% collect() # load into memory
gbif_names <- as.data.table(gbif_names)

# Organize by taxonomic status
gbif_accepted_names <- gbif_names[taxonomicStatus == "accepted"]
gbif_synonymous_names <- gbif_names[taxonomicStatus == "synonym"]

# Check for ambigious mapping in gbif on accepted names
gbif_non_uniq_accepted_names <- gbif_accepted_names %>% 
  group_by(scientificName) %>% 
  mutate(multipleAcceptableUsageIDs = n_distinct(acceptedNameUsageID)) %>% 
  filter(multipleAcceptableUsageIDs >= 2)

# There are non-existent acceptedNameUsageIDs cluttering in the synonyms, 
# first remove these, then check for multiple cases
gbif_accepted_usage_ids <- unique(gbif_accepted_names$acceptedNameUsageID)
gbif_synonymous_names <- gbif_synonymous_names[acceptedNameUsageID %in% gbif_accepted_usage_ids]
gbif_non_uniq_synonymous_names <- gbif_synonymous_names %>% 
  group_by(scientificName) %>% 
  mutate(multipleAcceptableUsageIDs = n_distinct(acceptedNameUsageID)) %>% 
  filter(multipleAcceptableUsageIDs >= 2)

# these non-unique cases are true issues likely related to authorship. We're going to choose to drop them
gbif_accepted_names_drop <- gbif_accepted_names[!scientificName %in% gbif_non_uniq_accepted_names$scientificName]
gbif_synonymous_names_drop <- gbif_synonymous_names[!scientificName %in% gbif_non_uniq_synonymous_names$scientificName]

# create relational table for gbif taxonomy 
gbif_relations <- merge(gbif_accepted_names_drop, gbif_synonymous_names_drop, by = "acceptedNameUsageID", all.x = TRUE)
gbif_relations <- gbif_relations %>% 
  rename(acceptedName = scientificName.x,
         synonym = scientificName.y) %>% 
  mutate(catalogID = "gbif")

## Deal with infraspecifics in gbif
gbif_relations <- gbif_relations %>% 
  mutate(word_count = str_count(acceptedName, "\\s+") + 1) 

# Duplicate these infraspecific instances into the synonymy so we know that they exist
gbif_relations <- gbif_relations %>% 
  rename(accepted_name_taxon_rank = taxonRank.x,
         synonym_name_taxon_rank = taxonRank.y,
         accepted_taxon_status = taxonomicStatus.x,
         synonym_taxon_status = taxonomicStatus.y) %>% 
  bind_rows(gbif_relations %>% 
              filter(word_count > 2) %>% 
              distinct(acceptedName, .keep_all = TRUE) %>% 
              mutate(synonym = acceptedName)) %>% 
  mutate(acceptedNameParent = ifelse(word_count > 2,
                                     stringr::word(acceptedName, 1, 2),
                                     acceptedName)
  )
# Filter name relations in gbif down to North American plants 
gbif_df <- coal_df %>% 
  filter(source == "gbif")
gbif_df <- gbif_df[, gbifTaxonomicStatus := ifelse(nameMatch == nameAligned, "ACCEPTED", "SYNONYM")]
gbif_df_accepted <- gbif_df[gbifTaxonomicStatus == "ACCEPTED"]
gbif_df_synonym <- gbif_df[gbifTaxonomicStatus == "SYNONYM"]
gbif_relations_na <- gbif_relations %>% 
  filter(acceptedName %in% gbif_df_accepted$nameMatch | synonym %in% gbif_df_synonym$nameMatch) %>% 
  select(-word_count) %>% 
  rename(taxonID = acceptedNameUsageID) %>% 
  select(taxonID, acceptedName, acceptedNameParent, synonym, catalogID)

### Combine North American Relational Tables
# Combine the relational datatables
relational_tb <- bind_rows(wcvp_relations_na, col_relations_na, itis_relations_na, gbif_relations_na)
# Clean up acceptedNames that do not have associated synonyms, basically just adds the acceptedName as a synonym as well for clarity
relational_tb <- relational_tb %>% 
  mutate(synonym = case_when(
    is.na(synonym) == TRUE ~ acceptedName, 
    is.na(synonym) == FALSE ~ synonym
  ))
### Deal with duplications across taxonomic backbones
# Write a fxn that flags duplicates for acceptedNames and synonyms
add_duplication_flags <- function(df) {
  df <- df %>%
    group_by(acceptedName) %>%
    mutate(a_dupes = n_distinct(catalogID)) %>%
    ungroup() %>%
    group_by(synonym) %>%
    mutate(s_breaks = n_distinct(catalogID)) %>%
    ungroup()
  return(df)
}
relational_tb_wflags <- add_duplication_flags(relational_tb)
priority_order <- c("wcvp", "col", "itis", "gbif")
# Filter by priority to remove duplicate names 
relational_tb_filtered <- relational_tb_wflags %>% 
  mutate(priority = case_when(
    a_dupes > 1 | s_breaks > 1  ~ match(catalogID, priority_order), 
    TRUE ~ 0)) %>% 
  group_by(acceptedName) %>% 
  filter(priority == min(priority))  %>% 
  ungroup() %>% 
  group_by(synonym) %>% 
  filter(priority == min(priority)) %>% 
  ungroup()

### Append Information for Multiple Mapping Statuses to this relational table 
# Merge acceptedNames and Synonyms with the multiple mapping information we built prior 
accepted_matched_names_wcvp <- filter(matched_names_wcvp, taxon_status == "Accepted")
accepted_half_merge <- merge(relational_tb_filtered, accepted_matched_names_wcvp, by.x = "accepted_powo_id", by.y = "powo_id", all.x = TRUE)
synonym_matched_names_wcvp <- filter(matched_names_wcvp, taxon_status != "Accepted")
synonym_half_merge <- merge(relational_tb_filtered, synonym_matched_names_wcvp, by.x = "synonym_powo_id", by.y = "powo_id", all.y = TRUE) 
merged_df <- rbind(accepted_half_merge, synonym_half_merge)

# Assure distinct values only present after all those merges
merged_d_df <- distinct(merged_df, accepted_plant_name_id, acceptedName, acceptedNameParent, accepted_powo_id, 
                        synonym, synonym_powo_id, synonym_taxon_authors, catalogID, accepted_taxon_authors, authorship,
                        multipleMappingsPossible, spacelessOurAuthorship, spacelessWCVPAuthorship, 
                        authorshipMatch, multMapAuthorshipMatch, multMapResolutionPossible, authorshipInAcceptedName)

# clean up naming on cols
merged_dc_df <- merged_d_df %>% 
  mutate(acceptedNameScientificNameAuthorship = case_when( # added this casewhen conditional to deal with some of the complexities concerning when things are actually multiple mapping and their accepted authorship.
    authorshipInAcceptedName == TRUE & !is.na(authorship) ~ authorship, 
    authorshipInAcceptedName == FALSE ~ accepted_taxon_authors, 
    is.na(authorshipInAcceptedName) ~ accepted_taxon_authors,
    authorshipInAcceptedName == TRUE & is.na(authorship) ~ accepted_taxon_authors)) %>% 
  rename(acceptedNameMultMapsPossible = multipleMappingsPossible,
         acceptedNameSpacelessWCVPAuthorship = spacelessWCVPAuthorship, 
         acceptedNameAuthorshipMatch = authorshipMatch,
         acceptedNameMultMapResolutionPossible =  multMapResolutionPossible, 
         synonymNameWCVPAuthorship = synonym_taxon_authors) %>% 
  mutate(acceptedNameScientificNameAuthorship = gsub(" ", "", acceptedNameScientificNameAuthorship)) %>% 
  mutate(acceptedNameScientificNameAuthorship = tolower(acceptedNameScientificNameAuthorship)) 

### Denote whether multiple mappings are possible for Synonyms  
multiple_s_mappings <- merged_dc_df %>% 
  group_by(synonym) %>% 
  mutate(numDiffAcceptedNames = length(unique((acceptedName))))
syn_mult_maps <- filter(merged_dc_df, multiple_s_mappings$numDiffAcceptedNames > 1)
merged_dc_df$synonymMultMapPossible <- ifelse(merged_dc_df$synonym %in% syn_mult_maps$synonym, TRUE, FALSE)

### Some other housekeeping 
# remove caps and spaces to synonym authorship to match the format of our accepted name authors
merged_dc_df <- merged_dc_df %>% 
  mutate(synonymNameSpacelessWCVPAuthorship = tolower(gsub(" ", "", synonymNameWCVPAuthorship))) %>% 
  select(-accepted_taxon_authors, -authorship, -spacelessOurAuthorship, -multMapAuthorshipMatch, -acceptedNameMultMapResolutionPossible, -authorshipInAcceptedName)

final_df <- merged_dc_df %>% 
  filter(!acceptedName == "Senna × floribunda")

### Report 
final_df %>% 
  mutate(copy = case_when(
    acceptedName == synonym ~ TRUE,
    acceptedName != synonym ~ FALSE
  )) %>% 
  filter(copy == FALSE) %>% 
  distinct(synonym, synonymNameWCVPAuthorship) %>% 
  nrow()

final_df %>% 
  summarize(acceptedNameLength = length(unique(acceptedName)), 
            acceptedNameParentLength = length(unique(acceptedNameParent)))

fwrite(final_df, "/home/jt-miller/Gurlab/BoCP/data/processed/finalized_name_alignment_wcvp.csv")

