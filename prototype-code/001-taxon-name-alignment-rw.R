# 001-taxon-name-alignment Rewrite
# Author: JT Miller 
# Date: 03-27-2024
# Project: BoCP 

## Load Dependencies 
library(data.table)
library(tidyverse)
library(rgnparser)
library(taxadb)

## Load Data
raw_bocp_names <- fread("/home/jt-miller/Gurlab/BoCP/data/raw/USACANMEX_plantNamesFieldsAdded.csv",encoding = "UTF-8")
wfo_backbone <- fread("/home/jt-miller/Soltis-lab/taxonomic-harmonization-resources/taxonomic-backbones/WFO_Backbone/classification.csv")

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

### Align with WorldFloraOnline Taxonomy 

## Using the wfo backbone, merge on the canonicalFull name 
matched_names <- merge(bocp_parsed_names_df2, wfo_backbone, by.x = "canonicalfull", by.y = "scientificName", all.x = TRUE)
matched_names <- matched_names %>% # Simplify to select fields of interest
  select(nameMatch = canonicalfull, concatenatedName, TropicosNameID, authorship, wfoTaxonID = taxonID, taxonRank, scientificNameAuthorship, family, genus, specificEpithet, taxonomicStatus, references) 

## Identify Multiple Mapping Cases based on underlying Authorship
wfo_backbone_m <- wfo_backbone %>%
  group_by(scientificName) %>%
  mutate(num_different_paths = n_distinct(taxonomicStatus))

multiple_mappings <- wfo_backbone_m %>%
  filter(num_different_paths > 1)

matched_names$multipleMappingsPossible <- ifelse(matched_names$nameMatch %in% multiple_mappings$scientificName, TRUE, FALSE)

matched_name_summary <- matched_names %>%
  group_by(taxonomicStatus, multipleMappingsPossible) %>%
  count() %>%
  rename(count = n)

ggplot(matched_name_summary, aes(x = taxonomicStatus, y = count, fill = multipleMappingsPossible)) +
  geom_col(width = 0.5) +
  geom_text(aes(label = count), vjust = -0.5, position = position_stack(vjust = 1.0)) +
  scale_fill_manual(values = c("TRUE" = "darkred", "FALSE" = "darkgrey"), name = "Multiple Name\nStatus Based On\n Authority") +  # Set colors for TRUE and FALSE
  theme_bw() +
  ggtitle("Taxon Status of North American Plant Name List According to the WorldFloraOnline Backbone") +
  theme(legend.title = element_text(hjust = 0.5),
        legend.text.align = 0.5)

## Identify 'recoverable' names based on authorship present in our dataset 
### rm spaces, then preform an exact match 
# Remove spacing in all strings in 'column1' and 'column2'
matched_names <- matched_names[, c("spacelessOurAuthorship", "spacelessWFOAuthorship") := lapply(.SD, function(x) gsub(" ", "", x)), .SDcols = c("authorship", "scientificNameAuthorship")]
matched_names$authorshipMatch <- matched_names$spacelessOurAuthorship == matched_names$spacelessWFOAuthorship

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
  mutate(numDiffTaxStatsWithSameAuthor = n_distinct(taxonomicStatus)) %>% # num of distinct paths
  ungroup() %>% 
  filter(numDiffTaxStatsWithSameAuthor == 1) # Filter to the ones that are usable 

matched_names4 <- matched_names3 %>% 
  mutate(multMapResolutionPossible = case_when(
    nameMatch %in% resolved_df$nameMatch & multMapAuthorshipMatch == TRUE ~ TRUE,
    multipleMappingsPossible == FALSE ~ NA, # if we dont need this leave as NA
    TRUE ~ FALSE # else make it false...
  ))

# Remove those redundant names
matched_names5 <- matched_names4 %>% 
  ungroup() %>% 
  filter(multMapResolutionPossible == TRUE & multipleMappingsPossible != FALSE | multipleMappingsPossible == FALSE)

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
  group_by(taxonomicStatus) %>% 
  summarize(n = n())

# Split datasets based on whether we can get a resolution via WFO 
matched_names_wfo <- matched_names6 %>% 
  filter(taxonomicStatus %in% c("Accepted", "Synonym")) %>% 
  mutate(source = "wfo") %>% 
  rename(taxonomicSourceID = wfoTaxonID)

unmatched_names <- matched_names6 %>% 
  filter(taxonomicStatus == "Unchecked" | is.na(taxonomicStatus))

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

# Keep aligned name, append to wfo harmonized names 
other_catalog_matched_names <- coal_df %>% 
  select(-col_id, -itis_id, -gbif_id, -col_harmonizedName, -itis_harmonizedName, -gbif_harmonizedName) %>% 
  rename(taxonomicStatus = wfoTaxonomicStatus)

matched_names_wfo <- matched_names_wfo %>% 
  mutate(nameAligned = nameMatch)

full_match_df <- rbind(matched_names_wfo, other_catalog_matched_names)

# Remove records that lack a source of resolution
full_match_df <- full_match_df %>% 
  filter(!is.na(source))

### Summarize the taxon name alignment 
# Check to see how many names align per taxonomic source, plus how many names will not align
full_match_df %>% 
  group_by(source) %>% 
  summarize(n = n())
# Check to see how many multiple mapping names we are capable of recovering 
full_match_df %>% 
  filter(multipleMappingsPossible == TRUE) %>% 
  group_by(multMapResolutionPossible) %>% 
  summarize(recoveredMultMaps = n())

# create a variable for storing these failed to resolve multiple mapping names 
mult_mapping_check <- matched_names4 %>% 
  filter(multipleMappingsPossible) %>% 
  group_by(nameMatch, multMapResolutionPossible)

mult_mapping_check <- mult_mapping_check %>% 
  group_by(nameMatch) %>% 
  mutate(min_resolved = any(multMapResolutionPossible)) %>% 
  ungroup()

mult_mapping_passes <- mult_mapping_check %>% 
  filter(min_resolved == TRUE) %>% 
  distinct(nameMatch, .keep_all = TRUE)

mult_mapping_fails <- mult_mapping_check %>% 
  filter(min_resolved == FALSE) %>% 
  distinct(nameMatch, .keep_all = TRUE)

nrow(mult_mapping_passes) / (nrow(mult_mapping_passes) + nrow(mult_mapping_fails))


# Check how many names we lost total 
name_beg <- length(unique(bocp_parsed_names_df2$canonicalfull))
name_end <- length(unique(full_match_df$nameMatch))

name_beg - name_end # Expected, these are the multiple mappings that were impossible to resolve with current methods. 

### Attach Synonyms to our sets of names


# Split according to required authority 
wfo_names <- full_match_df[source == "wfo"]
wfo_accepted_names <- wfo_names[taxonomicStatus == "Accepted"]
wfo_synonym_names <- wfo_names[taxonomicStatus == "Synonym"]

col_names <- full_match_df[source == "col"]
full_match_df2 <- full_match_df %>%
  select(!taxonomicStatus)
itis_names <- full_match_df[source == "itis"]
gbif_names <- full_match_df[source == "gbif"]

# Build synonym lists based on authority 
wfo_catalogue_accepted <- wfo_names[t]





### Create a dataframe that contains the accepted name and its associated synoynms 
# rebuild the taxonomic field for Accepted names and Synonym names as they are still currently only aligned according to WFO

full_match_df3 <- full_match_df2 %>% 
  mutate(taxonomicStatus = case_when(
    nameMatch == nameAligned ~ "Accepted", 
    nameMatch != nameAligned ~ "Synonym"
  ))

accepted_df <- full_match_df3 %>% 
  filter(taxonomicStatus == "Accepted")
synonym_df <- full_match_df3 %>% 
  filter(taxonomicStatus == "Synonym")
accepted_df <- full_match_df3[taxonomicStatus == "Accepted"] # Create an accepted name table
synonym_df <- full_match_df3[taxonomicStatus == "Synonym"] # Create a synonym list table 
