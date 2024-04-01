# 001-taxon-name-alignment Rewrite
# Author: JT Miller 
# Date: 03-27-2024
# Project: BoCP 

## Load Dependencies 
library(data.table)
library(tidyverse)
library(rgnparser)

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







