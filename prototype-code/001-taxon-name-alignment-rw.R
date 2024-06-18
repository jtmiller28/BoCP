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

# Create some simple summary reports 
num_of_instances_w_multiple_maps_possible <- matched_names %>% 
  filter(multipleMappingsPossible == TRUE) %>% 
  nrow()
num_distinct_names_w_multiple_maps_possible <- matched_names %>% 
  filter(multipleMappingsPossible == TRUE) %>% 
  distinct(nameMatch) %>%  
  nrow()


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

# Report number of names recovered using Authorship as mapping 
mult_map_recovery_name_num <- matched_names4 %>% 
  filter(multMapResolutionPossible) %>% 
  distinct(nameMatch) %>% 
  nrow()

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

  
  
### Attach Synonyms: WFO and other catalogues
## 6/18/2024 edit: Add in subspecific designations to parent taxon when they are accepted names within WFO's taxonomy
# First find the accepted subspecific assignments
wfo_backbone_subspecific <- wfo_backbone %>% 
  filter(taxonRank %in% c("variety", "subspecies", "form", "subvariety", "subform")) %>% 
  filter(taxonomicStatus == "Accepted")
# Find what their parent taxon is 
wfo_backbone_subspecific <- wfo_backbone_subspecific %>% 
  mutate(parentTaxon = stringr::word(scientificName, 1, 2)) %>% 
  mutate(taxonomicStatus = "Synonym") # manually change these to Synonyms.

updated_wfo_backbone <- merge(wfo_backbone, wfo_backbone_subspecific, by.x = "scientificName", by.y = "parentTaxon")
updated_wfo_backbone <- updated_wfo_backbone %>% 
  mutate(statusMatch = ifelse(taxonomicStatus.x == taxonomicStatus.y, TRUE, FALSE)) %>% 
  filter(statusMatch == FALSE) %>%  # grab the cases we designed
  mutate(acceptedNameUsageID.y = ifelse(acceptedNameUsageID.y == "", taxonID.x, acceptedNameUsageID.y)) %>%  
  rename(acceptedNameUsageID = acceptedNameUsageID.y, parentTaxon = scientificName, scientificName = scientificName.y, 
         scientificNameAuthorship = scientificNameAuthorship.y, taxonomicStatus = taxonomicStatus.y,
         taxonID = taxonID.x) %>% 
  select(taxonID, scientificName, scientificNameAuthorship, parentTaxon, taxonomicStatus,  acceptedNameUsageID, taxonID) 

test_mult_maps <- updated_wfo_backbone %>%
  group_by(scientificName) %>%
  mutate(num_different_paths = n_distinct(taxonomicStatus)) 

test_multiple_mappings <- test_mult_maps  %>% # tis all good! 
  filter(num_different_paths > 1)

wfo_backbone2 <- wfo_backbone %>% 
  filter(!scientificName %in% updated_wfo_backbone$scientificName)

wfo_backbone2 <- wfo_backbone2 %>% 
  mutate(parentTaxon = NA) %>% 
  select(taxonID, scientificName, scientificNameAuthorship, parentTaxon, taxonomicStatus,  acceptedNameUsageID) %>% 
  rbind(updated_wfo_backbone)
# WFO 
wfo_accepted_mapping<- wfo_backbone2[taxonomicStatus == "Accepted"]
wfo_synonym_mapping <- wfo_backbone2[taxonomicStatus == "Synonym"]
non_unique_accepted_names <- wfo_accepted_mapping[, .N, by = "scientificName"]# Identify non unique names 
# Create a simplified version of each respective dt, then attach synonyms
accepted_mapping <- wfo_accepted_mapping[, c("taxonID", "scientificName"), with = FALSE]
synonym_mapping <- wfo_synonym_mapping[, c("scientificName", "acceptedNameUsageID"), with = FALSE]

# Append acceptedNameUsageID with that of its parent taxon. wf
# Subset the backbone for our needs 
wfo_df_accepted <- matched_names_wfo[taxonomicStatus == "Accepted"]
wfo_df_synonym <- matched_names_wfo[taxonomicStatus == "Synonym"]
wfo_relations <- merge(accepted_mapping, synonym_mapping, by.x = "taxonID", by.y = "acceptedNameUsageID", all.x = TRUE)
wfo_relations <- wfo_relations %>% 
  rename(acceptedName = scientificName.x, synonyms = scientificName.y) 
wfo_relations_f <- wfo_relations %>%  
  filter(acceptedName %in% wfo_df_accepted$nameMatch | synonyms %in% wfo_df_synonym$nameMatch) %>% 
  mutate(catalogID = "wfo") %>% 
  distinct(taxonID, acceptedName, synonyms, .keep_all = TRUE)
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
relational_tb <- bind_rows(wfo_relations_f, col_relations_f, itis_relations_f, gbif_relations_f )

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
priority_order <- c("wfo", "col", "itis", "gbif")

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

# There are cases where there are repeated synonyms for the same accepted name, this is presumably due to synonyms having multiple mapping paths depending on their underlying authorship. Lets resolve this. 
final_df2 <- merge(final_df, matched_names_wfo, by.x = "taxonID", by.y = "taxonomicSourceID", all.x = TRUE)
final_df2 <- select(final_df2, names(final_df), scientificNameAuthorship, multipleMappingsPossible, spacelessWFOAuthorship, authorshipMatch, multMapResolutionPossible)
final_df3 <- final_df2 %>% # rename for clarity that our multiple mapping resolution was depedent on the acceptedNames, now we need to deal with the synonyms
  rename(acceptedNameScientificNameAuthorship = scientificNameAuthorship,
         acceptedNameMultMapsPossible = multipleMappingsPossible,
         acceptedNameSpacelessWFOAuthorship = spacelessWFOAuthorship, 
         acceptedNameAuthorshipMatch = authorshipMatch,
         acceptedNameMultMapResolutionPossible =  multMapResolutionPossible)

# merge synonymy, create spaceless authorship
wfo_backbone_s <- wfo_backbone %>% 
  filter(scientificName %in% final_df3$synonyms) %>% 
  group_by(scientificName) %>% 
  mutate(numDiffAcceptedPaths = n_distinct(acceptedNameUsageID)) %>% 
  ungroup()

multiple_s_mappings <- wfo_backbone_s %>% 
  filter(numDiffAcceptedPaths > 1)

multiple_s_mappings <- multiple_s_mappings %>%  
  mutate(synonymNameSpacelessWFOAuthorship = gsub(" ", "", scientificNameAuthorship)) # remove spacing 

  
final_df4 <- final_df3 %>%
  mutate(synonymMultipleMapsPossible = case_when(
    catalogID != "wfo" ~ NA,
    synonyms %in% multiple_s_mappings$scientificName ~ TRUE,
    TRUE ~ FALSE
  ))

multiple_s_mappings_fields <- multiple_s_mappings %>% 
  select(taxonID, scientificName, synonymNameSpacelessWFOAuthorship, acceptedNameUsageID) %>% 
  rename(synonym = scientificName) %>% 
  mutate(source = "df1")

final_df5 <- final_df4 %>% 
  rename(acceptedNameUsageID = taxonID) %>% 
  rename(synonym = synonyms) %>% 
  mutate(source = "df2")

final_df6 <- merge(final_df5, multiple_s_mappings_fields, by = c("synonym", "acceptedNameUsageID"), all.x = TRUE)

final_df7 <- final_df6 %>%  select(acceptedName, synonym, acceptedNameUsageID, catalogID, a_dupes, s_breaks, priority, acceptedNameSpacelessWFOAuthorship, synonymNameSpacelessWFOAuthorship, acceptedNameMultMapResolutionPossible, synonymMultipleMapsPossible, source.x, source.y)

final_df8 <- final_df7 %>% 
  distinct(acceptedName, synonym, acceptedNameUsageID, catalogID, a_dupes, s_breaks, priority, acceptedNameSpacelessWFOAuthorship, synonymNameSpacelessWFOAuthorship, acceptedNameMultMapResolutionPossible, synonymMultipleMapsPossible,.keep_all = TRUE)

length(unique(final_df8$acceptedName))

final_synonym_count <- final_df8 %>% 
  mutate(copy = case_when(
    acceptedName == synonym ~ TRUE,
    acceptedName != synonym ~ FALSE
  )) %>% 
  filter(copy == FALSE) %>% 
  distinct(synonym, synonymNameSpacelessWFOAuthorship) %>% 
  nrow()

final_multiple_map_on_synonym_count <- final_df8 %>% 
  mutate(copy = case_when(
    acceptedName == synonym ~ TRUE,
    acceptedName != synonym ~ FALSE
  )) %>% 
  filter(copy == FALSE & synonymMultipleMapsPossible == TRUE) %>% 
  distinct(synonym, synonymNameSpacelessWFOAuthorship) %>% 
  nrow()

final_df9 <- final_df8 %>% 
  select(-source.x, -source.y, -a_dupes, -s_breaks, -priority)


fwrite(final_df9, "/home/jt-miller/Gurlab/BoCP/data/processed/finalized_name_alignment.csv")

  
  
  














































