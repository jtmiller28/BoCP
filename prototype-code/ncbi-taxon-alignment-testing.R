### Testing taxonomy alignment to NCBI

# Load Packages 
library(ape)
library(data.table)
library(taxadb)
library(tidyverse)
library(purrr)
library(rgnparser)

# Load data
raw_bocp_names <- fread("./data/raw/USACANMEX_plantNamesFieldsAdded.csv",encoding = "UTF-8")
full_tree <- ape::read.tree("./data/processed/tree-outputs/RAxML_bestTree.POAOUT_rax.rr.pr")

# plot intial tree
plot.phylo(full_tree)

# load in ncbi's taxonomic db
ncbi_names <- taxa_tbl("ncbi") # load in ncbi's names according to taxadb's last pull
ncbi_names <- ncbi_names %>%  collect()
ncbi_names <- ncbi_names %>% 
  mutate(acceptedNameUsageIDLess = gsub("NCBI:", "", acceptedNameUsageID)) # standardize with the notation is Stephen's tree

# Standardize names in the North American list #################################################################################################################################################################################################################
# Set up pathing to use gnparser
my_path <- Sys.getenv("PATH") # grab our path
Sys.setenv(PATH = paste0(my_path, "/home/millerjared/gnparser"))
## Collapse all mentions of the author into the species to then let gnparser take care of seperation.
raw_bocp_names$concatenatedName <- apply(raw_bocp_names[, c("Name", "ExteriorNameAuthorship", "ExteriorNameAuthorship2", "V5")], 1, function(x) paste(x[x != ""], collapse = ", ")) # Collapse fields into one by adding a space when there is a string present
bocp_names <- raw_bocp_names[, .(TropicosNameID = NameID, concatenatedName)] # Select only fields of interest 
bocp_parsed_names <- gn_parse_tidy(bocp_names$concatenatedName)  # break up names
bocp_parsed_names_df <- merge(bocp_names, bocp_parsed_names, by.x = 'concatenatedName', by.y = 'verbatim')
## Identify special cases in the data 
copied_names <- bocp_parsed_names_df %>% # Find the cases of non-unique names 
  group_by(concatenatedName) %>%
  mutate(n = n()) %>%
  filter(n > 1)
## Remove hybrids 
hybrids <- bocp_parsed_names_df %>% # Grab a list of names that are hybrids
  filter(grepl("×", concatenatedName)) 
bocp_parsed_names_df <- bocp_parsed_names_df %>%
  filter(!concatenatedName %in% hybrids$concatenatedName)
################################################################################################################################################################################################################################################################

# Join these names to NCBI's taxonomy to infer where they sit ##################################################################################################################################################################################################
## Parse ncbi's taxonomy because they set up their db oddly (specificEpithet is the translated name, not verbatim)
ncbi_parsed_names <- gn_parse_tidy(ncbi_names$scientificName)
ncbi_parsed_names_sl <- select(ncbi_parsed_names, canonicalfull, verbatim)
ncbi_parsed_names_f <- filter(ncbi_parsed_names_sl, !is.na(canonicalfull)) # removal of things that are not parsable. 
ncbi_names2 <- merge(ncbi_names, ncbi_parsed_names_f, by.x = "scientificName", by.y = "verbatim")
## Do a left-join with our parent list and NCBI's taxonomy 
matched_names <- merge(bocp_parsed_names_df, ncbi_names2, by.x = "canonicalfull", by.y = "canonicalfull", all.x = TRUE)
## Remove Metazoa kindgom as those are not plants
matched_names <- filter(matched_names, kingdom == "Viridiplantae" | is.na(kingdom)) # retain NAs for now as Im not sure what these are...
## Select and Rename for standardization
matched_names <- matched_names %>% 
  select(nameMatch = canonicalfull, concatenatedName, TropicosNameID, authorship, scientificName, ncbiTaxonRank = taxonRank, 
         ncbiTaxonomicStatus = taxonomicStatus, ncbiAcceptedNameUsageID = acceptedNameUsageID, 
         ncbiAcceptedNameUsageIDLess = acceptedNameUsageIDLess, ncbiKingdom = kingdom, ncbiPhylum = phylum,
         ncbiClass = class, ncbiOrder = order, ncbiFamily = family, ncbiGenus = genus, ncbiMappedSpecies = specificEpithet)
################################################################################################################################################################################################################################################################

# Summaries + Multiple Mapping situation ############################################################################################################################################################################################################################
matched_names_m <- matched_names %>%
  group_by(nameMatch) %>%
  mutate(num_different_paths = n_distinct(ncbiAcceptedNameUsageIDLess))

matched_names_m <- matched_names %>%
  group_by(nameMatch) %>%
  mutate(num_different_paths = n_distinct(ncbiAcceptedNameUsageIDLess))

multiple_mappings <- matched_names_m %>%
  filter(num_different_paths > 1) 

matched_names$multipleMappingsPossible <- ifelse(matched_names$nameMatch %in% multiple_mappings$nameMatch, TRUE, FALSE)

summary_report <- matched_names %>% 
  filter(multipleMappingsPossible == TRUE) %>% 
  distinct(ncbiMappedSpecies, .keep_all = TRUE)

# Testing removal of erroneous names: cf stands for uncertain designations, which shouldnt be in our alignment. Remove these. Numeric values in the name are not helpful/useful so removal of these as well. 
# Filter rows without " cf. " or numeric values
filtered_names <- summary_report[!grepl(" cf\\. ", summary_report$ncbiMappedSpecies) & !grepl("\\d", summary_report$ncbiMappedSpecies), ]

true_m_mappings <- filtered_names %>% 
  group_by(nameMatch) %>% 
  mutate(n = n()) %>% 
  arrange(desc(n))

report_true_m_mappings <- filter(true_m_mappings, n > 1)

fwrite(true_m_mappings, "./data/processed/report_true_ncbi_m_maps_plants.csv")

# Remove names that have true multiple mappings
matched_names <- filter(matched_names, !ncbiMappedSpecies %in% true_m_mappings$ncbiMappedSpecies)
# Remove names with associated uncertainty designations + just numeric nonesense. 
unc_names <- summary_report[grepl(" cf\\. ", summary_report$ncbiMappedSpecies) & grepl("\\d", summary_report$ncbiMappedSpecies), ]
matched_names <- matched_names[!ncbiMappedSpecies %in% unc_names$ncbiMappedSpecies]
################################################################################################################################################################################################################################################################

### There are about 20% of the names that did not obtain a mapping using NCBI taxonomy, to at least address some proportion of them, we will try using wcvp to see if any of them are synoynms of more well known names to the ncbi db #########################
unmatched_names <- matched_names[is.na(ncbiMappedSpecies)] # filter for these uncertain species
unmatched_names <- select(unmatched_names, canonicalfull = nameMatch, concatenatedName, TropicosNameID, authorship)
wcvp_backbone <- fread("./data/raw/wcvp_names.csv") # load wcvp's db 

# check to see if we can find an alignment within POWO/WCVP
rematched_names <- merge(unmatched_names, wcvp_backbone, by.x = "canonicalfull", by.y = "taxon_name", all.x = TRUE)

rematched_names %>% 
  filter(is.na(taxon_status)) %>% 
  nrow() # 207 names cannot find a match.

rematched_names <- rematched_names %>% # Simplify to select fields of interest
  select(nameMatch = canonicalfull, concatenatedName, TropicosNameID, authorship, 
         wcvp_accepted_plant_name_id = accepted_plant_name_id, taxon_rank, wcvp_taxon_authors = taxon_authors, 
         family, genus, specificEpithet = species, taxon_status, reviewed, powo_id) 

# set aside NAs (actually skip this, messes with things in the future)
# names_not_in_NCBI_or_POWO <- filter(rematched_names, is.na(taxon_status))
# rematched_names <- filter(rematched_names, !nameMatch %in% names_not_in_NCBI_or_POWO$nameMatch)

## Check for multiple mapping names within the wcvp/powo backbone so we can isolate these names for further resolution via authorship
wcvp_backbone_m <- wcvp_backbone %>%
  group_by(taxon_name) %>%
  mutate(num_different_paths = n_distinct(accepted_plant_name_id)) # adjusting this to match wcvp 

multiple_mappings <- wcvp_backbone_m %>%
  filter(num_different_paths > 1) # it would seem that most of these cases are actually subspecific. 

# attach a boolean value to denote whether names in this subset of found names contain multiple mapping matches
rematched_names$multipleMappingsPossible <- ifelse(rematched_names$nameMatch %in% multiple_mappings$taxon_name, TRUE, FALSE)

# check to see how many names fall under this criteria for these rematched subset of nams
matched_name_summary <- rematched_names %>%
  group_by(taxon_status, multipleMappingsPossible) %>%
  count() %>%
  rename(count = n)

## Identify 'recoverable' names based on authorship present in our dataset 
### rm spaces & capitilization, then perform an exact match 
rematched_names <- rematched_names[, c("spacelessOurAuthorship", "spacelessWCVPAuthorship") := 
                                 lapply(.SD, function(x) tolower(gsub(" ", "", x))), 
                               .SDcols = c("authorship", "wcvp_taxon_authors")]
rematched_names$authorshipMatch <- rematched_names$spacelessOurAuthorship == rematched_names$spacelessWCVPAuthorship

## Create new a new field called multMapAuthorship match that requires authorship to match when multipleMappings is Possible
rematched_names <- rematched_names %>% 
  mutate(multMapAuthorshipMatch = case_when(
    multipleMappingsPossible & authorshipMatch == TRUE ~ TRUE, 
    multipleMappingsPossible & authorshipMatch != TRUE ~ FALSE))

## Check to make sure that mapping is 1:1 while accounting for authorship in multiple mappings.
resolvable_names_df <- rematched_names %>% 
  filter(multipleMappingsPossible == TRUE) %>% 
  filter(multMapAuthorshipMatch == TRUE) %>% 
  group_by(nameMatch) %>% 
  mutate(numDiffTaxStatsWithSameAuthor = n_distinct(taxon_status)) %>% # num of distinct paths
  ungroup() %>% 
  filter(numDiffTaxStatsWithSameAuthor == 1) # Filter to the ones that are usable 

## Create a field called multMapResolutionPossible that denotes whether in case there is multiple maps if AuthorshipMatch is possible, else other cases
rematched_names <- rematched_names %>% 
  mutate(multMapResolutionPossible = case_when(
    nameMatch %in% resolvable_names_df$nameMatch & multMapAuthorshipMatch == TRUE ~ TRUE,
    multipleMappingsPossible == FALSE ~ NA, # if we dont need this leave as NA
    TRUE ~ FALSE # else make it false...
  ))

## Report number of names recovered using Authorship as mapping 
rematched_names %>% 
  filter(multMapResolutionPossible) %>% 
  distinct(nameMatch) %>% 
  nrow()

## Remove the multiple matching names that cannot be resolved via authorship
rematched_names <- rematched_names %>% 
  filter(multMapResolutionPossible == TRUE & multipleMappingsPossible == TRUE | multipleMappingsPossible == FALSE) %>% 
  mutate(authorshipInAcceptedName = ifelse(taxon_status == "Accepted", TRUE, FALSE))

## Create a summary of alignment status 
alignment_summary <- rematched_names %>% 
  group_by(taxon_status) %>% 
  summarize(n = n()) 
print(alignment_summary)

matched_names_wcvp <- rematched_names %>% 
  mutate(source = "wcvp") %>% 
  filter(!taxon_status == "Unplaced" & !is.na(taxon_status))

### Split the dataset. Taxa that have any form of alignment outside of Unplaced and NA are ok for wcvp.
## Unplaced names will be thrown away as the creators of wcvp have issues with these names 
## NA names will be fed through the catalogs COL, ITIS, and GBIF in order to assign aligmnet 
rematched_names_wcvp <- rematched_names %>% 
  mutate(source = "wcvp") %>% 
  filter(!taxon_status == "Unplaced" & !is.na(taxon_status))

names_not_in_NCBI_or_POWO <- filter(rematched_names, is.na(taxon_status)) 

## For unmatched names, harmonize against COL, ITIS, and GBIF backbones
alternative_catalog_matches <- names_not_in_NCBI_or_POWO %>%  
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

################################################################################################################################################################################################################################################################





# join and format the tip.label 
full_tree$tip.label <- map_chr(full_tree$tip.label, ~ {
  match_row <- df %>% filter(acceptedNameUsageIDLess == .x)
  
  if(nrow(match_row) > 0){
    paste0(match_row$scientificName, ":", .x)
  } else{
    .x
  }
})
