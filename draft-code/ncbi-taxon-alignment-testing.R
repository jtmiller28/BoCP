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
matched_names <- matched_names %>% rename(original_authorship = authorship)
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

## Create relational table that denotes where each name should map to 
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

wcvp_relations <- wcvp_relations %>% 
  mutate(word_count = str_count(acceptedName, "\\s+") + 1)

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

accepted_matched_names_wcvp <- filter(rematched_names_wcvp, taxon_status == "Accepted")
accepted_half_merge <- merge(wcvp_relations, accepted_matched_names_wcvp, by.x = "accepted_powo_id", by.y = "powo_id", all.x = TRUE)
synonym_matched_names_wcvp <- filter(rematched_names_wcvp, taxon_status != "Accepted")
synonym_half_merge <- merge(wcvp_relations, synonym_matched_names_wcvp, by.x = "synonym_powo_id", by.y = "powo_id", all.y = TRUE) 
merged_df <- rbind(accepted_half_merge, synonym_half_merge)
merged_df_less <- select(merged_df, acceptedName, synonym)

# clean up naming on cols
merged_df_t <- merged_df %>% 
  rename(acceptedNameMultMapsPossible = multipleMappingsPossible,
         acceptedNameSpacelessWCVPAuthorship = spacelessWCVPAuthorship, 
         acceptedNameAuthorshipMatch = authorshipMatch,
         acceptedNameMultMapResolutionPossible =  multMapResolutionPossible, 
         synonymNameWCVPAuthorship = synonym_taxon_authors)

merged_df_f <- merged_df_t %>% 
  filter(nameMatch %in% rematched_names_wcvp$nameMatch) %>% 
  distinct(.keep_all = TRUE)

# Merged_df denotes all relations and will be used to build nameAligned 
rematched_names_wcvp_synonyms <- filter(rematched_names_wcvp, taxon_status == "Synonym")
rematched_names_wcvp_synonyms <- rename(rematched_names_wcvp_synonyms, synonym = "nameMatch")
rematched_names_wcvp_synonyms_aligned <- merge(rematched_names_wcvp_synonyms, merged_df_less, by = "synonym", all.x = TRUE)
rematched_names_wcvp_synonyms_aligned <- distinct(rematched_names_wcvp_synonyms_aligned, .keep_all = TRUE)
rematched_names_wcvp_synonyms_aligned <- rename(rematched_names_wcvp_synonyms_aligned,
                                                nameMatch = synonym, nameAligned = acceptedName)

# To do -> fix these multiple mapping names, but for now throw them out since they make up few data points
rematched_names_wcvp_synonyms_aligned <- rematched_names_wcvp_synonyms_aligned %>% 
  group_by(nameMatch) %>% 
  mutate(n = n_distinct(nameAligned)) %>% 
  mutate(synAuthorCheck = case_when(
    n > 1 & tolower(gsub(" ", "", authorship)) == tolower(gsub(" ", "", wcvp_taxon_authors)) ~ TRUE, 
    n == 1 ~ TRUE,
    n > 1 & tolower(gsub(" ", "", authorship)) == tolower(gsub(" ", "", wcvp_taxon_authors)) ~ FALSE,
  )) %>% 
  filter(synAuthorCheck != FALSE)

#### Check if rematched names aligned acceptedName can be found in ncbi
ncbi_names3 <- ncbi_names2[!grepl(" cf\\. ", ncbi_names2$specificEpithet) & !grepl("\\d", ncbi_names2$specificEpithet), ]
rematched_names_with_ncbi <- merge(rematched_names_wcvp_synonyms_aligned, ncbi_names3, by.x = "nameAligned", by.y = "canonicalfull", all.x = TRUE)
rematched_names_with_ncbi <- rematched_names_with_ncbi %>% 
  distinct(.keep_all = TRUE)


rematched_names_with_ncbi %>% 
  filter(is.na(acceptedNameUsageID)) %>%  
  nrow()


test <- rematched_names_with_ncbi %>% 
  filter(!is.na(acceptedNameUsageID)) %>%  
  group_by(nameAligned) %>% 
  mutate(n_dups = n_distinct(acceptedNameUsageID))

# remove things that are mult mapping here, as there isnt really a way to deal with them

rem_syn_ncbi <- test %>% 
  filter(!n_dups > 1)

rem_names <- rem_syn_ncbi %>% 
  ungroup() %>% 
  select(nameMatch, concatenatedName, TropicosNameID, original_authorship = authorship, scientificName, ncbiTaxonRank = taxonRank, 
         ncbiTaxonomicStatus = taxonomicStatus, ncbiAcceptedNameUsageID = acceptedNameUsageID, 
         ncbiAcceptedNameUsageIDLess = acceptedNameUsageIDLess, ncbiKingdom = kingdom, ncbiPhylum = phylum,
         ncbiClass = class, ncbiOrder = order, ncbiFamily = family.y, ncbiGenus = genus.y, ncbiMappedSpecies = specificEpithet.y) %>% 
  mutate(multipleMappingsPossible = FALSE)

all_matched_names <- rbind(matched_names, rem_names)

# basically just count things up from here, then add up the discrepencies. 

n_not_aligned <- all_matched_names %>% 
  filter(is.na(ncbiAcceptedNameUsageID)) %>%
  distinct(nameMatch, .keep_all = TRUE) %>% 
  nrow()

n_aligned <- all_matched_names %>% 
  filter(!is.na(ncbiAcceptedNameUsageID)) %>% 
  distinct(ncbiMappedSpecies, .keep_all = TRUE) %>%  
  nrow()

(n_not_aligned/(n_not_aligned + n_aligned)) * 100 # ~ 36% of the names cant be aligned with ncbi 


# create a named vector to map IDs to names
id_to_name <- setNames(all_matched_names$ncbiMappedSpecies, all_matched_names$ncbiAcceptedNameUsageIDLess)

# replace the IDs in tip.label with corresponding names 
poa_tree <- full_tree
poa_tree$tip.label <- id_to_name[poa_tree$tip.label]

plot(poa_tree)

plot.phylo(poa_tree, cex=0.7, label.offset = 0.0005, no.margin=TRUE)
nodelabels(cex=0.7)


saveRDS(all_matched_names, file = "/home/millerjared/blue_guralnick/millerjared/BoCP/outputs/ncbi_matched_names_temp.rds")

saveRDS()



### Denote whether multiple mappings are possible for Synonyms  
multiple_s_mappings <- rematched_names_wcvp_synonyms_aligned %>% 
  group_by(nameMatch) %>% 
  mutate(numDiffAcceptedNames = length(unique((nameAligned)))) %>% 
  filter(numDiffAcceptedNames > 1)
syn_mult_maps <- filter(rematched_names_wcvp_synonyms_aligned, nameMatch %in% multiple_s_mappings$nameMatch)
rematched_names_wcvp_synonyms_aligned$synonymMultMapPossible <- ifelse(rematched_names_wcvp_synonyms_aligned$nameMatch %in% syn_mult_maps$nameMatch, TRUE, FALSE)

## create alignment info 
# aligned_rematched_wcvp <- rematched_names_wcvp %>% 
#   left_join(wcvp_relations_long, by = "nameMatch") %>% 
#   mutate(nameAligned = ifelse(relationType == "synonym", acceptedName, nameMatch))

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

### Coalesce the names via priorty col > itis > gbif
coal_df <- alternative_catalog_matches %>% 
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

## Designate whether the name match is a synonym or acceptedName
coal_df <- coal_df %>% mutate(taxon_status = ifelse(nameMatch == nameAligned, "Accepted", "Synonym"))

##

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
