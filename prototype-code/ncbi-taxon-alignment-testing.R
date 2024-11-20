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
ncbi_names2 <- merge(ncbi_names, ncbi_parsed_names_sl, by.x = "scientificName", by.y = "verbatim")
## Do a left-join with our parent list and NCBI's taxonomy 
matched_names <- merge(bocp_parsed_names_df, ncbi_names2, by.x = "canonicalfull", by.y = "canonicalfull", all.x = TRUE)
## Remove Metazoa kindgom as those are not plants
matched_names <- filter(matched_names, kingdom == "Viridiplantae" | is.na(kingdom)) # retain NAs for now as Im not sure what these are...
## Select and Rename for standardization
matched_names <- matched_names %>% 
  select(canonicalfull, concatenatedName, TropicosNameID, authorship, scientificName, ncbiTaxonRank = taxonRank, 
         ncbiTaxonomicStatus = taxonomicStatus, ncbiAcceptedNameUsageID = acceptedNameUsageID, 
         ncbiAcceptedNameUsageIDLess = acceptedNameUsageIDLess, ncbiKingdom = kingdom, ncbiPhylum = phylum,
         ncbiClass = class, ncbiOrder = order, ncbiFamily = family, ncbiGenus = genus, ncbiMappedSpecies = specificEpithet)
################################################################################################################################################################################################################################################################

# Summaries + Multiple Mapping situation ############################################################################################################################################################################################################################
test <- matched_names %>%
  group_by(canonicalfull) %>%
  mutate(num_different_paths = n_distinct(ncbiAcceptedNameUsageIDLess))

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
