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




