# Generate a list of all plant species in FL for Zoe's work
library(data.table)
library(sf)
library(tidyverse)

# First, for the names that arent derived from NCBI, we can use the WCVP geo database to figure out whats in FL
name_alignment <- fread("/blue/guralnick/millerjared/PlantSweepeR/data/processed/wcvp-ncbi-alignment-na.csv")
wcvp_geo <- fread("/blue/guralnick/millerjared/PlantSweepeR/data/raw/wcvp_distribution_2025_update.csv")
wcvp_name_ids <- fread("/blue/guralnick/millerjared/PlantSweepeR/data/raw/wcvp_names_2025_update.csv") # the backbone because we need to link up ids to names

# from wcvp-north-america-delim
na_string <- c("ALA", "ABT", "ASK", "ARI", "ARK", "BRC", "CAL", "COL", "CNT", "DEL",
               "GEO", "FLA", "IDA", "IOW", "ILL","INI", "KAN", "MAN", "LOU", "KTY", "MAI",
               "MNT", "MIN", "MIC", "MAS", "MSO", "MSI", "MRY", "MXC", "MXE", "MXG", "MXN",
               "MXS", "MXT", "NDA", "NCA", "NBR", "NEV", "NEB", "NFL", "NUN", "NSC", 
               "NWT", "NWM", "OHI", "NWJ", "NWH", "NWY", "ORE", "ONT", "OKL", "PEN", "PEI", "QUE",
               "RHO", "SAS", "SDA", "SCA", "TEX", "TEN", "UTA", "VRG", "VER", "WAS", "WIS", "WDC", "WVA",
               "YUK", "WYO", "LAB")
# grab plants that are in NA
wcvp_geo_NA <- wcvp_geo %>% filter(grepl(paste(na_string, collapse = "|"), area_code_l3))
name_alignment_wcvp <- name_alignment %>% filter(nameStatus == "Accepted" & source == "wcvp")

# filter to just accepted names
wcvp_name_ids <- wcvp_name_ids %>% 
  filter(taxon_status == "Accepted")

# affix wcvp_geo and wcvp_names togther
wcvp_info <- merge(wcvp_name_ids, wcvp_geo_NA, by.x = "plant_name_id", by.y = "plant_name_id")

# combine this table with name aligment for wcvp
NA_names_geo <- merge(name_alignment_wcvp, wcvp_info, by.x = "alignedName", by.y = "taxon_name") # use aligned name, if any fall into parent then it'll by default be called 

# summarize parent name with aligned name, where if any of the aligned names that fall within a parent are in NA, it is a NA species
NA_names_parent_geo <- NA_names_geo %>% 
  group_by(alignedParentName) %>% 
  mutate(NA_species = ifelse(any(area_code_l3 %in% na_string), TRUE, FALSE)) %>% 
  mutate(NA_native = ifelse(any(introduced == 0), TRUE, FALSE)) %>% 
  mutate(NA_introduced = ifelse(any(introduced == 1), TRUE, FALSE)) %>% 
  ungroup() %>% 
  select(alignedName, alignedParentName, area, NA_species, NA_native, NA_introduced) %>% 
  distinct()

mult_map_check <- NA_names_parent_geo %>% 
  group_by(alignedName) %>% 
  summarize(n = length(unique(alignedParentName))) # no mult maps, proceed as normal

fwrite(NA_names_parent_geo, "/blue/guralnick/millerjared/BoCP/data/processed/NA-names-2025-update.csv")
