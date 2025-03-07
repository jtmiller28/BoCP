### Title: North American Taxa Filter
### Author: JT Miller
### Date: 02-07-2025

### Purpose: Take the name alignment and filter down to just North American (NA) taxa

## Load libraries
library(data.table)
library(tidyverse)
library(sf)

## Load in name alignment, wcvp geographic info, and the ids as we need to append these back to link them
name_alignment <- fread("/home/millerjared/blue_guralnick/millerjared/BoCP/outputs/wcvp-ncbi-alignment.csv")
wcvp_geo <- fread("/home/millerjared/blue_guralnick/millerjared/BoCP/data/raw/wcvp_distribution.csv")
wcvp_name_ids <- fread("/home/millerjared/blue_guralnick/millerjared/BoCP/data/raw/wcvp_names.csv")
bot_regions <- read_sf("./data/raw/level3-wgsrpd/level3.shp")

## Figure out what botanical level 3 regions are in North America 
bot_regions_test <- bot_regions[361:369, ]
ggplot() + 
  geom_sf(bot_regions, mapping = aes()) + 
  geom_sf(bot_regions_test, mapping = aes(fill = "red")) +
  geom_sf_text(data = bot_regions_test, aes(label = LEVEL3_COD), 
               size = 2, color = "black")
  

## Keep the following string for filtering
na_string <- c("ALA", "ABT", "ASK", "ARI", "ARK", "BRC", "CAL", "COL", "CNT", "DEL",
               "GEO", "FLA", "IDA", "IOW", "ILL","INI", "KAN", "MAN", "LOU", "KTY", "MAI",
                "MNT", "MIN", "MIC", "MAS", "MSO", "MSI", "MRY", "MXC", "MXE", "MXG", "MXN",
               "MXS", "MXT", "NDA", "NCA", "NBR", "NEV", "NEB", "NFL", "NUN", "NSC", 
               "NWT", "NWM", "OHI", "NWJ", "NWH", "NWY", "ORE", "ONT", "OKL", "PEN", "PEI", "QUE",
               "RHO", "SAS", "SDA", "SCA", "TEX", "TEN", "UTA", "VRG", "VER", "WAS", "WIS", "WDC", "WVA",
               "YUK", "WYO", "LAB") 
ggplot() + 
  geom_sf(bot_regions, mapping = aes()) +
  geom_sf(filter(bot_regions, LEVEL3_COD %in% na_string), mapping = aes(fill = "red")) +
  ggtitle("North American Extent for BoCP Implementation")

na_bot_regions <- filter(bot_regions, LEVEL3_COD %in% na_string)
## filter ids to only accepted names, then retain only names that include NA
wcvp_name_ids <- wcvp_name_ids %>% 
  filter(taxon_status == "Accepted")
# affix wcvp_geo and wcvp_names togther
wcvp_info <- merge(wcvp_name_ids, wcvp_geo, by.x = "plant_name_id", by.y = "plant_name_id")
# organize by names that occur in NA
wcvp_info <- wcvp_info %>% 
  group_by(plant_name_id) %>% 
  mutate(inNA = ifelse(any(area_code_l3 %in% na_bot_regions$LEVEL3_COD), TRUE, FALSE))

wcvp_info_summary <- wcvp_info %>% 
  select(taxon_name, inNA) %>% 
  distinct()
name_alignment_wcvp <- name_alignment %>% 
  filter(source == "wcvp")
name_alignment_ncbi <- name_alignment %>% 
  filter(source == "ncbi") %>% 
  filter(nameStatus == "Accepted")
wcvp_alignment_summary <- name_alignment_wcvp %>% 
  select(alignedName, alignedParentName) %>% 
  distinct()


# now merge with our wcvp_ncbi_alignment
alignment_w_geo <- merge(wcvp_alignment_summary
                         ,wcvp_info_summary, by.x = "alignedName", by.y = "taxon_name", all.x = TRUE)
alignment_w_geo <- alignment_w_geo %>% 
  group_by(alignedParentName) %>% 
  mutate(parentInNA = ifelse(any(inNA == TRUE), TRUE, FALSE)) %>% 
  ungroup() %>% 
  filter(alignedParentName != "")

na_name_alignment <- alignment_w_geo %>% 
  filter(parentInNA == TRUE)

# check new names to see if these make sense. 
old_alignment <- fread("./data/processed/finalized_name_alignment_wcvp.csv")

old_alignment <- distinct(old_alignment, acceptedNameParent)

difference <- setdiff(na_name_alignment$alignedParentName, old_alignment$acceptedNameParent)
difference2 <- setdiff(old_alignment$acceptedNameParent, na_name_alignment$alignedParentName)
