wcvp_taxonomy_asparagales 

wcvp_geo <- fread("/blue/guralnick/millerjared/PlantSweepeR/data/raw/wcvp_distribution_2025_update.csv")
wcvp_name_ids <- fread("/blue/guralnick/millerjared/PlantSweepeR/data/raw/wcvp_names_2025_update.csv")
bot_regions <- read_sf("/blue/guralnick/millerjared/PlantSweepeR/data/raw/level3-wgsrpd/level3.shp")
wcvp_backbone <- fread("/blue/guralnick/millerjared/PlantSweepeR/data/raw/wcvp_names_2025_update.csv")
## Botanical Regions strings that are Yukatan through Canada
## Keep the following string for filtering
na_string <- c("ALA", "ABT", "ASK", "ARI", "ARK", "BRC", "CAL", "COL", "CNT", "DEL",
               "GEO", "FLA", "IDA", "IOW", "ILL","INI", "KAN", "MAN", "LOU", "KTY", "MAI",
               "MNT", "MIN", "MIC", "MAS", "MSO", "MSI", "MRY", "MXC", "MXE", "MXG", "MXN",
               "MXS", "MXT", "NDA", "NCA", "NBR", "NEV", "NEB", "NFL", "NUN", "NSC", 
               "NWT", "NWM", "OHI", "NWJ", "NWH", "NWY", "ORE", "ONT", "OKL", "PEN", "PEI", "QUE",
               "RHO", "SAS", "SDA", "SCA", "TEX", "TEN", "UTA", "VRG", "VER", "WAS", "WIS", "WDC", "WVA",
               "YUK", "WYO", "LAB")
# plotted for clarity
ggplot() + 
  geom_sf(bot_regions, mapping = aes()) +
  geom_sf(filter(bot_regions, LEVEL3_COD %in% na_string), mapping = aes(fill = "red")) +
  ggtitle("North American Extent for BoCP Implementation")

# apply filter to shapefile regions to only include shapes from NA 
na_bot_regions <- filter(bot_regions, LEVEL3_COD %in% na_string)

## filter ids to only accepted names, then retain only names that include NA
wcvp_name_ids <- wcvp_name_ids %>% 
  filter(taxon_status == "Accepted")

# affix wcvp_geo and wcvp_names togther
wcvp_info <- merge(wcvp_name_ids, wcvp_geo, by.x = "plant_name_id", by.y = "plant_name_id")
wcvp_info <- wcvp_info %>% select(plant_name_id, taxon_name, area_code_l3, introduced, extinct, location_doubtful)

wcvp_asp_info <- wcvp_info %>% 
  filter(taxon_name %in% wcvp_taxonomy_asparagales$alignedName) %>% 
  filter(introduced == 0) %>%  # means native
  filter(area_code_l3 %in% na_string) %>% 
  distinct(taxon_name)

# if these are in the taxonomy, use that

wcvp_taxonomy_asparagales_natives <- wcvp_taxonomy_asparagales %>% 
  filter(name %in% wcvp_asp_info$taxon_name)

length(unique(wcvp_taxonomy_asparagales_natives$alignedParentName))

wcvp_taxonomy_asparagales_natives <- wcvp_taxonomy_asparagales_natives %>% 
  select(alignedParentName, wcvpAlignedFamily) %>% 
  rename(family = wcvpAlignedFamily) %>% 
  left_join(asparagales, by = "family")
