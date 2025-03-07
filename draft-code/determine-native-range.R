### A prototype script for determining native ranges via the wcvp designations. 
# Author: JT Miller 
# Date: 08-08-2024
# Project: BoCP 

# load libraries 
library(data.table)
library(dplyr)
library(sf)
library(ggplot2)
library(lwgeom)

sf_use_s2(FALSE)
# load data 
wcvp_distribution_tb <- fread("./data/raw/wcvp_distribution.csv")
wcvp_name_tb <- fread("./data/raw/wcvp_names.csv")
na_plants_alignment <- fread("./data/processed/finalized_name_alignment_wcvp.csv")
level3_regions <- read_sf("./data/raw/level3-wgsrpd/level3.shp")
# select for relevant info 
wcvp_names <- wcvp_name_tb %>% 
  select(plant_name_id, taxon_name, taxon_authors, taxon_status)

# merge names with distribution info.
wcvp_names_w_dist <- merge(wcvp_names, wcvp_distribution_tb, by = "plant_name_id")

na_wcvp_names_w_dist <- wcvp_names_w_dist %>% 
  filter(taxon_name %in% na_plants_alignment$acceptedName)

na_wcvp_names_w_dist <- select(na_wcvp_names_w_dist, plant_name_id, taxon_name, continent, region_code_l2, region, area_code_l3,
                               area, introduced, extinct, location_doubtful)

# Load in occurrence data, denote whether native vs non-native. 
## Build a vector to read in the names by 
accepted_name_v <- unique(na_plants_alignment$acceptedNameParent)
accepted_name_filestyle_v <- gsub(" ", "-", accepted_name_v)
distinct_cleans <- list.files("./data/processed/distinct-clean-data/")
distinct_cleans <- gsub(".csv", "", distinct_cleans)
accepted_name_filestyle_v <- accepted_name_filestyle_v[accepted_name_filestyle_v %in% distinct_cleans]
accepted_name_reg_style <- gsub("-", " ", accepted_name_filestyle_v)
accepted_name_v <- accepted_name_v[accepted_name_v %in% accepted_name_reg_style]
for(i in 1:length(accepted_name_v)){
if(file.exists(paste0("./data/processed/full-spatial-clean-data/", accepted_name_filestyle_v[i], ".csv"))){
  occur_data <- fread(paste0("./data/processed/full-spatial-clean-data/", accepted_name_filestyle_v[i], ".csv"))
} else{
  print(paste("No occur data for", accepted_name_v[[i]]))
} # end of else no data. 

## make the occurrence data spatial 
occur_data_sf <- st_as_sf(occur_data, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)
occur_data_sf <- st_transform(occur_data_sf, crs = st_crs(level3_regions))

## intersect the data by region. 
occur_data_regioned <- st_join(occur_data_sf, level3_regions, join = st_intersects) # note the warning of bbox potentially invalid is due to level3_regions having a lon of -180.00003, the effect is so minimal im electing to ignore this.

## Determine if the occurrence is within native range
native_status_for_name <- na_wcvp_names_w_dist[taxon_name == accepted_name_v[i]]
occur_data_statuses <- occur_data_regioned %>% 
  rename(area_code_l3 = LEVEL3_COD) %>% 
  left_join(native_status_for_name, by = "area_code_l3") %>% 
  mutate(introduced = ifelse(!is.na(introduced), introduced, 2)) %>%  # Note that novel occurrences will be noted as 2. 
  mutate(nativeRangeStatus = ifelse(introduced == 0, TRUE, FALSE)) %>% 
  select(-taxon_name)

## Make a summary table 
regioned_data_summary <- data.frame(
  species = accepted_name_v[i], 
  nativeRangeCount = nrow(subset(occur_data_statuses, occur_data_statuses$introduced == 0)),
  introducedRangeCount = nrow(subset(occur_data_statuses, occur_data_statuses$introduced == 1)),
  undocumentedRangeCount = nrow(subset(occur_data_statuses, occur_data_statuses$introduced == 2))
)

## Write out this data 
fwrite(occur_data_statuses, paste0("./data/processed/regioned-data/", accepted_name_filestyle_v[i], ".csv"))
fwrite(regioned_data_summary, "./data/processed/regioned_summaries/region-summary-tb.csv" , append = TRUE, col.names = !file.exists("./data/processed/regioned_summaries/region-summary-tb.csv"))
## Make a nice map featuring the distribution of the plant 
occur_data_bbox <- st_bbox(occur_data_statuses)
true_native_regions <- native_status_for_name %>% 
  filter(introduced == 0)
true_introduced_regions <- native_status_for_name %>% 
  filter(introduced == 1)

region_fills <- level3_regions %>% 
  mutate(regionStatus = case_when(
    LEVEL3_COD %in% true_native_regions$area_code_l3 ~ "native range", 
    LEVEL3_COD %in% true_introduced_regions$area_code_l3 ~ "introduced range",
    TRUE ~ "undocumented"
  ))
# Calculate the current bounding box and region fills bounding box
occur_data_bbox <- st_bbox(occur_data_statuses)
region_fills_bbox <- st_bbox(region_fills)

if(nrow(occur_data) > 1){
# Calculate the actual area of region_fills and the intersection with occur_data_bbox
area_region_fills <- st_area(st_union(region_fills))
area_intersection <- st_area(st_intersection(st_union(region_fills), st_as_sfc(occur_data_bbox)))

# Convert areas to numeric for proportion calculation
area_region_fills_numeric <- as.numeric(area_region_fills)
area_intersection_numeric <- as.numeric(area_intersection)

# Proportion of region_fills shown in the current occurrence data bounding box
proportion_shown <- area_intersection_numeric / area_region_fills_numeric

hull <- st_convex_hull(occur_data_statuses)
hull_bbox <- st_bbox(hull)
# Check if the proportion is less than 20%
if (proportion_shown < 0.20) {
  
  hull <- st_convex_hull(occur_data_statuses)
  hull_b <- st_buffer(hull, 10)
  hull_bbox <- st_bbox(hull_b)
  # Adjust the bounding box to cover 1/3 of the region_fills area, centered on original bbox
  # x_center <- (occur_data_bbox[1] + occur_data_bbox[3]) / 2
  # y_center <- (occur_data_bbox[2] + occur_data_bbox[4]) / 2
  # 
  # x_range <- region_fills_bbox[3] - region_fills_bbox[1]
  # y_range <- region_fills_bbox[4] - region_fills_bbox[2]
  # 
  # new_x_range <- x_range / 3
  # new_y_range <- y_range / 3
  # 
  # # Update the bbox to the new adjusted range, centered on the original bbox center
  # occur_data_bbox[1] <- x_center - (new_x_range / 2)
  # occur_data_bbox[3] <- x_center + (new_x_range / 2)
  # occur_data_bbox[2] <- y_center - (new_y_range / 2)
  # occur_data_bbox[4] <- y_center + (new_y_range / 2)
}
color_mapping <- c("native range" = "darkgreen", "introduced range" = "goldenrod", "undocumented" = "grey")
ggplot() + 
  geom_sf(region_fills, mapping = aes(fill = regionStatus)) +
  scale_fill_manual(values = color_mapping) +
  geom_sf(occur_data_statuses, mapping = aes()) + 
  coord_sf(xlim = c(hull_bbox[1], hull_bbox[3]), ylim = c(hull_bbox[2], hull_bbox[4])) +
  theme_bw() +
  ggtitle(paste("Occurrence Records for", accepted_name_v[i])) + 
  xlab("Longitude") + 
  ylab("Latitude") +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(fill=guide_legend(title="Region Status"))

ggsave(paste0("./outputs/cleaned-data-range-maps/", accepted_name_filestyle_v[i], ".jpg"), width = 10, height = 10)
} # end of if there is more than 1 occurrence. 
}