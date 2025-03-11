### A script for determining native ranges via the wcvp designations. 
# Author: JT Miller 
# Date: 09-25-2024
# Project: BoCP 

# Load Libraries
library(data.table)
library(dplyr)
library(sf)
library(ggplot2)

# set dir 
setwd("/home/millerjared/blue_guralnick/millerjared/BoCP/")

# Set up array logic
start_num <- as.numeric(Sys.getenv("START_NUM"))
task_id <- as.numeric(start_num)
part <- paste0("part", task_id)

# Must be turned off for map subsetting to operate correctly 
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
  filter(taxon_name %in% na_plants_alignment$acceptedName) # use acceptedName as this will create the full range for the parentTaxon

# select for relevant info
na_wcvp_names_w_dist <- select(na_wcvp_names_w_dist, plant_name_id, taxon_name, continent, region_code_l2, region, area_code_l3,
                               area, introduced, extinct, location_doubtful)

# determine parentTaxon using the name alignment's criteria 
na_wcvp_parentTaxon_w_dist <- na_wcvp_names_w_dist %>% 
  left_join(na_plants_alignment, by = join_by(taxon_name == acceptedName)) %>% 
  rename(parentTaxon = acceptedNameParent) %>% 
  select(names(na_wcvp_names_w_dist), parentTaxon)

# build a vector of names to index by, then index according to array tasks
accepted_name_v <- unique(na_plants_alignment$acceptedNameParent)

# index the vector for what names we've already finished 
names_done <- list.files("./data/processed/regioned-data/all/")
names_done <- gsub(".csv", "", names_done)
names_done <- gsub("-", " ", names_done)
unfinished_accepted_names <- setdiff(accepted_name_v, names_done)
accepted_name <- unfinished_accepted_names[task_id] # index by task ID
print(paste("working on accepted name:", accepted_name))
accepted_name_filestyle <- gsub(" ", "-", accepted_name)

# Function to write errors to a internal script error
log_error <- function(error_message) {
  log_file <- "./logs/internal-native-range-error/"
  write(paste(Sys.time(), " - ERROR: ", error_message, sep=""), file = log_file, append = TRUE)
}
tryCatch({
  if(file.exists(paste0("./data/processed/coord-clean-data/", accepted_name_filestyle, ".csv"))){
    
    # load the data
    occur_data <- fread(paste0("./data/processed/coord-clean-data/", accepted_name_filestyle, ".csv"))
    ## make the occurrence data spatial 
    occur_data_sf <- st_as_sf(occur_data, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) # Note we could use the roundedLon & roundedLat values, but we're buffering so much this shouldnt matter.
    ## buffer the occurrences to account for some wiggleroom on native/introduced/novel statuses
    crs_aea <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m" # albers equal area centered on NA
    occur_data_reproj <- st_transform(occur_data_sf, crs = crs_aea) # reproj to Mercator proj so that we're in a projected crs for buffering
    occur_data_b <- st_buffer(occur_data_reproj, dist = 10000)
    occur_data_b <- st_transform(occur_data_b, crs = 4326)
    
    # determine the occurrence's native/introdced/novel statuses 
    statuses_for_name <- na_wcvp_parentTaxon_w_dist[parentTaxon == accepted_name, ]
    ## deal with collasping the subspecific names into the parentTaxon concerning distribution
    statuses_for_name <- statuses_for_name %>% 
      distinct(across(-c(plant_name_id, taxon_name)), .keep_all = TRUE) %>% # Take distinct values for parentTaxon rather than the taxon
      group_by(parentTaxon, area_code_l3) %>% # group up values of the parentTaxon and the location
      summarize(across(everything(), ~ first(.)), # Priortize natives vs introduced values per grouping  
                introduced = min(introduced)) %>% # take the 0 value instead of 1 if it exists
      ungroup()
    
    ## join occur data with the regions
    occur_data_regioned <- st_join(occur_data_b, level3_regions, join = st_intersects) # note the warning of bbox potentially invalid is due to level3_regions having a lon of -180.00003, the effect is so minimal im electing to ignore this.
    
    ## and join statuses with occur data, note that this deals with the duplications of records created by the buffering
    occur_data_statuses <- occur_data_regioned %>% 
      rename(area_code_l3 = LEVEL3_COD) %>% 
      left_join(statuses_for_name, by = "area_code_l3") %>% 
      mutate(introduced = ifelse(!is.na(introduced), introduced, 2)) %>%  # Note that novel occurrences will be noted as 2. 
      group_by(ID) %>% # locate duplicates where present
      summarize(across(everything(), ~ first(.)), # find the lowest value to take in priority 
                introduced = min(introduced)) %>% # we're interested if it becomes native, so take the min (0) being native, (1) being introduced, (2) being novel    ungroup() %>% 
      mutate(nativeRangeStatus = ifelse(introduced == 0, TRUE, FALSE))
    
    ## Make a summary table 
    regioned_data_summary <- data.frame(
      parentTaxon = accepted_name, 
      nativeRangeCount = nrow(subset(occur_data_statuses, occur_data_statuses$introduced == 0)),
      introducedRangeCount = nrow(subset(occur_data_statuses, occur_data_statuses$introduced == 1)),
      undocumentedRangeCount = nrow(subset(occur_data_statuses, occur_data_statuses$introduced == 2))
    )
    
    
    # Produce mapping for visual inspection
    
    ## determine the natives vs introduced for color fill choice
    true_native_regions <- statuses_for_name %>% 
      filter(introduced == 0)
    true_introduced_regions <- statuses_for_name %>% 
      filter(introduced == 1)
    
    region_fills <- level3_regions %>% 
      mutate(regionStatus = case_when(
        LEVEL3_COD %in% true_native_regions$area_code_l3 ~ "native range", 
        LEVEL3_COD %in% true_introduced_regions$area_code_l3 ~ "introduced range",
        TRUE ~ "undocumented"
      ))
    
    ## determine the zoom ratio for mapping purposes
    ### calculate the current bounding box and region fills bounding box
    occur_data_bbox <- st_bbox(occur_data_statuses)
    region_fills_bbox <- st_bbox(region_fills)
    
    if(nrow(occur_data) > 1){
      ### calculate the actual area of region_fills and the intersection with occur_data_bbox
      area_region_fills <- st_area(st_union(region_fills))
      area_intersection <- st_area(st_intersection(st_union(region_fills), st_as_sfc(occur_data_bbox)))
      
      ### convert areas to numeric for proportion calculation
      area_region_fills_numeric <- as.numeric(area_region_fills)
      area_intersection_numeric <- as.numeric(area_intersection)
      
      ### proportion of region_fills shown in the current occurrence data bounding box
      proportion_shown <- area_intersection_numeric / area_region_fills_numeric
      
      hull <- st_convex_hull(occur_data_statuses)
      hull_bbox <- st_bbox(hull)
      ### check if the proportion is less than 20%
      if (proportion_shown < 0.20) {
        
        hull <- st_convex_hull(occur_data_statuses)
        hull_b <- st_buffer(hull, 10)
        hull_bbox <- st_bbox(hull_b)
      }
      
      ## and map
      color_mapping <- c("native range" = "darkgreen", "introduced range" = "goldenrod", "undocumented" = "grey")
      ggplot() + 
        geom_sf(region_fills, mapping = aes(fill = regionStatus)) +
        scale_fill_manual(values = color_mapping) +
        geom_sf(occur_data_statuses, mapping = aes(), fill = "darkorange", alpha = 0.9) + 
        coord_sf(xlim = c(hull_bbox[1], hull_bbox[3]), ylim = c(hull_bbox[2], hull_bbox[4])) +
        theme_bw() +
        ggtitle(paste("Occurrence Records for", accepted_name)) + 
        xlab("Longitude") + 
        ylab("Latitude") +
        theme(plot.title = element_text(hjust = 0.5)) +
        guides(fill=guide_legend(title="Region Status"))
      
      ggsave(paste0("./outputs/occur-data-range-status-maps/", accepted_name_filestyle, ".jpg"), width = 10, height = 10)
    } # end of if there is more than one occur record make the map...
    
    ## remove the geom col before writing out these data
    occur_data_statuses <- occur_data_statuses %>% 
      st_drop_geometry()
    
    # Write out this data 
    fwrite(subset(occur_data_statuses, occur_data_statuses$introduced == 0), paste0("./data/processed/regioned-data/native/", accepted_name_filestyle, ".csv"))
    fwrite(subset(occur_data_statuses, occur_data_statuses$introduced == 1), paste0("./data/processed/regioned-data/introduced/", accepted_name_filestyle, ".csv"))
    fwrite(occur_data_statuses, paste0("./data/processed/regioned-data/all/", accepted_name_filestyle, ".csv"))
    fwrite(regioned_data_summary, "./data/processed/cleaning-summaries/native-range-summary-tb.csv" , append = TRUE, col.names = !file.exists("./data/processed/cleaning-summaries/native-range-summary-tb.csv"))
    
    
  } else{ # If the name does not exist at this step, include this in the final output table. 
    print(paste("No occur data for", accepted_name))
    regioned_data_summary <- data.frame(
      parentTaxon = accepted_name, 
      nativeRangeCount = NA,
      introducedRangeCount = NA,
      undocumentedRangeCount = NA
    )
    fwrite(regioned_data_summary, "./data/processed/cleaning-summaries/native-range-summary-tb.csv" , append = TRUE, col.names = !file.exists("./data/processed/cleaning-summaries/native-range-summary-tb.csv"))
  }
  
}, error = function(e) {
  # Log any errors
  log_error(conditionMessage(e))
})
print(part)
