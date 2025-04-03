# 005 Spatially flagging and rounding the dataset
# Author: JT Miller 
# Date: 03-12-2024
# Project: BoCP 

## Purpose: Taking the dataset we want to 1) where decimal degrees data is provided, round coordinate percision, 2) where datum is provided, flag non-WGS84 datums, 3) where missing coordinates or impossible coordinates are provided flag

# load libraries 
library(data.table)
library(tidyverse)
library(sf)

# Set up array logic
start_num <- as.numeric(Sys.getenv("START_NUM"))
task_id <- as.numeric(start_num)
part <- paste0("part", task_id)

# read the taxonomic harmonized list of names + synonyms
name_list <- readRDS("/blue/guralnick/millerjared/BoCP/data/processed/name_list.rds")

# check for names that have already been retrieved
list_found_names <- list.files("/blue/guralnick/millerjared/BoCP/data/processed/flagged-species-occs/")
list_found_names <- gsub("-", " ", list_found_names)
list_found_names <- gsub(".csv", "", list_found_names)
# retrieve accepted names 
accepted_name_retrieval_list <- list()
for(i in 1:length(name_list)){
  accepted_names_retrieval <- name_list[[i]][1]
  accepted_name_retrieval_list[[i]] <- accepted_names_retrieval
}
accepted_name_retrieval_vector <- do.call(rbind, accepted_name_retrieval_list)
# check diff between dir names and total names to call
names_not_found <- setdiff(accepted_name_retrieval_vector, list_found_names)
# extract 
# Extract elements where any name in the list matches a target name
name_list_update <- name_list[sapply(name_list, function(x) any(x %in% names_not_found))]
# using taskid as our index, extract names to read 
names_to_retrieve <- name_list_update[[task_id]]
accepted_name <- names_to_retrieve[1]
accepted_name_file_style <- gsub(" ", "-", accepted_name)
print(paste("reading in data for", accepted_name))
# read in the species occurrence table 
occ_data <- fread(paste0("/blue/guralnick/millerjared/BoCP/data/processed/species-occs/", accepted_name_file_style, ".csv"))
print(paste("Standardizing coord fields for", accepted_name))
# 1) Where decimal degrees geocoordinates are provided, create a new field called roundedLongitude & roundedLatitude, these are to avoid overpercision when modeling. 
occ_data <- occ_data %>% 
  mutate(decimalLatitude = as.character(decimalLatitude), # due to input errors, we'll first change the class to character so that commas may be removed. 
         decimalLongitude = as.character(decimalLongitude)) %>% 
  mutate(decimalLatitude = gsub(",", ".", decimalLatitude),
         decimalLongitude = gsub(",", ".", decimalLongitude)) %>% 
  mutate(decimalLatitude = as.numeric(decimalLatitude),
         decimalLongitude = as.numeric(decimalLongitude)) %>% 
  mutate(roundedLatitude = ifelse(!is.na(decimalLatitude) & decimalLatitude != 0.00, round(decimalLatitude, 2), NA),
         roundedLongitude = ifelse(!is.na(decimalLongitude) & decimalLongitude != 0.00, round(decimalLongitude, 2), NA))

# 2) Where datum is provided, flag non-WGS84 data
occ_data <- occ_data %>% 
  mutate(wgs84Datum = ifelse(is.na(geodeticDatum) | 
                               toupper(geodeticDatum) == "WGS84" | 
                               toupper(geodeticDatum) == "WORLD GEODETIC SYSTEM 1984" |
                               geodeticDatum == "" | 
                               toupper(geodeticDatum) == "NO DISPONIBLE" | 
                               toupper(geodeticDatum) == "EPSG:4326", TRUE, FALSE))
print(paste("Flagging coordinate Issues for", accepted_name))
# 3) flag coordinates for either not being provided OR being impossible
occ_data <- occ_data %>% 
  mutate(coordinateIssue = ifelse(is.na(decimalLatitude) | 
                                    is.na(decimalLongitude) |
                                    decimalLatitude == 0.00 |
                                    decimalLongitude == 0.00 | 
                                    decimalLatitude > 90 |
                                    decimalLatitude < -90 |
                                    decimalLongitude > 180 | 
                                    decimalLongitude < -180 |
                                    decimalLatitude == "" | 
                                    decimalLongitude == "", TRUE, FALSE)) %>% 
  mutate(coordUncertaintyHigh = ifelse(!is.na(coordinateUncertaintyInMeters) & coordinateUncertaintyInMeters > 1000, TRUE, FALSE))
print(paste("flagging true coordinates witheld for", accepted_name))
# 4) flag data for locality information withheld (we're referencing Natalie's gatoRs implementation of remove_skewed()  indicating "Coordinate uncertainty increased")
occ_data <- occ_data %>% 
  mutate(informationWithheld = as.character(informationWithheld)) %>% # when there is nothing in this field for all records in the csv it gets converted to a logical output
  mutate(informationWithheld = ifelse(is.na(informationWithheld), "", informationWithheld)) %>% # NAs will cause issues if unaccounted for in the proceeding line of code. 
  mutate(trueCoordsWithheld = ifelse(
                                     toupper(informationWithheld) == "SPECIFIC LOCALITY INFORMATION MAY BE WITHHELD" |
                                     toupper(informationWithheld) == "PRECISE LOCATION INFORMATION AVAILABLE ON REQUEST" |
                                     toupper(informationWithheld) == "LOCATION INFORMATION NOT GIVEN FOR ENDANGERED SPECIES" |
                                     toupper(informationWithheld) == "SENSITIVE LOCATION DATA WITHHELD" |
                                     toupper(informationWithheld) == "LA UBICACIÓN NO ES PROVISTA PARA ESPECIES AMENAZADAS" |
                                     toupper(informationWithheld) == "LOCALITY AND PRECISE COORDINATES WITHHELD TO PROTECT SENSITIVE DATA." |
                                     toupper(informationWithheld) == "LOCALITY, COLLECTORS, FULL EVENT DATE, AND OTHER FIELDS WITHHELD TO PROTECT SENSITIVE INFORMATION." |
                                     toupper(informationWithheld) == "LOCATION DATA AVAILABLE TO QUALIFIED RESEARCHERS ON REQUEST." |
                                     toupper(informationWithheld) == "PRECISE LOCALITY INFORMATION NOT PROVIDED FOR ENDANGERED SPECIES." |
                                     toupper(informationWithheld) == "LOCALITY DATA HIDDEN DUE TO SENSITITIVE STATUS OF TAXON" |
                                     toupper(informationWithheld) == "MASK PART ATTRIBUTE LOCATION" |
                                     toupper(informationWithheld) == "LOCATION INFORMATION NOT GIVEN PRECISELY : SENSITIVE TAXON" |
                                     toupper(informationWithheld) == "SE OCULTARON COORDENADAS. SE TRATA DE ESPECIE SENSIBLE. SI REQUIERE MÁS INFORMACIÓN CONTÁCTESE CON EL CREADOR DEL RECURSO." |
                                     toupper(informationWithheld) == "OCCURRENCE IS PLACED IN THE GRID SQUARE CENTROID; REAL COORDINATES, HABITAT DETAILS AND VOUCHER INFORMATION (IF PRESENT) ARE OBSCURED." |
                                     toupper(informationWithheld) == "ORIGINAL LOCATIONS AVAILABLE UPON REQUEST"  |
                                     toupper(informationWithheld) == "LOCATION, COORDINATES, INDIVIDUALS, LIFESTAGE, BEHAVIOUR ETC" |
                                     toupper(informationWithheld) == "BELOW COUNTY-LEVEL LOCALITY DATA" |
                                     toupper(informationWithheld) == "ASK FOR PRECISE COORDINATES" |
                                     toupper(informationWithheld) == "OCCURRENCE BLURRED" |
                                     toupper(informationWithheld) == "LOCALITY AND EXACT COORDINATES NOT GIVEN FOR ENDANGERED SPECIES." |
                                     toupper(informationWithheld) == "LAS COORDENADAS (VERBATIMCOORDINATES, VERBATIMLATITUDE, VERBATIMLONGITUDE), ALTITUD (VERBATIMELEVATION) Y EXPOSICIÓN (ASPECT) DE LAS PARCELAS DE MONITOREO Y SUS ESPECIES, HAN SIDO RESTRINGIDAS A SOLICITUD DE LOS ADMINISTRADORES DE ESTAS ÁREAS PROTEGIDAS EN FAVOR DE LA CONSERVACIÓN DE LAS MISMAS." |
                                     toupper(informationWithheld) == "LOCATION INFORMATION NOT GIVEN FOR SPECIES AT RISK"|
                                     toupper(informationWithheld) == "SUB-COUNTRY 3 AND BELOW SUPPRESSED"|
                                     grepl("DETAIL LOCATION DATA UNDISCLOSED", toupper(informationWithheld)) |
                                     grepl("ASK ABOUT DECIMALLATITUDE & DECIMALLONGITUDE", toupper(informationWithheld)) |
                                     grepl("OCCURRENCE IS PLACED IN THE GRID SQUARE CENTROID", toupper(informationWithheld)) |
                                     grepl("MASK COORDINATES", toupper(informationWithheld)) |
                                     grepl("MASK LOCALITY", toupper(informationWithheld)) |
                                     grepl("SOME LOCALITY DATA IS WITHHELD", toupper(informationWithheld)) |
                                     grepl("LOCALITY DETAILS WITHHELD", toupper(informationWithheld)) |
                                     grepl("COORDINATES HIDDEN", toupper(informationWithheld)) |
                                     grepl("PRECISE LOCATION INFORMATION WITHHELD", toupper(informationWithheld)) |
                                     grepl("COORDINATES HIDDEN", toupper(informationWithheld)) |
                                     grepl("LOCALITY NOT GIVEN", toupper(informationWithheld)) |
                                     grepl("GEOGRAPHIC INFORMATION GENERALIZED DURING AGGREGATION AT THE REQUEST OF THE PRODUCER", toupper(informationWithheld)) |
                                     grepl("FIELD VALUES REDACTED.*(DECIMALLATITUDE|DECIMALLONGITUDE)", toupper(informationWithheld)),
                                     TRUE, FALSE)) %>% 
  mutate(coordUncertaintyInc = ifelse(
                                     grepl("COORDINATE UNCERTAINTY INCREASED", toupper(informationWithheld)),
                                     TRUE, FALSE))

# Range Status #######################################################
print(paste("Preparing range status field for", accepted_name))
## using wcvp's info tables, assign native/introduced/undocumented status for spatial placement of occurrences.
wcvp_distribution_tb <- fread("/blue/guralnick/millerjared/BoCP/data/raw/wcvp_distribution.csv")
wcvp_name_tb <- fread("/blue/guralnick/millerjared/BoCP/data/raw/wcvp_names.csv")
na_plants_alignment <- fread("/blue/guralnick/millerjared/BoCP/data/processed/wcvp-ncbi-alignment-na.csv")
level3_regions <- read_sf("/blue/guralnick/millerjared/BoCP/data/raw/level3-wgsrpd/level3.shp")
## Must be turned off for map subsetting to operate correctly 
sf_use_s2(FALSE)
## filter down alignment to capture all the alignedNames under our alignedParentName. Conditionally disclude accepted names that are built only by ncbi 
ncbi_accepted_names <- na_plants_alignment %>% filter(source == "ncbi" & nameStatus == "Accepted")
if(!accepted_name %in% ncbi_accepted_names$alignedParentName){
na_plants_alignment <- na_plants_alignment %>% 
  filter(alignedParentName == accepted_name)
## select for relevant info & filter to the relevant name info for this extraction
wcvp_names <- wcvp_name_tb %>% 
  select(plant_name_id, taxon_name, taxon_authors, taxon_status) %>% 
  filter(taxon_name %in% na_plants_alignment$alignedName) # grab all possible 
## merge to grab geographic info according to plant_name_id
wcvp_names_w_dist <- merge(wcvp_names, wcvp_distribution_tb, by = "plant_name_id")
## select for relevant info
wcvp_names_w_dist <- select(wcvp_names_w_dist, plant_name_id, taxon_name, continent, region_code_l2, region, area_code_l3,
                               area, introduced, extinct, location_doubtful)
## determine parentTaxon using the name alignment's criteria 
wcvp_parentTaxon_w_dist <- wcvp_names_w_dist %>% 
  left_join(na_plants_alignment, by = join_by(taxon_name == alignedName)) %>% 
  rename(parentTaxon = alignedParentName) %>% 
  select(names(wcvp_names_w_dist), parentTaxon) %>% 
  distinct()
# Spatial Subsetting for relevant Data
occ_data <- occ_data %>% 
    mutate(tempID = 1:n()) # use for end merge.
## remove records that do not have the attributes needed for looking at them spatially
occ_spatial_data <- occ_data %>% 
    filter(wgs84Datum == TRUE) %>%  # necessary for projecting correctly
    filter(coordinateIssue == FALSE) %>%  # needs to be clean and correct coords
    filter(trueCoordsWithheld == FALSE) # only want true clean coords 
## make the data spatial 
occ_data_sf <- st_as_sf(occ_spatial_data, coords = c("roundedLongitude", "roundedLatitude"), crs = 4326, remove = FALSE)
## reproj in equal area
crs_aea <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m" # albers equal area centered on NA
occ_data_reproj <- st_transform(occ_data_sf, crs = crs_aea) # reproj to Mercator proj so that we're in a projected crs for buffering
## buffer points to account for some wiggle room in regards to the ranges we're categorizing by
occ_data_b <- st_buffer(occ_data_reproj, dist = 10000) # allow 10km of wiggle room
occ_data_b <- st_transform(occ_data_b, crs = 4326) # reproj back so we can label
## Summarize the possible botanical region codes for a parentTaxon (the alignedName)
statuses_for_name <- wcvp_parentTaxon_w_dist %>% 
    distinct(across(-c(plant_name_id, taxon_name)), .keep_all = TRUE) %>% # Take distinct values for parentTaxon rather than the taxon
    group_by(parentTaxon, area_code_l3) %>% # group up values of the parentTaxon and the location
    summarize(across(everything(), ~ first(.)), # Priortize natives vs introduced values per grouping  
              introduced = min(introduced)) %>% # take the 0 value instead of 1 if it exists
    ungroup()
## join occ data with these regions. 
occ_data_regioned <- st_join(occ_data_b, level3_regions, join = st_intersects)
## and join statuses with occur data, note that this deals with the duplications of records created by the buffering
occ_data_statuses <- occ_data_regioned %>% 
    rename(area_code_l3 = LEVEL3_COD) %>% 
    left_join(statuses_for_name, by = "area_code_l3") %>% 
    mutate(introduced = ifelse(!is.na(introduced), introduced, 2)) %>%  # Note that novel occurrences will be noted as 2. 
    group_by(tempID) %>% # locate duplicates where present
    summarize(across(everything(), ~ first(.)), # find the lowest value to take in priority 
              introduced = min(introduced)) %>% # we're interested if it becomes native, so take the min (0) being native, (1) being introduced, (2) being novel    ungroup() %>% 
    mutate(wcvpRangeStatus = case_when(
      introduced == 0 ~ "native",
      introduced == 1 ~ "introduced", 
      introduced == 2 ~ "undocumented"
    ))
  
# Build Range Map Visual
print(paste("Building plot visuals for", accepted_name))
## create fills for native, introduced, and undocumented ranges
true_native_regions <- statuses_for_name %>% 
    filter(introduced == 0) # we're looking at all ranges for taxon
true_introduced_regions <- statuses_for_name %>% 
    filter(introduced == 1)
region_fills <- level3_regions %>% # create fills accordingly 
    mutate(regionStatus = case_when(
      LEVEL3_COD %in% true_native_regions$area_code_l3 ~ "native range", 
      LEVEL3_COD %in% true_introduced_regions$area_code_l3 ~ "introduced range",
      TRUE ~ "undocumented"
    ))
## create a plotting view, allowing for some conditional buffering based on the point occurrences extent
occ_data_bbox <- st_bbox(occ_data_statuses) ## determine the zoom ratio for mapping purposes
region_fills_bbox <- st_bbox(region_fills) ### calculate the current bounding box and region fills bounding box
occ_native_data_bbox <- st_bbox(filter(occ_data_statuses, wcvpRangeStatus == "native"))
  # conditional to adjust for plotting complexity...
  if(nrow(occ_data_sf) > 1){
    ### calculate the actual area of region_fills and the intersection with occur_data_bbox
    area_region_fills <- st_area(st_union(region_fills))
    area_intersection <- st_area(st_intersection(st_union(region_fills), st_as_sfc(occ_data_bbox)))
    
    ### convert areas to numeric for proportion calculation
    area_region_fills_numeric <- as.numeric(area_region_fills)
    area_intersection_numeric <- as.numeric(area_intersection)
    
    ### proportion of region_fills shown in the current occurrence data bounding box
    proportion_shown <- area_intersection_numeric / area_region_fills_numeric
    
    hull <- st_convex_hull(occ_data_statuses)
    hull_bbox <- st_bbox(hull)
    ### check if the proportion is less than 20%
    if (proportion_shown < 0.20) {
      
      hull <- st_convex_hull(occ_data_statuses)
      hull_b <- st_buffer(hull, 10)
      hull_bbox <- st_bbox(hull_b)
    }
    
    ## and map
    color_mapping <- c("native range" = "darkgreen", "introduced range" = "goldenrod", "undocumented" = "grey")
    ggplot() + 
      geom_sf(region_fills, mapping = aes(fill = regionStatus)) +
      scale_fill_manual(values = color_mapping) +
      geom_sf(occ_data_statuses, mapping = aes(), fill = "darkorange", alpha = 0.9) + 
      coord_sf(xlim = c(hull_bbox[1], hull_bbox[3]), ylim = c(hull_bbox[2], hull_bbox[4])) +
      theme_bw() +
      ggtitle(paste("Occurrence Records for", accepted_name)) + 
      xlab("Longitude") + 
      ylab("Latitude") +
      theme(plot.title = element_text(hjust = 0.5)) +
      guides(fill=guide_legend(title="Region Status"))
    ggsave(paste0("/blue/guralnick/millerjared/BoCP/outputs/occur-data-range-status-maps/full-occ-range-map/", accepted_name_file_style, ".jpg"), width = 10, height = 10)
# Native range adjusted map
    # Extract bounding box
    xmin <- occ_native_data_bbox["xmin"]
    xmax <- occ_native_data_bbox["xmax"]
    ymin <- occ_native_data_bbox["ymin"]
    ymax <- occ_native_data_bbox["ymax"]
    
    # Apply buffer based on latitude values
    lat_buffer <- ifelse(abs(ymax) - abs(ymin) > 10, 0.5, 3)
    lon_buffer <- ifelse(abs(xmin) - abs(xmax) > 10, 0.5, 3)
    
    # Expand bounding box
    expanded_bbox <- c(
      xmin - lon_buffer,  # Longitude buffer (constant)
      ymin - lat_buffer,  # Latitude buffer (conditional)
      xmax + lon_buffer,  
      ymax + lat_buffer
    )
    
    # Plot with adjusted bounding box
    ggplot() +
      geom_sf(region_fills, mapping = aes(fill = regionStatus)) +
      scale_fill_manual(values = color_mapping) +
      geom_sf(occ_data_statuses, mapping = aes(), fill = "darkorange", alpha = 0.9) + 
      coord_sf(xlim = c(expanded_bbox[1], expanded_bbox[3]), 
               ylim = c(expanded_bbox[2], expanded_bbox[4])) +
      theme_bw() +
      ggtitle(paste("Adjusted Window Occurrence Records for", accepted_name)) + 
      xlab("Longitude") + 
      ylab("Latitude") +
      theme(plot.title = element_text(hjust = 0.5)) +
      guides(fill=guide_legend(title="Region Status"))
    
    ggsave(paste0("/blue/guralnick/millerjared/BoCP/outputs/occur-data-range-status-maps/adjusted-native-occ-range-map/", accepted_name_file_style, ".jpg"), width = 10, height = 10)
  } # end of if there is more than one occur record make the map...
  
  ## remove the geom col before writing out these data
occ_data_statuses <- occ_data_statuses %>% 
    st_drop_geometry() %>% 
    select(tempID, wcvpRangeStatus)
  
  occ_data <- occ_data %>% 
    left_join(occ_data_statuses, by = "tempID") %>% 
    distinct() %>% 
    select(-tempID)
  
} else{ # deal with the cases of ncbi only name
  occ_data <- occ_data %>% 
    mutate(wcvpRangeStatus = NA) # for ncbi only names
}




# reorder df for flags to be at the end
cols <- colnames(occ_data)
cols <- cols[cols != "taxonomicExactMatch"]  # Remove 'taxonomicExactMatch' from its original position
index <- which(cols == "roundedLongitude")   # Find position of 'roundedLongitude'
new_order <- append(cols, "taxonomicExactMatch", after = index)  # Insert at correct position

occ_data <- occ_data[, ..new_order]  # Reorder the data.table
# write out the data 
print(paste("Writing data for", accepted_name))
fwrite(occ_data, paste0("/blue/guralnick/millerjared/BoCP/data/processed/flagged-species-occs/", accepted_name_file_style, ".csv"))
