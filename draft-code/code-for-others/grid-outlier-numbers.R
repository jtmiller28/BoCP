### Title: Grid Outlier Numbers
### Author: JT Miller
### Date: 04/18/2025

## Load Libraries 
library(tidyverse)
library(data.table)
library(sf)


## Set up Array tasks for SLURM scheduler
start_num <- as.numeric(Sys.getenv("START_NUM"))
task_id <- as.numeric(start_num)
part <- paste0("part", task_id)

### A simple script for determining how much data we're throwing out by discluding outliers from the grids 
## pull out species for this task 
species_files <- list.files("/blue/guralnick/millerjared/BoCP/data/processed/fully-flagged-data/")
# find finished names
finished <- fread("/blue/guralnick/millerjared/BoCP/data/processed/grid-project/sp-occs-outliers.csv")
finished <- gsub(" ", "-", finished$species)
finished <- paste0(finished, ".csv")
unfinished <- setdiff(species_files, finished)
species_file_style <-unfinished[task_id]
species_rds_style <- gsub(".csv", ".rds", species_file_style)
species <- gsub("-", " ", species_file_style)
species <- gsub(".csv", "", species)

## There exists data where BoCP does not have any records, g

## Read in relevant data
# point occurrences from BoCP 
point_occs <- fread(paste0("/blue/guralnick/millerjared/BoCP/data/processed/fully-flagged-data/", species_file_style))

point_occs <- point_occs %>% 
  filter(taxonomicExactMatch == TRUE) %>% 
  filter(wgs84Datum == TRUE) %>% 
  filter(coordinateIssue == FALSE) %>% 
  filter(trueCoordsWithheld == FALSE) %>% 
  filter(wcvpRangeStatus == "native" | wcvpRangeStatus == "introduced") %>% 
  filter(validRecord == TRUE) %>% 
  filter(equalLatLon == FALSE) %>% 
  filter(zeroCoords == FALSE) %>% 
  filter(capitalCoord == FALSE) %>% 
  filter(centroidCoord == FALSE) %>% 
  filter(inOceanCoord == FALSE) %>% 
  filter(inGBIFHeadquarters == FALSE) %>% 
  filter(inInstitutionBounds == FALSE)
if(nrow(point_occs) > 0){
  # Deduplicate data where applicable 
  point_occs <- point_occs %>%
    filter(!amongAggDuplicate | is.na(AggDuplicateGroupID)) %>% # Keep non-duplicates
    bind_rows( # Add in the lowest-ranked duplicates
      point_occs %>%
        filter(amongAggDuplicate) %>%
        group_by(AggDuplicateGroupID) %>%
        filter(AggDuplicateRank == min(AggDuplicateRank, na.rm = TRUE)) %>%
        slice(1) %>%  # if tied for lowest rank, just take the first
        ungroup()
    )
  
  point_occs <- point_occs %>%
    filter(!specimenDuplicate | is.na(specimenDuplicateGroupID)) %>% # Keep non-duplicates
    bind_rows( # Add in the lowest-ranked duplicates
      point_occs %>%
        filter(specimenDuplicate) %>%
        group_by(specimenDuplicateGroupID) %>%
        filter(specimenDuplicateRank == min(specimenDuplicateRank, na.rm = TRUE)) %>%
        slice(1) %>%  # if tied for lowest rank, just take the first
        ungroup()
    )
}

if(nrow(point_occs) > 0){
  point_occs_50 <- point_occs  %>% 
    filter(coordinateUncertaintyInMeters <= 25000 | is.na(coordinateUncertaintyInMeters)) # we're just going with half as a rule
  point_occs_25 <- point_occs %>% 
    filter(coordinateUncertaintyInMeters <= 12500 | is.na(coordinateUncertaintyInMeters))
}
# finally assign to sf object 
if(nrow(point_occs) > 0){
if(nrow(point_occs_50) > 0){
  point_occs_50 <- st_as_sf(point_occs_50, coords = c("roundedLongitude", "roundedLatitude"), crs = 4326)
} else {
  point_occs_50 <- NULL
}
if(nrow(point_occs_25) > 0){
  point_occs_25 <- st_as_sf(point_occs_25, coords = c("roundedLongitude", "roundedLatitude"), crs = 4326)
} else {
  point_occs_25 <- NULL
}
} else{
  point_occs_50 <- NULL
  point_occs_25 <- NULL
}

## Load in botanical regions 
bot_regions <- read_sf("/blue/guralnick/millerjared/BoCP/data/raw/level3-wgsrpd/level3.shp")

## grab relevant regions for Israel's analysis (contiguous US and Canada)
na_string <- c("ALA", "ABT", "ASK", "ARI", "ARK", "BRC", "CAL", "COL", "CNT", "DEL",
               "GEO", "FLA", "IDA", "IOW", "ILL","INI", "KAN", "MAN", "LOU", "KTY", "MAI",
               "MNT", "MIN", "MIC", "MAS", "MSO", "MSI", "MRY",  "NDA", "NCA", "NBR", "NEV", "NEB", "NFL", "NUN", "NSC", 
               "NWT", "NWM", "OHI", "NWJ", "NWH", "NWY", "ORE", "ONT", "OKL", "PEN", "PEI", "QUE",
               "RHO", "SAS", "SDA", "SCA", "TEX", "TEN", "UTA", "VRG", "VER", "WAS", "WIS", "WDC", "WVA",
               "WYO", "LAB", "YUK") # removed "MXE","MXN","MXC",  "MXG", "MXS", "MXT", 

# Filter regions to remove
na_regions <- filter(bot_regions, LEVEL3_COD %in% na_string)


## First lets work with the Nature Serve (NS) data, which is in a flavor of Albers Equal Area grids that I've extracted the exact projection using QGIS tools on the file called Element_occurrences_01_2023Update provided by Israel
crs_102008 <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 
               +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"

# reproj our na regions to this projection 
na_regions_aea <- st_transform(na_regions, crs = crs_102008)

# create a 50x50km grid that matches the extent of the NS data, extent was retrieved manually by reviewing the same file mentioned above in QGIS
x_min <- -5105000.0000000000000000
y_min <- -2905000.0000000000000000
x_max<- 3045000.0001220712438226
y_max <- 4645000.0001220703125000
ns_bbox <- st_bbox(c(xmin = x_min, ymin = y_min, xmax = x_max, ymax = y_max), crs = crs_102008)
ns_grid <- st_make_grid(ns_bbox, cellsize = c(50000, 50000), what = "polygons", square = TRUE)

# clip the grid to only include the maximum spatial extent of na_regions_aea 
ns_grid_sf <- st_sf(geometry = ns_grid)
clipped_ns_grid <- st_intersection(ns_grid_sf, na_regions_aea)

# take occ data & subset by this clipped grids extent boundaries 

if(!is.null(point_occs_50)){
  point_occs_50 <- st_transform(point_occs_50, crs = st_crs(clipped_ns_grid))
  point_occs_50 <- point_occs_50[clipped_ns_grid, ]
}

if(!is.null(point_occs_25)){
  point_occs_25 <- st_transform(point_occs_25, crs = st_crs(clipped_ns_grid))
  point_occs_25 <- point_occs_25[clipped_ns_grid, ]
}

# report outliers as a table for this intended study region
outlier_report <- data.frame(
  species = species,
  num_all_occs_50 = ifelse(!is.null(point_occs_50), nrow(point_occs_50), 0),
  num_non_outliers_50 = ifelse(!is.null(point_occs_50), nrow(filter(point_occs_50, distOutlier == FALSE)), 0),
  num_outliers_50 = ifelse(!is.null(point_occs_50), nrow(filter(point_occs_50, distOutlier == TRUE)), 0),
  num_all_occs_25 = ifelse(!is.null(point_occs_25), nrow(point_occs_25), 0), 
  num_non_outliers_25 = ifelse(!is.null(point_occs_25), nrow(filter(point_occs_25, distOutlier == FALSE)), 0),
  num_outliers_25 = ifelse(!is.null(point_occs_25), nrow(filter(point_occs_25, distOutlier == TRUE)), 0))
  fwrite(outlier_report, file = "/home/millerjared/blue_guralnick/millerjared/BoCP/data/processed/grid-project/sp-occs-outliers.csv", 
         append = TRUE, col.names = !file.exists("/home/millerjared/blue_guralnick/millerjared/BoCP/data/processed/grid-project/sp-occs-outliers.csv"))
  

