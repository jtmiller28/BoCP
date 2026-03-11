### Title: Clean Outliers To Geographic Scope 
### Author: JT Miller
### Date: 05/02/2025

### A script to clean up the output rasters to fit the geographic limits of the project: Canada + USA. 

# load libraries
library(sf)
library(tidyverse)
library(data.table)
library(terra)

# Set up Array Logic 
start_num <- as.numeric(Sys.getenv("START_NUM"))
task_id <- as.numeric(start_num)
part <- paste0("part", task_id)

## pull out species for this task
# species_files <- list.files("/blue/guralnick/millerjared/BoCP/data/processed/fully-flagged-data/")
# finished_files <- list.files("/blue/guralnick/millerjared/BoCP/data/processed/grid-project/worldwide/raster-intro-native-50x50km-masked/")
# finished_files <- gsub("-50km-masked.tif", ".csv", finished_files)
# unfinished_files <- setdiff(species_files, finished_files)
# saveRDS(unfinished_files, "/blue/guralnick/millerjared/BoCP/data/processed/grid-project/worldwide/unfinished-intro-native-mask-names.rds")
  unfinished_files <- readRDS("/blue/guralnick/millerjared/BoCP/data/processed/grid-project/worldwide/unfinished-intro-native-mask-names.rds")
  species_file_style <- unfinished_files[task_id]
  species_tif50_style <- gsub(".csv", "-50km.tif", species_file_style)
  species_tif25_style <- gsub(".csv", "-25km.tif", species_file_style)
  species <- gsub("-", " ", species_file_style)
  species <- gsub(".csv", "", species)
  
  sp_50_rast <- terra::rast(paste0("/blue/guralnick/millerjared/BoCP/data/processed/grid-project/worldwide/raster-intro-native-50x50km/", species_tif50_style))
  sp_25_rast <- terra::rast(paste0("/blue/guralnick/millerjared/BoCP/data/processed/grid-project/worldwide/raster-intro-native-25x25km/", species_tif25_style))
  
  # reproj species rasters into an easier viewable mollweide center
  moll_proj_center <- "+proj=moll +lon_0=-100 +datum=WGS84 +units=m +no_defs"
  sp_50_rast <- project(sp_50_rast, crs(moll_proj_center), method = "near")
  sp_25_rast <- project(sp_25_rast, crs(moll_proj_center), method = "near")
  ### build clipped grid for delimiting these rasters
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
  # Transform the na regions to the mollweide view proj
  na_regions <- st_transform(na_regions, crs(moll_proj_center))
  # clip via cropping to region of interest
  sp_50_rast_cropped <- crop(sp_50_rast, na_regions)
  sp_25_rast_cropped <- crop(sp_25_rast, na_regions)
  # Mask to the exact dilimitation of region of interest
  sp_50_rast_masked <- mask(sp_50_rast_cropped, na_regions)
  sp_25_rast_masked <- mask(sp_25_rast_cropped, na_regions)
  # reform spatRaster into dataframe for plotting
  sp_50_df <- as.data.frame(sp_50_rast_masked, xy = TRUE)
  sp_50_df$presence <- as.factor(sp_50_df$presence)
  sp_25_df <- as.data.frame(sp_25_rast_masked, xy = TRUE)
  sp_25_df$presence <- as.factor(sp_25_df$presence)
  # Visualize using raster 
  ggplot() +
    geom_tile(sp_50_df, mapping = aes(x = x, y = y, fill = presence)) +
    scale_fill_viridis_d(na.value = "transparent", name = "presence") +
    geom_sf(na_regions, fill = NA, color = "black", size = 0.3, mapping = aes()) +
    theme_minimal() + 
    ggtitle("50x50km Presence Cells for ", species) + 
    xlab("Longitude") + 
    ylab("Latitude")
  ggsave(paste0("/blue/guralnick/millerjared/BoCP/data/processed/grid-project/worldwide/raster-species-PA-visuals/species-intro-native-50x50/", gsub(".csv", "", species_file_style), "-50x50.png"), width = 10, height = 10)
  ggplot() +
    geom_tile(sp_25_df, mapping = aes(x = x, y = y, fill = presence)) +
    scale_fill_viridis_d(na.value = "transparent", name = "presence") +
    geom_sf(na_regions, fill = NA, color = "black", size = 0.3, mapping = aes()) +
    theme_minimal() + 
    ggtitle("25x25km Presence Cells for ", species) + 
    xlab("Longitude") + 
    ylab("Latitude")
  ggsave(paste0("/blue/guralnick/millerjared/BoCP/data/processed/grid-project/worldwide/raster-species-PA-visuals/species-intro-native-25x25/", gsub(".csv", "", species_file_style), "-25x25.png"), width = 10, height = 10)
  # write out the raster object
  writeRaster(sp_50_rast_masked, paste0("/blue/guralnick/millerjared/BoCP/data/processed/grid-project/worldwide/raster-intro-native-50x50km-masked/", gsub(" ", "-", species), "-50km-masked.tif"), overwrite = TRUE)
  writeRaster(sp_25_rast_masked, paste0("/blue/guralnick/millerjared/BoCP/data/processed/grid-project/worldwide/raster-intro-native-25x25km-masked/", gsub(" ", "-", species), "-25km-masked.tif"), overwrite = TRUE)
  # add species to the r_df, then write out for community level visualization.
  # first run only, create a background
  # sp_50_df$species <- "background" # run only first time
  # sp_25_df$species <- "background" # run only first time
  sp_50_df$species <- species # comment off for first run
  sp_25_df$species <- species # comment off for first run
  # remove absences to create a less bulky file
  sp_50_df <- sp_50_df %>% filter(presence == 1) # comment off for first background run
  sp_25_df <- sp_25_df %>% filter(presence == 1) # comment off for first background run
  fwrite(sp_50_df, file = "/home/millerjared/blue_guralnick/millerjared/BoCP/data/processed/grid-project/worldwide/raster-species-PA-visuals/community-richness-50x50/community-sp-PA-makeup-intro-native.csv", 
         append = TRUE, col.names = !file.exists("/home/millerjared/blue_guralnick/millerjared/BoCP/data/processed/grid-project/worldwide/raster-species-PA-visuals/community-richness-50x50/community-sp-PA-makeup-intro-native.csv"))
  fwrite(sp_25_df, file = "/home/millerjared/blue_guralnick/millerjared/BoCP/data/processed/grid-project/worldwide/raster-species-PA-visuals/community-richness-25x25/community-sp-PA-makeup-intro-native.csv", 
         append = TRUE, col.names = !file.exists("/home/millerjared/blue_guralnick/millerjared/BoCP/data/processed/grid-project/worldwide/raster-species-PA-visuals/community-richness-25x25/community-sp-PA-makeup-intro-native.csv"))