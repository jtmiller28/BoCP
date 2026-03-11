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
species_files <- list.files("/blue/guralnick/millerjared/BoCP/data/processed/fully-flagged-data/")
finished_files <- list.files("/blue/guralnick/millerjared/BoCP/data/processed/grid-project/worldwide/raster-50x50km-masked/")
finished_files <- gsub("-50km-masked.tif", ".csv", finished_files)
unfinished_files <- setdiff(species_files, finished_files)
species_file_style <- unfinished_files[task_id]
species_tif50_style <- gsub(".csv", "-50km.tif", species_file_style)
species_tif25_style <- gsub(".csv", "-25km.tif", species_file_style)
species <- gsub("-", " ", species_file_style)
species <- gsub(".csv", "", species)

sp_50_rast <- terra::rast(paste0("/blue/guralnick/millerjared/BoCP/data/processed/grid-project/worldwide/raster-50x50km/", species_tif50_style))
sp_25_rast <- terra::rast(paste0("/blue/guralnick/millerjared/BoCP/data/processed/grid-project/worldwide/raster-25x25km/", species_tif25_style))

# reproj species rasters into an easier viewable mollweide center
moll_proj_center <- "+proj=moll +lon_0=-100 +datum=WGS84 +units=m +no_defs"
sp_50_rast <- project(sp_50_rast, crs(moll_proj_center))
sp_25_rast <- project(sp_25_rast, crs(moll_proj_center))
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

### Extract extent and calc dimensions
bot_regions <- terra::vect("/blue/guralnick/millerjared/BoCP/data/raw/level3-wgsrpd/level3.shp")
### project to mollweide 
moll_proj <- "+proj=moll +datum=WGS84 +units=m +no_defs"
bot_regions_moll <- project(bot_regions, moll_proj)
### Extract extent and calc dimensions
moll_ext <- ext(bot_regions_moll)
### Set res to 50,000 m = 50km & 25,000 m = 25km
res_50_km <- 50000 # meters
res_25_km <- 25000 # meters
### Compute the num of rows and cols 
ncol_50 <- ceiling((xmax(moll_ext) - xmin(moll_ext)) / res_50_km)
nrow_50 <- ceiling((ymax(moll_ext) - ymin(moll_ext)) / res_50_km)
ncol_25 <- ceiling((xmax(moll_ext) - xmin(moll_ext)) / res_25_km)
nrow_25 <- ceiling((ymax(moll_ext) - ymin(moll_ext)) / res_25_km)
### Generate the template rasters for 50x50km and 25x25km mollweide rasters of earth
moll_raster_50 <- terra::rast(ncols = ncol_50, nrows = nrow_50, extent = moll_ext, crs = moll_proj)
moll_raster_25 <- terra::rast(ncols = ncol_25, nrows = nrow_25, extent = moll_ext, crs = moll_proj)
### Create grid cells of the same dimensions
moll_50_bbox <- st_bbox(bot_regions_moll)
moll_50_grid <- st_make_grid(moll_50_bbox, 
                             cellsize = c(49991.21, 49916.21))
moll_25_bbox <- st_bbox(bot_regions_moll)
moll_25_grid <- st_make_grid(moll_25_bbox, 
                             cellsize = c(24995.61, 24993.21))
## First lets work with the Nature Serve (NS) data, which is in a flavor of Albers Equal Area grids that I've extracted the exact projection using QGIS tools on the file called Element_occurrences_01_2023Update provided by Israel
# crs_102008 <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 
#                +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"

# reproj our na regions to this projection 
#na_regions_aea <- st_transform(na_regions, crs = crs_102008)

# create a 50x50km grid that matches the extent of the NS data, extent was retrieved manually by reviewing the same file mentioned above in QGIS
# x_min <- -5105000.0000000000000000
# y_min <- -2905000.0000000000000000
# x_max<- 3045000.0001220712438226
# y_max <- 4645000.0001220703125000
# ns_bbox <- st_bbox(c(xmin = x_min, ymin = y_min, xmax = x_max, ymax = y_max), crs = crs_102008)
# ns_grid <- st_make_grid(ns_bbox, cellsize = c(50000, 50000), what = "polygons", square = TRUE)

# clip the grid to only include the maximum spatial extent of na_regions_aea 
# ns_grid_sf <- st_sf(geometry = ns_grid)
# clipped_ns_grid <- st_intersection(ns_grid_sf, na_regions_aea)
grid_50_sf <- st_sf(geometry = moll_50_grid)
grid_25_sf <- st_sf(geometry = moll_25_grid)
clipped_50_grid <- st_intersection(grid_50_sf, na_regions)
# rematch the grid to the mollweide proj
clipped_ns_grid <- st_transform(clipped_ns_grid, moll_proj_center)

# convert clipped nature serve grid to SpatVector 
clipped_ns_grid_vec <- vect(clipped_ns_grid)

# reproj the clipped nature serve grid to match that of the rasters
clipped_ns_grid_vec <- project(clipped_ns_grid_vec, crs(sp_50_rast))

# clip via cropping
sp_50_rast_cropped <- crop(sp_50_rast, clipped_ns_grid_vec)
sp_25_rast_cropped <- crop(sp_25_rast, clipped_ns_grid_vec)

# Mask to the exact dilimitation of shapes boundary
sp_50_rast_masked <- mask(sp_50_rast_cropped, clipped_ns_grid_vec)
sp_25_rast_masked <- mask(sp_25_rast_cropped, clipped_ns_grid_vec)


# plot to check
r_df <- as.data.frame(sp_50_rast_masked, xy = TRUE, na.rm = TRUE)
r_df$presence <- as.factor(r_df$presence)
na_regions_rastProj <- st_transform(na_regions_aea, crs = crs(sp_50_rast_cropped))
ggplot() + 
  geom_raster(r_df, mapping = aes(x=x, y=y, fill = presence)) + 
  scale_fill_viridis_d(na.value = "transparent", name = "presence") + 
  geom_sf(na_regions_rastProj, fill = NA, color = "black", size = 0.3, mapping = aes()) +
  theme_minimal() + 
  ggtitle("50x50km Presence Cells for ", species) + 
  xlab("Longitude") + 
  ylab("Latitude")

ggsave(paste0("/blue/guralnick/millerjared/BoCP/data/processed/grid-project/worldwide/raster-species-PA-visuals/species-50x50/", gsub(".csv", "", species_file_style), "-50x50.png"), width = 10, height = 10)

# write out the raster object
writeRaster(sp_50_rast_masked, paste0("/blue/guralnick/millerjared/BoCP/data/processed/grid-project/worldwide/raster-50x50km-masked/", gsub(" ", "-", species), "-50km-masked.tif"), overwrite = TRUE)
writeRaster(sp_25_rast_masked, paste0("/blue/guralnick/millerjared/BoCP/data/processed/grid-project/worldwide/raster-25x25km-masked/", gsub(" ", "-", species), "-25km-masked.tif"), overwrite = TRUE)

# add species to the r_df, then write out for community level visualization.
# first run only, create a background
#r_df$species <- "background" # run only first time
r_df$species <- species

# remove absences to create a less bulky file
r_df <- r_df %>% filter(presence == 1) # comment off for first background run
fwrite(r_df, file = "/home/millerjared/blue_guralnick/millerjared/BoCP/data/processed/grid-project/worldwide/raster-species-PA-visuals/community-richness-50x50/community-sp-PA-makeup.csv", 
       append = TRUE, col.names = !file.exists("/home/millerjared/blue_guralnick/millerjared/BoCP/data/processed/grid-project/worldwide/raster-species-PA-visuals/community-richness-50x50/community-sp-PA-makeup.csv"))

# plot to check
r_df <- as.data.frame(sp_25_rast_masked, xy = TRUE, na.rm = TRUE)
r_df$presence <- as.factor(r_df$presence)
na_regions_rastProj <- st_transform(na_regions_aea, crs = st_crs(sp_25_rast_cropped))
ggplot() + 
  geom_raster(r_df, mapping = aes(x=x, y=y, fill = presence)) + 
  scale_fill_viridis_d(na.value = "transparent", name = "presence") + 
  geom_sf(na_regions_rastProj, fill = NA, color = "black", size = 0.3, mapping = aes()) +
  theme_minimal() + 
  ggtitle("25x25km Presence Cells for ", species) + 
  xlab("Longitude") + 
  ylab("Latitude")

ggsave(paste0("/blue/guralnick/millerjared/BoCP/data/processed/grid-project/worldwide/raster-species-PA-visuals/species-25x25/", gsub(".csv", "", species_file_style), "-25x25.png"), width = 10, height = 10)

# add species to the r_df, then write out for community level visualization.
#r_df$species <- "background" # run only first time as background layer
r_df$species <- species
# remove absences to create a less bulky file
r_df <- r_df %>% filter(presence == 1) # comment off for first background run

fwrite(r_df, file = "/home/millerjared/blue_guralnick/millerjared/BoCP/data/processed/grid-project/worldwide/raster-species-PA-visuals/community-richness-25x25/community-sp-PA-makeup.csv", 
       append = TRUE, col.names = !file.exists("/home/millerjared/blue_guralnick/millerjared/BoCP/data/processed/grid-project/worldwide/raster-species-PA-visuals/community-richness-25x25/community-sp-PA-makeup.csv"))


