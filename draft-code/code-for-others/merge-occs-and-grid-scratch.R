### Title: Merge Occ and Grid data
### Author: JT Miller
### Date: 04/18/2025

# Purpose, find how to effectively integrate israel's datasets and our point occurrences at 50x50km and 25x25km, and provide figures to show how well modeling might work given these data

## Load Libraries 
library(tidyverse)
library(data.table)
library(sf)

## Load in botanical regions 
bot_regions <- read_sf("./data/raw/level3-wgsrpd/level3.shp")

## grab relevant regions for Israel's analysis (contiguous US and Canada)
na_string <- c("ALA", "ABT", "ASK", "ARI", "ARK", "BRC", "CAL", "COL", "CNT", "DEL",
               "GEO", "FLA", "IDA", "IOW", "ILL","INI", "KAN", "MAN", "LOU", "KTY", "MAI",
               "MNT", "MIN", "MIC", "MAS", "MSO", "MSI", "MRY",  "NDA", "NCA", "NBR", "NEV", "NEB", "NFL", "NUN", "NSC", 
               "NWT", "NWM", "OHI", "NWJ", "NWH", "NWY", "ORE", "ONT", "OKL", "PEN", "PEI", "QUE",
               "RHO", "SAS", "SDA", "SCA", "TEX", "TEN", "UTA", "VRG", "VER", "WAS", "WIS", "WDC", "WVA",
               "WYO", "LAB", "YUK") # removed "MXE","MXN","MXC",  "MXG", "MXS", "MXT", 

# Filter regions to remove
na_regions <- filter(bot_regions, LEVEL3_COD %in% na_string)

# make plot thats more visually easy to see as an extent plot 
ggplot() + 
  geom_sf(data = bot_regions) + 
  geom_sf(data = na_regions, fill = "darkred", color = "black") + 
  ggtitle("USA and Canada for Endemic Pixel Project")

## First lets work with the Nature Serve (NS) data, which is in a flavor of Albers Equal Area grids that I've extracted the exact projection using QGIS tools on the file called Element_occurrences_01_2023Update provided by Israel
crs_102008 <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 
               +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"

# reproj our na regions to this projection 
na_regions_aea <- st_transform(na_regions, crs = crs_102008)

# plot to check if things are askew 
ggplot() + 
  geom_sf(data = na_regions_aea, fill = "darkred", color = "black") + 
  ggtitle("USA and Canada for Endemic Pixel Project")

# create a 50x50km grid that matches the extent of the NS data, extent was retrieved manually by reviewing the same file mentioned above in QGIS
x_min <- -5105000.0000000000000000
y_min <- -2905000.0000000000000000
x_max<- 3045000.0001220712438226
y_max <- 4645000.0001220703125000
ns_bbox <- st_bbox(c(xmin = x_min, ymin = y_min, xmax = x_max, ymax = y_max), crs = crs_102008)
ns_grid <- st_make_grid(ns_bbox, cellsize = c(50000, 50000), what = "polygons", square = TRUE)

# plot to see if this looks about right 
ggplot() + 
  geom_sf(data = na_regions_aea, fill = "darkred", color = "black") + 
  geom_sf(data = ns_grid, fill = "black", alpha = 0.7) + 
  ggtitle("USA and Canada with 50x50km grid \n for Endemic Pixel Project")

# now load in Israel's NS centroid data per species 
nature_serve_centroids <- read_sf("/blue/guralnick/millerjared/BoCP/data/raw/grid-data/na-rare-plant-grid-data/NA spatial phylogenetics/tns_points.shp")

# reproj these to the AEA proj 
nature_serve_centroids <- st_set_crs(nature_serve_centroids, crs_102008)

# create a test set of names (determined prior by a summary of 10-20 unique centroid points)
test_names <- c("Tumamoca macdougalii", "Heuchera eastwoodiae", "Lachnanthes caroliniana", 
                "Agalinis filicaulis", "Carex dasycarpa", "Rhexia salicifolia", "Parnassia asarifolia",
                "Viburnum lentago", "Echinacea pallida", "Trillium sessile")

# for(i in 1:length(test_names)){
#   nature_serve_test <- nature_serve_centroids %>%  # extract species level centroid pts
#     filter(SNAME == test_names[i]) %>% 
#     st_transform(crs = st_crs(crs_102008))
#   bbox <- st_bbox(nature_serve_test) # create a bbox to view their extent
#   ggplot() +  
#     geom_sf(na_regions_aea, mapping = aes()) + 
#     geom_sf(ns_grid, mapping = aes(), alpha = 0.2) +
#     geom_sf(nature_serve_test, mapping = aes()) +
#     coord_sf(xlim = c(bbox["xmin"] - 10000, bbox["xmax"] + 10000), # add slight buffers so its viewable
#              ylim = c(bbox["ymin"] - 10000, bbox["ymax"] + 10000)) +
#     ggtitle(paste0("Nature Serve 50x50km presence data for ", test_names[i]))
#   ggsave(paste0("/blue/guralnick/millerjared/BoCP/outputs/example-grid-extractions/", test_names[i], "-NS-sample-zoom.png"), width = 10, height = 10)
#   
#   ggplot() +  
#     geom_sf(na_regions_aea, mapping = aes()) + 
#     geom_sf(ns_grid, mapping = aes(), alpha = 0.1) +
#     geom_sf(nature_serve_test, mapping = aes()) +
#     coord_sf(xlim = c(bbox["xmin"] - 1000000, bbox["xmax"] + 1000000), # add slight buffers so its viewable
#              ylim = c(bbox["ymin"] - 1000000, bbox["ymax"] + 1000000)) +
#     ggtitle(paste0("Nature Serve 50x50km presence data for ", test_names[i]))
#   ggsave(paste0("/blue/guralnick/millerjared/BoCP/outputs/example-grid-extractions/", test_names[i], "-NS-sample-.png"), width = 10, height = 10)
# }

## Redo these grids as presence absence cells
library(sf)
library(terra)

# reformat our grid as an sf object
ns_grid_sf <- st_sf(id = 1:length(ns_grid), geometry = ns_grid)

# test case 
nature_serve_test <- nature_serve_centroids %>%  # extract species level centroid pts
  filter(SNAME == test_names[1]) 

## make a test 25x25 grid to make sure things are operating correctly
ns_extent <- st_bbox(ns_grid_sf)
ns_grid_25 <- st_make_grid(ns_extent, 
                           crs = st_crs(ns_grid),
                           cellsize = c(25000,25000), 
                           square = TRUE)
# preform a spatial join to determine presence within a cell
joined <- st_join(nature_serve_test, ns_grid_sf, join = st_within)

# get p/a if the grid contains a centroid 
presence <- unique(joined$id)
sp_grid <- ns_grid_sf # make a copy for this particular sp
sp_grid <- sp_grid %>% # label presence absence 
  dplyr::mutate(presence = ifelse(id %in% presence, 1, 0)) %>% 
  dplyr::mutate(presence = as.factor(presence))

# export as raster and shapefiles to give flexibility down the line
## raster
# create an empty raster template, contains extent of original grid at 50x50km, assuure crs also matches
r_template <- rast(ext(sp_grid), resolution = 50000, crs = st_crs(sp_grid)$wkt)
# maps presence to raster cells, choose fun = "max" to take 1 if available
r_presence <- rasterize(vect(sp_grid), r_template, field = "presence", fun = "max")
# validate whether this worked correctly via a visual
# Convert raster to points
r_points <- as.points(r_presence)
r_points_sf <- st_as_sf(r_points)
bbox <- st_bbox(filter(r_points_sf, presence == 1))
# ggplot
ggplot() +
  geom_sf(data = nature_serve_test, color = "red", size = 2) +
  geom_sf(data = r_points_sf, aes(color = as.factor(presence)), size = 1) +
  geom_sf(data = ns_grid_sf, fill = NA, color = "grey") +
  scale_color_manual(values = c("white", "darkgreen")) +
  coord_sf(xlim = c(bbox["xmin"] - 100000, bbox["xmax"] + 100000), # add slight buffers so its viewable
           ylim = c(bbox["ymin"] - 100000, bbox["ymax"] + 100000)) +
  ggtitle(paste0("Presence/Absence raster vs. points for ", test_names[1])) +
  theme_minimal()

# Looks about right. Write it out as a raster & Shapefile 
# write out 
writeRaster(r_presence, paste0("/blue/guralnick/millerjared/BoCP/data/processed/grid-project/raster-50x50km/raster-PA-", test_names[1], ".tif"), overwrite = TRUE)
## shapefile
st_write(sp_grid, paste0("/blue/guralnick/millerjared/BoCP/data/processed/grid-project/shp-50x50km/", test_names[1], "-PA", ".shp"), append = FALSE)

### Now try a resample over a template.
r_50 <- r_presence # create a var name for tracking 
template_25 <- rast(ext = ext(r_50), resolution = 25000, crs = crs(r_50))
r_50_to_25 <- resample(r_50, template_25, method = "max")

# check via plotting to see if this worked correctly
r_points <- as.points(r_50_to_25)
r_points_sf <- st_as_sf(r_points)
bbox <- st_bbox(filter(r_points_sf, presence == 1))
ggplot() +
  geom_sf(data = nature_serve_test, color = "red", size = 2) +
  geom_sf(data = r_points_sf, aes(color = as.factor(presence)), size = 1) +
  geom_sf(data = ns_grid_25, fill = NA, color = "grey") +
  scale_color_manual(values = c("white", "darkgreen")) +
  coord_sf(xlim = c(bbox["xmin"] - 100000, bbox["xmax"] + 100000), # add slight buffers so its viewable
           ylim = c(bbox["ymin"] - 100000, bbox["ymax"] + 100000)) +
  ggtitle(paste0("Presence/Absence raster vs. points for ", test_names[1])) +
  theme_minimal()

# as intended, 50x50km cells are divided up. 

# check to see if canada data works properly
# read in both sets of 10x10km canada centroid data 
canada1_centroids <- read_sf("/blue/guralnick/millerjared/BoCP/data/raw/grid-data/na-rare-plant-grid-data/NA spatial phylogenetics/canada1.shp")
canada2_centroids <- read_sf("/blue/guralnick/millerjared/BoCP/data/raw/grid-data/na-rare-plant-grid-data/NA spatial phylogenetics/canada2.shp")
# standardize fields for both datasets
canada1_centroids <- canada1_centroids %>% 
  select(NNAME) %>% 
  rename(SNAME = NNAME)
canada2_centroids <- canada2_centroids %>% 
  select(GNAME) %>% 
  rename(SNAME = GNAME)
# combine datasets
canada_centroids <- rbind(canada1_centroids, canada2_centroids)
# bring in canada grid
canada_grid <- read_sf("./data/raw/grid-data/na-rare-plant-grid-data/NA spatial phylogenetics/NAm_cell_template_10km/NAm_cell_template_10km.shp")
# standardize that grid
canada_grid <- canada_grid %>% 
  select(-Id, -GRID_ID) %>% 
  mutate(id = 1:n())
# target grid size 
can_extent <- st_bbox(canada_grid)
can_25_grid <- st_make_grid(can_extent, # make a 25x25km grid using the same extent, so that we may plot the transformation later as a check
                            cellsize = c(25000, 25000),
                            crs = st_crs(canada_grid), 
                            square = TRUE)
# combine grid and canada centroid data, on a test name
can_test <- filter(canada_centroids, SNAME == "Solidago fallax")
joined <- st_join(can_test, canada_grid, join = st_intersects)
# get p/a if the grid contains a centroid 
presence <- unique(joined$id)
sp_grid <- canada_grid # make a copy for this particular sp
sp_grid <- sp_grid %>% # label presence absence 
  dplyr::mutate(presence = ifelse(id %in% presence, 1, 0)) %>% 
  dplyr::mutate(presence = as.factor(presence))
# create an empty raster template, contains extent of original grid at 50x50km, assuure crs also matches
r_template <- rast(ext(sp_grid), resolution = 10000, crs = st_crs(sp_grid)$wkt)
# maps presence to raster cells, choose fun = "max" to take 1 if available
r_presence <- rasterize(vect(sp_grid), r_template, field = "presence", fun = "max")
# validate whether this worked correctly via a visual
# Convert raster to points
r_points <- as.points(r_presence)
r_points_sf <- st_as_sf(r_points)
bbox <- st_bbox(filter(r_points_sf, presence == 1))
# ggplot
ggplot() +
  geom_sf(data = can_test, color = "red", size = 2) +
  geom_sf(data = r_points_sf, aes(color = as.factor(presence)), size = 1) +
  geom_sf(data = canada_grid, fill = NA, color = "grey") +
  scale_color_manual(values = c("white", "darkgreen")) +
  coord_sf(xlim = c(bbox["xmin"] - 100000, bbox["xmax"] + 100000), # add slight buffers so its viewable
           ylim = c(bbox["ymin"] - 100000, bbox["ymax"] + 100000)) +
  ggtitle(paste0("Presence/Absence raster vs. points for ")) +
  theme_minimal()

# Looks about right. Write it out as a raster & Shapefile (do later)
r_10 <- r_presence # create a var name for tracking 
template_25 <- rast(ext = ext(r_10), resolution = 25000, crs = crs(r_10))
r_10_to_25 <- resample(r_10, template_25, method = "max")

# check via plotting to see if this worked correctly
r_points <- as.points(r_10_to_25)
r_points_sf <- st_as_sf(r_points)
bbox <- st_bbox(filter(r_points_sf, presence == 1))
ggplot() +
  geom_sf(data = can_test, color = "red", size = 2) +
  geom_sf(data = r_points_sf, aes(color = as.factor(presence)), size = 1) +
  geom_sf(data = can_25_grid, fill = NA, color = "grey") +
  scale_color_manual(values = c("white", "darkgreen")) +
  coord_sf(xlim = c(bbox["xmin"] - 100000, bbox["xmax"] + 100000), # add slight buffers so its viewable
           ylim = c(bbox["ymin"] - 100000, bbox["ymax"] + 100000)) +
  ggtitle(paste0("Presence/Absence raster vs. points for ", test_names[1])) +
  theme_minimal()

### Now combine these methods, accounting for the merged extent.
# make a test case, something with low amount of pts for visualizing
ns_summary <- nature_serve_centroids %>% st_drop_geometry() %>% group_by(SNAME) %>% summarize(n = n())
can_summary <- canada_centroids %>% st_drop_geometry() %>% group_by(SNAME) %>% summarize(n = n())
test <- ns_summary %>% filter(SNAME %in% can_summary$SNAME) %>% filter(n <= 10)
test <- merge(test, can_summary, by = "SNAME")
test <- test %>% filter(n.y <= 10)
# we'll use Alchemilla filicaulis for this test case.

## Set up nature serve data
crs_102008 <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 
               +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
# create a 50x50km grid that matches the extent of the NS data, extent was retrieved manually by reviewing the same file mentioned above in QGIS
x_min <- -5105000.0000000000000000
y_min <- -2905000.0000000000000000
x_max<- 3045000.0001220712438226
y_max <- 4645000.0001220703125000
ns_bbox <- st_bbox(c(xmin = x_min, ymin = y_min, xmax = x_max, ymax = y_max), crs = crs_102008)
ns_grid <- st_make_grid(ns_bbox, cellsize = c(50000, 50000), what = "polygons", square = TRUE)
# reformat our grid as an sf object
ns_grid_sf <- st_sf(id = 1:length(ns_grid), geometry = ns_grid)
# now load in Israel's NS centroid data per species 
nature_serve_centroids <- read_sf("/blue/guralnick/millerjared/BoCP/data/raw/grid-data/na-rare-plant-grid-data/NA spatial phylogenetics/tns_points.shp")
# reproj these to the AEA proj 
nature_serve_centroids <- st_set_crs(nature_serve_centroids, crs_102008)
nature_serve_test <- nature_serve_centroids %>%  # extract species level centroid pts
  filter(SNAME == "Alchemilla filicaulis") 
# preform a spatial join to determine presence within a cell
joined <- st_join(nature_serve_test, ns_grid_sf, join = st_within)
# get p/a if the grid contains a centroid 
presence <- unique(joined$id)
sp_grid <- ns_grid_sf # make a copy for this particular sp
sp_grid <- sp_grid %>% # label presence absence 
  dplyr::mutate(presence = ifelse(id %in% presence, 1, 0)) %>% 
  dplyr::mutate(presence = as.factor(presence))
# create an empty raster template, contains extent of original grid at 50x50km, assuure crs also matches
r_template <- rast(ext(sp_grid), resolution = 50000, crs = st_crs(sp_grid)$wkt)
# maps presence to raster cells, choose fun = "max" to take 1 if available
r_presence <- rasterize(vect(sp_grid), r_template, field = "presence", fun = "max")
r_50 <- r_presence # create a var name for tracking 

## Set up Canada data 
canada1_centroids <- read_sf("/blue/guralnick/millerjared/BoCP/data/raw/grid-data/na-rare-plant-grid-data/NA spatial phylogenetics/canada1.shp")
canada2_centroids <- read_sf("/blue/guralnick/millerjared/BoCP/data/raw/grid-data/na-rare-plant-grid-data/NA spatial phylogenetics/canada2.shp")
# standardize fields for both datasets
canada1_centroids <- canada1_centroids %>% 
  select(NNAME) %>% 
  rename(SNAME = NNAME)
canada2_centroids <- canada2_centroids %>% 
  select(GNAME) %>% 
  rename(SNAME = GNAME)
# combine datasets
canada_centroids <- rbind(canada1_centroids, canada2_centroids)
# bring in canada grid
canada_grid <- read_sf("./data/raw/grid-data/na-rare-plant-grid-data/NA spatial phylogenetics/NAm_cell_template_10km/NAm_cell_template_10km.shp")
# standardize that grid
canada_grid <- canada_grid %>% 
  select(-Id, -GRID_ID) %>% 
  mutate(id = 1:n())
# combine grid and canada centroid data, on a test name
can_test <- filter(canada_centroids, SNAME == "Alchemilla filicaulis")
joined <- st_join(can_test, canada_grid, join = st_intersects)
# get p/a if the grid contains a centroid 
presence <- unique(joined$id)
sp_grid <- canada_grid # make a copy for this particular sp
sp_grid <- sp_grid %>% # label presence absence 
  dplyr::mutate(presence = ifelse(id %in% presence, 1, 0)) %>% 
  dplyr::mutate(presence = as.factor(presence))
# create an empty raster template, contains extent of original grid at 10x10km, assuure crs also matches
r_template <- rast(ext(sp_grid), resolution = 10000, crs = st_crs(sp_grid)$wkt)
# maps presence to raster cells, choose fun = "max" to take 1 if available
r_presence <- rasterize(vect(sp_grid), r_template, field = "presence", fun = "max")
r_10 <- r_presence # rename for variable clarity

## Combine the extents of these rasters
# first check if extents match
crs(r_50) == crs(r_10) # FALSE, therefore we need to reproj
r_50_proj <- project(r_50, r_10, method = "near") # choose the larger grid proj so we can maintain finer grain projection as is. 
crs(r_50_proj) == crs(r_10)
# check how this changed the data:
r_points <- as.points(r_50_proj)
r_points_sf <- st_as_sf(r_points)
bbox <- st_bbox(filter(r_points_sf, presence == 1))
ns_test_proj <- st_transform(nature_serve_test, crs = crs(r_10))
ns_grid_sf_proj <- st_transform(ns_grid_sf, crs = crs(r_10))
# ggplot
ggplot() +
  geom_sf(data = ns_test_proj, color = "red", size = 2) +
  geom_sf(data = r_points_sf, aes(color = as.factor(presence)), size = 1) +
  geom_sf(data = ns_grid_sf_proj, fill = NA, color = "grey") +
  scale_color_manual(values = c("white", "darkgreen")) +
  coord_sf(xlim = c(bbox["xmin"] - 100000, bbox["xmax"] + 100000), # add slight buffers so its viewable
           ylim = c(bbox["ymin"] - 100000, bbox["ymax"] + 100000)) +
  ggtitle(paste0("Presence/Absence raster vs. points for ")) +
  theme_minimal()

# combine the extents of both projections 
combined_extent <- ext(r_50_proj) + ext(r_10)

# create templates of our desired resolutions based on those combined extents 
template_50 <- rast(ext = combined_extent, resolution = 50000, crs = crs(r_10))
template_25 <- rast(ext = combined_extent, resolution = 25000, crs = crs(r_10))

# now run the resampling
r_10_to_50 <- resample(r_10, template_50, method = "max")
r_10_to_25 <- resample(r_10, template_25, method = "max")

# check results 
r_points <- as.points(r_10_to_25)
r_points_sf <- st_as_sf(r_points)
bbox <- st_bbox(filter(r_points_sf, presence == 1))
combined_extent_25_grid <- st_make_grid(combined_extent, 
                                        crs = crs(r_10_to_25), 
                                        cellsize = c(25000, 25000))
can_test_proj <- st_transform(can_test, crs = crs(r_10_to_25))
ggplot() +
  geom_sf(data = can_test_proj, color = "red", size = 2) +
  geom_sf(data = r_points_sf, aes(color = as.factor(presence)), size = 1) +
  geom_sf(data = combined_extent_25_grid, fill = NA, color = "grey") +
  scale_color_manual(values = c("white", "darkgreen")) +
  coord_sf(xlim = c(bbox["xmin"] - 100000, bbox["xmax"] + 100000), # add slight buffers so its viewable
           ylim = c(bbox["ymin"] - 100000, bbox["ymax"] + 100000)) +
  theme_minimal()
# Passes all the visuals, proceed with NS downsizing 
r_50_to_50 <- resample(r_50_proj, template_50, method = "max")
r_50_to_25 <- resample(r_50_proj, template_25, method = "max")

# check results 
r_points <- as.points(r_50_to_25)
r_points_sf <- st_as_sf(r_points)
bbox <- st_bbox(filter(r_points_sf, presence == 1))
combined_extent_25_grid <- st_make_grid(combined_extent, 
                                        crs = crs(r_50_to_25), 
                                        cellsize = c(25000, 25000))
ns_test_proj <- st_transform(nature_serve_test, crs = crs(r_50_to_25))
ggplot() +
  geom_sf(data = ns_test_proj, color = "red", size = 2) +
  geom_sf(data = r_points_sf, aes(color = as.factor(presence)), size = 1) +
  geom_sf(data = combined_extent_25_grid, fill = NA, color = "grey") +
  scale_color_manual(values = c("white", "darkgreen")) +
  coord_sf(xlim = c(bbox["xmin"] - 100000, bbox["xmax"] + 100000), # add slight buffers so its viewable
           ylim = c(bbox["ymin"] - 100000, bbox["ymax"] + 100000)) +
  theme_minimal()

# Passes, lets figure out a way to now add them together. 
# First, lets create a visual of what they look like when together on a 25x25
r_points_ns <- as.points(r_50_to_25)
r_points_sf_ns <- st_as_sf(r_points_ns)
bbox_ns <- st_bbox(filter(r_points_sf_ns, presence == 1))
combined_extent_25_grid <- st_make_grid(combined_extent, 
                                        crs = crs(r_50_to_25), 
                                        cellsize = c(25000, 25000))
ns_test_proj <- st_transform(nature_serve_test, crs = crs(r_50_to_25))

r_points_can <- as.points(r_10_to_25)
r_points_sf_can <- st_as_sf(r_points_can)
bbox_can <- st_bbox(filter(r_points_sf_can, presence == 1))
can_test_proj <- st_transform(can_test, crs = crs(r_10_to_25))
# force factor levels for plotting shennaningans 
r_points_sf_can$presence <- factor(r_points_sf_can$presence, levels = c(0,1))
r_points_sf_ns$presence  <- factor(r_points_sf_ns$presence, levels = c(0,1))
ggplot() +
  geom_sf(data = ns_test_proj, color = "steelblue", size = 2) +
  geom_sf(data = can_test_proj, color = "goldenrod", size = 1.5) +
  geom_sf(data = filter(r_points_sf_can, presence == 1), color = "red", size = 1, alpha = 0.5) +
  geom_sf(data = filter(r_points_sf_ns, presence == 1), color = "green", size = 1, alpha = 0.5) +
  geom_sf(data = combined_extent_25_grid, fill = NA, color = "grey") +
  coord_sf(xlim = c(bbox["xmin"] - 100000, bbox["xmax"] + 100000), # add slight buffers so its viewable
           ylim = c(bbox["ymin"] - 100000, bbox["ymax"] + 100000)) +
  theme_minimal()

# now, attempt to combine these data 
combined_25 <- app(c(r_50_to_25, r_10_to_25), fun = max, na.rm = TRUE)
names(combined_25) <- "presence"
r_points_c <- as.points(combined_25)
r_points_sf_c <- st_as_sf(r_points_c)
bbox_c <- st_bbox(filter(r_points_sf_c, presence == 1))

ggplot() +
  geom_sf(data = r_points_sf_c, aes(color = as.factor(presence)), size = 1) +
  geom_sf(data = combined_extent_25_grid, fill = NA, color = "grey") +
  scale_color_manual(values = c("white", "darkgreen")) +
  coord_sf(xlim = c(bbox_c["xmin"] - 100000, bbox_c["xmax"] + 100000), # add slight buffers so its viewable
           ylim = c(bbox_c["ymin"] - 100000, bbox_c["ymax"] + 100000)) +
  theme_minimal()

## Success, now just check the last case of going from 50 to 50 and 10 to 50
# Second, lets create a visual of what they look like when together on a 50x50
r_points_ns <- as.points(r_50_to_50)
r_points_sf_ns <- st_as_sf(r_points_ns)
bbox_ns <- st_bbox(filter(r_points_sf_ns, presence == 1))
combined_extent_50_grid <- st_make_grid(combined_extent, 
                                        crs = crs(r_50_to_50), 
                                        cellsize = c(50000, 50000))
ns_test_proj <- st_transform(nature_serve_test, crs = crs(r_50_to_50))

r_points_can <- as.points(r_10_to_50)
r_points_sf_can <- st_as_sf(r_points_can)
bbox_can <- st_bbox(filter(r_points_sf_can, presence == 1))
can_test_proj <- st_transform(can_test, crs = crs(r_10_to_50))
# force factor levels for plotting shennaningans 
r_points_sf_can$presence <- factor(r_points_sf_can$presence, levels = c(0,1))
r_points_sf_ns$presence  <- factor(r_points_sf_ns$presence, levels = c(0,1))
ggplot() +
  geom_sf(data = ns_test_proj, color = "steelblue", size = 2, alpha = 0.5) +
  geom_sf(data = can_test_proj, color = "goldenrod", size = 1.5, alpha = 0.5) +
  geom_sf(data = filter(r_points_sf_can, presence == 1), color = "red", size = 1, alpha = 0.9) +
  geom_sf(data = filter(r_points_sf_ns, presence == 1), color = "green", size = 1, alpha = 0.5) +
  geom_sf(data = combined_extent_50_grid, fill = NA, color = "grey") +
  coord_sf(xlim = c(bbox["xmin"] - 100000, bbox["xmax"] + 100000), # add slight buffers so its viewable
           ylim = c(bbox["ymin"] - 100000, bbox["ymax"] + 100000)) +
  theme_minimal()

# combine data 
combined_50 <- app(c(r_50_to_50, r_10_to_50), fun = max, na.rm = TRUE)
names(combined_50) <- "presence"
r_points_c <- as.points(combined_50)
r_points_sf_c <- st_as_sf(r_points_c)
bbox_c <- st_bbox(filter(r_points_sf_c, presence == 1))

ggplot() +
  geom_sf(data = r_points_sf_c, aes(color = as.factor(presence)), size = 1) +
  geom_sf(data = combined_extent_50_grid, fill = NA, color = "grey") +
  scale_color_manual(values = c("white", "darkgreen")) +
  coord_sf(xlim = c(bbox_c["xmin"] - 100000, bbox_c["xmax"] + 100000), # add slight buffers so its viewable
           ylim = c(bbox_c["ymin"] - 100000, bbox_c["ymax"] + 100000)) +
  theme_minimal()

## Also success. 

## Last step, now we need to add occ data from BoCP where applicable, and fill cells that can be filled 
bocp_data_test <- fread("/blue/guralnick/millerjared/BoCP/data/processed/fully-flagged-data/Alchemilla-filicaulis.csv")
test_df_filtered <- bocp_data_test %>% 
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


test_df_filtered <- test_df_filtered %>%
  filter(!amongAggDuplicate | is.na(AggDuplicateGroupID)) %>% # Keep non-duplicates
  bind_rows( # Add in the lowest-ranked duplicates
    test_df_filtered %>%
      filter(amongAggDuplicate) %>%
      group_by(AggDuplicateGroupID) %>%
      filter(AggDuplicateRank == min(AggDuplicateRank, na.rm = TRUE)) %>%
      slice(1) %>%  # if tied for lowest rank, just take the first
      ungroup()
  )

test_df_filtered <- test_df_filtered %>%
  filter(!specimenDuplicate | is.na(specimenDuplicateGroupID)) %>% # Keep non-duplicates
  bind_rows( # Add in the lowest-ranked duplicates
    test_df_filtered %>%
      filter(specimenDuplicate) %>%
      group_by(specimenDuplicateGroupID) %>%
      filter(specimenDuplicateRank == min(specimenDuplicateRank, na.rm = TRUE)) %>%
      slice(1) %>%  # if tied for lowest rank, just take the first
      ungroup()
  )
# additionally, depending on the resolution we want to filter out coordinate uncertainty that we can live with
test_df_50 <- test_df_filtered %>% 
  filter(coordinateUncertaintyInMeters <= 25000) # we're just going with half as a rule
test_df_25 <- test_df_filtered %>% 
  filter(coordinateUncertaintyInMeters <= 12500)
combined_extent_50_grid <- st_make_grid(combined_extent, 
                                        crs = crs(r_50_to_50), 
                                        cellsize = c(50000, 50000))
combined_extent_50_grid_sf <- st_sf(geometry = combined_extent_50_grid)
combined_extent_50_grid_sf <- combined_extent_50_grid_sf %>% 
  mutate(id = 1:n())
combined_extent_25_grid <- st_make_grid(combined_extent, 
                                        crs = crs(r_50_to_25), 
                                        cellsize = c(25000, 25000))
combined_extent_25_grid_sf <- st_sf(geometry = combined_extent_25_grid)
combined_extent_25_grid_sf <- combined_extent_25_grid_sf %>% 
  mutate(id = 1:n())
# convert these data to pts in raster format
#pts <- terra::vect(as.data.frame(test_df_50), geom = c("roundedLongitude", "roundedLatitude"), crs = "EPSG:4326")
#pts <- terra::project(pts, crs(template_50))
# rasterize using unique index per raster cell
#test_df_50$cell_id <- terra::cellFromXY(template_50, terra::geom(pts)[, c("x", "y")])
occ_pts <- st_as_sf(test_df_50, coords = c("roundedLongitude", "roundedLatitude"), crs = 4326)
occ_pts <- st_transform(occ_pts, crs = st_crs(combined_extent_50_grid))
joined <- st_join(occ_pts, combined_extent_50_grid_sf, join = st_intersects)
joined <- joined %>% filter(!is.na(id))
# get p/a if the grid contains a centroid 
presence <- unique(joined$id)
sp_grid <- combined_extent_50_grid_sf # make a copy for this particular sp
sp_grid <- sp_grid %>% # label presence absence 
  dplyr::mutate(presence = ifelse(id %in% presence, 1, 0)) %>% 
  dplyr::mutate(presence = as.factor(presence))
# create an empty raster template, contains extent of original grid at 10x10km, assuure crs also matches
r_template <- rast(ext(sp_grid), resolution = 50000, crs = st_crs(sp_grid)$wkt)
# maps presence to raster cells, choose fun = "max" to take 1 if available
r_presence <- rasterize(vect(sp_grid), r_template, field = "presence", fun = "max")
r_occs_50 <- r_presence # rename for variable clarity
# check how this changed the data:
r_points <- as.points(r_occs_50)
r_points_sf <- st_as_sf(r_points)
bbox <- st_bbox(filter(r_points_sf, presence == 1))
occ_pts_proj <- st_transform(occ_pts, crs = crs(combined_extent_50_grid_sf))
# ggplot
ggplot() +
  geom_sf(data = occ_pts_proj, color = "red", size = 2) +
  geom_sf(data = r_points_sf, aes(color = as.factor(presence)), size = 1) +
  geom_sf(data = combined_extent_50_grid_sf, fill = NA, color = "grey") +
  scale_color_manual(values = c("white", "darkgreen")) +
  coord_sf(xlim = c(bbox["xmin"] - 100000, bbox["xmax"] + 100000), # add slight buffers so its viewable
           ylim = c(bbox["ymin"] - 100000, bbox["ymax"] + 100000)) +
  ggtitle(paste0("Presence/Absence raster vs. points for ")) +
  theme_minimal()

# rename for var clarity
r_point <- r_presence 
ext(r_point) == ext(combined_extent) # check to make sure this is TRUE
crs(r_point) == crs(r_10) # Check to see if TRUE
crs(r_point) == crs(r_50_proj) # Check to see if TRUE

# now resample
r_point_to_50 <- resample(r_point, template_50, method = "max")

# check results
r_points <- as.points(r_point_to_50)
r_points_sf <- st_as_sf(r_points)
bbox <- st_bbox(filter(r_points_sf, presence == 1))
crs(r_points_sf) == crs(r_occs_50)

ggplot() +
  geom_sf(data = occ_pts_proj , color = "red", size = 1) +
  geom_sf(data = r_points_sf, aes(color = as.factor(presence)), size = 1) +
  geom_sf(data = combined_extent_50_grid, fill = NA, color = "grey") +
  scale_color_manual(values = c("white", "darkgreen")) +
  coord_sf(xlim = c(bbox["xmin"] - 100000, bbox["xmax"] + 100000), # add slight buffers so its viewable
           ylim = c(bbox["ymin"] - 100000, bbox["ymax"] + 100000)) +
  theme_minimal()

# Operating as intending, check r_point_to_25
occ_pts <- st_as_sf(test_df_25, coords = c("roundedLongitude", "roundedLatitude"), crs = 4326)
occ_pts <- st_transform(occ_pts, crs = st_crs(combined_extent_25_grid))
joined <- st_join(occ_pts, combined_extent_25_grid_sf, join = st_intersects)
joined <- joined %>% filter(!is.na(id))
# get p/a if the grid contains a centroid 
presence <- unique(joined$id)
sp_grid <- combined_extent_25_grid_sf # make a copy for this particular sp
sp_grid <- sp_grid %>% # label presence absence 
  dplyr::mutate(presence = ifelse(id %in% presence, 1, 0)) %>% 
  dplyr::mutate(presence = as.factor(presence))
# create an empty raster template, contains extent of original grid at 10x10km, assuure crs also matches
r_template <- rast(ext(sp_grid), resolution = 25000, crs = st_crs(sp_grid)$wkt)
# maps presence to raster cells, choose fun = "max" to take 1 if available
r_presence <- rasterize(vect(sp_grid), r_template, field = "presence", fun = "max")
r_occs_25 <- r_presence # rename for variable clarity
# check how this changed the data:
r_points <- as.points(r_occs_25)
r_points_sf <- st_as_sf(r_points)
bbox <- st_bbox(filter(r_points_sf, presence == 1))
occ_pts_proj <- st_transform(occ_pts, crs = crs(combined_extent_25_grid_sf))
# ggplot
ggplot() +
  geom_sf(data = occ_pts_proj, color = "red", size = 2) +
  geom_sf(data = r_points_sf, aes(color = as.factor(presence)), size = 1) +
  geom_sf(data = combined_extent_25_grid_sf, fill = NA, color = "grey") +
  scale_color_manual(values = c("white", "darkgreen")) +
  coord_sf(xlim = c(bbox["xmin"] - 100000, bbox["xmax"] + 100000), # add slight buffers so its viewable
           ylim = c(bbox["ymin"] - 100000, bbox["ymax"] + 100000)) +
  ggtitle(paste0("Presence/Absence raster vs. points for ")) +
  theme_minimal()

# rename for var clarity
r_point <- r_presence 
ext(r_point) == ext(combined_extent) # check to make sure this is TRUE
crs(r_point) == crs(r_10) # Check to see if TRUE
crs(r_point) == crs(r_50_proj) # Check to see if TRUE

# now resample
r_point_to_25 <- resample(r_point, template_25, method = "max")

# check results
r_points <- as.points(r_point_to_25)
r_points_sf <- st_as_sf(r_points)
bbox <- st_bbox(filter(r_points_sf, presence == 1))
crs(r_points_sf) == crs(occ_pts_proj)

ggplot() +
  geom_sf(data = occ_pts_proj , color = "red", size = 1) +
  geom_sf(data = r_points_sf, aes(color = as.factor(presence)), size = 1) +
  geom_sf(data = combined_extent_25_grid, fill = NA, color = "grey") +
  scale_color_manual(values = c("white", "darkgreen")) +
  coord_sf(xlim = c(bbox["xmin"] - 100000, bbox["xmax"] + 100000), # add slight buffers so its viewable
           ylim = c(bbox["ymin"] - 100000, bbox["ymax"] + 100000)) +
  theme_minimal()
