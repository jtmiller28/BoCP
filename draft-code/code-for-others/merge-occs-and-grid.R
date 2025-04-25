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

for(i in 1:length(test_names)){
  nature_serve_test <- nature_serve_centroids %>%  # extract species level centroid pts
    filter(SNAME == test_names[i]) %>% 
    st_transform(crs = st_crs(crs_102008))
  bbox <- st_bbox(nature_serve_test) # create a bbox to view their extent
  ggplot() +  
    geom_sf(na_regions_aea, mapping = aes()) + 
    geom_sf(ns_grid, mapping = aes(), alpha = 0.2) +
    geom_sf(nature_serve_test, mapping = aes()) +
    coord_sf(xlim = c(bbox["xmin"] - 10000, bbox["xmax"] + 10000), # add slight buffers so its viewable
             ylim = c(bbox["ymin"] - 10000, bbox["ymax"] + 10000)) +
    ggtitle(paste0("Nature Serve 50x50km presence data for ", test_names[i]))
  ggsave(paste0("/blue/guralnick/millerjared/BoCP/outputs/example-grid-extractions/", test_names[i], "-NS-sample-zoom.png"), width = 10, height = 10)
  
  ggplot() +  
    geom_sf(na_regions_aea, mapping = aes()) + 
    geom_sf(ns_grid, mapping = aes(), alpha = 0.1) +
    geom_sf(nature_serve_test, mapping = aes()) +
     coord_sf(xlim = c(bbox["xmin"] - 1000000, bbox["xmax"] + 1000000), # add slight buffers so its viewable
              ylim = c(bbox["ymin"] - 1000000, bbox["ymax"] + 1000000)) +
    ggtitle(paste0("Nature Serve 50x50km presence data for ", test_names[i]))
  ggsave(paste0("/blue/guralnick/millerjared/BoCP/outputs/example-grid-extractions/", test_names[i], "-NS-sample-.png"), width = 10, height = 10)
}

## Next lets use the same concept to delimit data for the 10x10km grid cells for canada. 
# import canada grid Israel has provided 
canada_grid <- read_sf("./data/raw/grid-data/na-rare-plant-grid-data/NA spatial phylogenetics/NAm_cell_template_10km/NAm_cell_template_10km.shp")

# projection should be identical, so lets just reset up the bounding box 
# x_min <- 1900000.0000000000000000
# y_min <- 760000.0000000000000000
# x_max<- 2670000.0000000000000000
# y_max <- 1330000.0000000000000000
# ca1_bbox <- st_bbox(c(xmin = x_min, ymin = y_min, xmax = x_max, ymax = y_max), crs = crs_102008)
# ca1_grid <- st_make_grid(ca1_bbox, cellsize = c(10000, 10000), what = "polygons", square = TRUE)
# read in the canada 10x10km centroid data
canada1_centroids <- read_sf("/blue/guralnick/millerjared/BoCP/data/raw/grid-data/na-rare-plant-grid-data/NA spatial phylogenetics/canada1.shp")
canada2_centroids <- read_sf("/blue/guralnick/millerjared/BoCP/data/raw/grid-data/na-rare-plant-grid-data/NA spatial phylogenetics/canada2.shp")
# create a test set of names that are known prior 
test_names <- c("Solidago brendae", "Lycopus laurentianus", "Elymus canadensis", 
                "Astragalus eucosmus", "Elatine americana", "Astragalus eucosmus", 
                "Stachys hispida", "Galearis rotundifolia", "Juncus greenei",
                "Polygonum achoreum")

for(i in 1:length(test_names)){
  can1_test <- canada2_centroids %>%  # extract species level centroid pts
    filter(GNAME == test_names[i]) %>% 
    st_transform(crs = st_crs(crs_102008))
  bbox <- st_bbox(can1_test) # create a bbox to view their extent
  ggplot() +  
    geom_sf(na_regions_aea, mapping = aes()) + 
    geom_sf(canada_grid, mapping = aes(), alpha = 0.2) +
    geom_sf(can1_test, mapping = aes()) +
    coord_sf(xlim = c(bbox["xmin"] - 10000, bbox["xmax"] + 10000), # add slight buffers so its viewable
             ylim = c(bbox["ymin"] - 10000, bbox["ymax"] + 10000)) +
    ggtitle(paste0("Canada1 10x10km presence data for ", test_names[i]))
  ggsave(paste0("/blue/guralnick/millerjared/BoCP/outputs/example-grid-extractions/", test_names[i], "-CAN1-sample-zoom.png"), width = 10, height = 10)
  
  ggplot() +  
    geom_sf(na_regions_aea, mapping = aes()) + 
    geom_sf(ca1_grid, mapping = aes(), alpha = 0.1) +
    geom_sf(can1_test, mapping = aes()) +
    coord_sf(xlim = c(bbox["xmin"] - 1000000, bbox["xmax"] + 1000000), # add slight buffers so its viewable
             ylim = c(bbox["ymin"] - 1000000, bbox["ymax"] + 1000000)) +
    ggtitle(paste0("Canda1 10x10km presence data for ", test_names[i]))
  ggsave(paste0("/blue/guralnick/millerjared/BoCP/outputs/example-grid-extractions/", test_names[i], "-CAN1-sample-.png"), width = 10, height = 10)
}

