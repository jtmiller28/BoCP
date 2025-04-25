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

# Ensure geometries are valid
na_regions_valid <- st_make_valid(na_regions)

# Use st_point_on_surface() for safer label placement
na_centroids <- st_point_on_surface(na_regions_valid)

# Extract coordinates
na_centroids <- na_centroids %>%
  mutate(
    lon = st_coordinates(.)[,1],
    lat = st_coordinates(.)[,2]
  )
# Plot with text labels
ggplot() +
  geom_sf(data = bot_regions) +
  geom_sf(data = na_regions, fill = "red", color = NA) +
  geom_text(data = na_centroids,
            aes(x = lon, y = lat, label = LEVEL3_COD),
            size = 1.5, color = "black") +
  ggtitle("North American Extent for BoCP Implementation")

# make plot thats more visually easy to see as an extent plot 
ggplot() + 
  geom_sf(data = bot_regions) + 
  geom_sf(data = na_regions, fill = "darkred", color = "black") + 
  ggtitle("USA and Canada for Endemic Pixel Project")


# Next we want to build the study region into a equal area grid: 
## first lets reproj into the projection native to what qgis said
crs_102008 <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 
               +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
#na_regions_aea <- st_transform(na_regions, crs = crs_aea)
na_regions_aea <- st_transform(na_regions, crs = crs_102008)
## plot to check and see if things look about right
ggplot() + 
  geom_sf(data = na_regions_aea, fill = "darkred", color = "black") + 
  ggtitle("USA and Canada for Endemic Pixel Project")

## create 50x50km grid 
grid_50km <- st_make_grid(na_regions_aea, 
                          cellsize = c(50000, 50000), 
                          square = TRUE) %>% 
  st_as_sf() %>% # convert to sf 
  mutate(cell_id = row_number()) # assign uniq id to each cell

# clip to shape 
grid_50km_clipped <- st_intersection(grid_50km, na_regions_aea)

## Create an alternative grid
x_min <- -5105000.0000000000000000
y_min <- -2905000.0000000000000000
x_max<- 3045000.0001220712438226
y_max <- 4645000.0001220703125000
ns_bbox <- st_bbox(c(xmin = x_min, ymin = y_min, xmax = x_max, ymax = y_max), crs = crs_102008)

ns_grid <- st_make_grid(ns_bbox, cellsize = c(50000, 50000), what = "polygons", square = TRUE)
## plot to check and see how that looks
ggplot() + 
  geom_sf(data = na_regions_aea, fill = "darkred", color = "black") + 
  geom_sf(data = ns_grid, fill = "black", alpha = 0.7) + 
  ggtitle("USA and Canada with 50x50km grid \n for Endemic Pixel Project")

## Load in Israels shapefile nature serve data
nature_serve_shape <- read_sf("/blue/guralnick/millerjared/BoCP/data/raw/grid-data/na-rare-plant-grid-data/NA spatial phylogenetics/tns_points.shp")
# missing prj file, so define coordinate projection provided by Israel via emial: Albers Equal Area ESPG:3310
nature_serve_shape <- st_set_crs(nature_serve_shape, crs_102008)
#nature_serve_shape <- st_transform(nature_serve_shape, crs = st_crs(na_regions_aea))
# take a subset of these data to better understand its structure. 
nature_serve_shapet <- nature_serve_shape %>% 
  group_by(SNAME) %>% 
  mutate(uniq_loc = length(unique(geometry))) %>% 
  select(SNAME, uniq_loc) %>% 
  filter(!is.na(SNAME))
#test_ns <- nature_serve_shape %>% group_by(SNAME) %>% distinct() %>% summarize(n = n())
nature_serve_test <- nature_serve_shape %>%  # create a test shape
  filter(SNAME == "Rhexia salicifolia") %>% 
  st_transform(crs = st_crs(na_regions_aea))
# define the extent of this test data set
bbox <- st_bbox(nature_serve_test)
ggplot() +  
  geom_sf(na_regions_aea, mapping = aes()) + 
  geom_sf(ns_grid, mapping = aes(), alpha = 0.2) +
  geom_sf(nature_serve_test, mapping = aes()) +
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
           ylim = c(bbox["ymin"], bbox["ymax"])) +
  ggtitle("Nature Serve 50x50km Centroid Presence Map Aeschynomene virginica")

# run as a for loop 
test_names <- c("Tumamoca macdougalii", "Heuchera eastwoodiae", "Lachnanthes caroliniana", 
                "Agalinis filicaulis", "Carex dasycarpa", "Rhexia salicifolia", "Parnassia asarifolia",
                "Viburnum lentago", "Echinacea pallida", "Trillium sessile")

for(i in 1:length(test_names)){
  nature_serve_test <- nature_serve_shape %>%  # create a test shape
    filter(SNAME == test_names[i]) %>% 
    st_transform(crs = st_crs(na_regions_aea))
  bbox <- st_bbox(nature_serve_test)
  ggplot() +  
    geom_sf(na_regions_aea, mapping = aes()) + 
    geom_sf(grid_50km, mapping = aes(), alpha = 0.2) +
    geom_sf(nature_serve_test, mapping = aes()) +
    coord_sf(xlim = c(bbox["xmin"] - 10000, bbox["xmax"] + 10000),
             ylim = c(bbox["ymin"] - 10000, bbox["ymax"] + 10000)) +
    ggtitle(paste0("Nature Serve 50x50km presence data for ", test_names[i]))
  ggsave(paste0("/blue/guralnick/millerjared/BoCP/outputs/example-grid-extractions/", test_names[i], "-NS-sample.png"), width = 10, height = 10)
}


## Do the same with the canada data
canada1_shape <- read_sf("/blue/guralnick/millerjared/BoCP/data/raw/grid-data/na-rare-plant-grid-data/NA spatial phylogenetics/canada1.shp")
test_ca <- canada1_shape %>% group_by(NNAME)%>% summarize(n = n()) %>% ungroup()
test_ca <- test_ca %>% group_by(NNAME) %>% 
  mutate(uniq_geoms = st_coordinates(geometry))

canada1_shape_t <- canada1_shape %>% 
  group_by(NNAME) %>% 
  mutate(uniq_loc = length(unique(geometry)))

canada1_shape_test <- canada1_shape %>% 
  filter(NNAME == "Juncus pylaei")

# transform na to canada's shapefile projection
na_regions_can <- st_transform(na_regions, crs = st_crs(canada1_shape))

## create 50x50km grid canada style
grid_50km <- st_make_grid(na_regions_can, 
                          cellsize = c(10000, 10000), 
                          square = TRUE) %>% 
  st_as_sf() %>% # convert to sf 
  mutate(cell_id = row_number()) # assign uniq id to each cell

# clip to shape 
grid_50km_clipped <- st_intersection(grid_50km, na_regions_can)

bbox <- st_bbox(canada1_shape_test)

ggplot() +  
  geom_sf(na_regions_can, mapping = aes()) + 
  geom_sf(grid_50km, mapping = aes(), alpha = 0.2) +
  geom_sf(canada1_shape_test, mapping = aes()) +
  coord_sf(xlim = c(bbox["xmin"] - 10000, bbox["xmax"] + 10000),
           ylim = c(bbox["ymin"] - 10000, bbox["ymax"] + 10000)) +
  ggtitle("Canada 50x50km presence data for Carex comosa")

### Create a loop to validate for rob 
test_names <- c("Solidago brendae", "Lycopus laurentianus", "Elymus canadensis", 
                "Astragalus eucosmus", "Elatine americana", "Astragalus eucosmus", 
                "Stachys hispida", "Galearis rotundifolia", "Juncus greenei",
                "Polygonum achoreum")

for(i in 1:length(test_names)){
canada1_shape_test <- canada1_shape %>% 
    filter(NNAME == test_names[i])
bbox <- st_bbox(canada1_shape_test)
ggplot() +  
  geom_sf(na_regions_can, mapping = aes()) + 
  geom_sf(grid_50km, mapping = aes(), alpha = 0.2) +
  geom_sf(canada1_shape_test, mapping = aes()) +
  coord_sf(xlim = c(bbox["xmin"] - 10000, bbox["xmax"] + 10000),
           ylim = c(bbox["ymin"] - 10000, bbox["ymax"] + 10000)) +
  ggtitle(paste0("Canada 10x10km presence data for ", test_names[i]))
ggsave(paste0("/blue/guralnick/millerjared/BoCP/outputs/example-grid-extractions/", test_names[i], "-CAN-sample.png"), width = 10, height = 10)
}


### Fxns 
# Load presence points shapefile
pts <- read_sf("/blue/guralnick/millerjared/BoCP/data/raw/grid-data/na-rare-plant-grid-data/NA spatial phylogenetics/tns_points.shp")
# missing prj file, so define coordinate projection provided by Israel via emial: Albers Equal Area ESPG:3310
pts <- st_set_crs(pts, crs_102008)

pts <- pts %>% filter(SNAME == "Rhexia salicifolia")

# Load and transform the North America shapefile
na <- na_regions %>% st_transform(crs_102008)

# create extent matching the intial grid (found via loading file in QGIS)
x_min <- -5026845.9361000061035156
x_max <- -1713353.4216001033782959
y_min <- 3023154.0638999938964844
y_max <- 4436646.5784001350402832
ns_bbox <- st_bbox(c(xmin = x_min, ymin = y_min, xmax = x_max, ymax = y_max), crs = crs_102008)
# Define grid cell size
cell_size <- 50000

# Compute range of cell centers (x and y)
x_seq <- seq(
  from = floor(na_bbox$xmin / cell_size) * cell_size + cell_size / 2,
  to   = ceiling(na_bbox$xmax / cell_size) * cell_size - cell_size / 2,
  by   = cell_size
)

y_seq <- seq(
  from = floor(na_bbox$ymin / cell_size) * cell_size + cell_size / 2,
  to   = ceiling(na_bbox$ymax / cell_size) * cell_size - cell_size / 2,
  by   = cell_size
)

# Create all combinations of center points
grid_centers <- expand.grid(x = x_seq, y = y_seq)

# Convert to sf points (same logic as original centering)
grid_points <- st_as_sf(grid_centers, coords = c("x", "y"), crs = crs_102008)

# Your original grid-cell builder function
make_grid_cell <- function(point, cell_size = 50000) {
  x <- st_coordinates(point)[1]
  y <- st_coordinates(point)[2]
  xmin <- x - cell_size / 2
  ymin <- y - cell_size / 2
  xmax <- x + cell_size / 2
  ymax <- y + cell_size / 2
  st_polygon(list(matrix(c(
    xmin, ymin,
    xmin, ymax,
    xmax, ymax,
    xmax, ymin,
    xmin, ymin
  ), ncol = 2, byrow = TRUE))) %>%
    st_sfc(crs = st_crs(point))
}

# Build grid cells using the above centroids
grid_cells <- do.call(rbind, lapply(1:nrow(grid_points), function(i) {
  st_sf(geometry = make_grid_cell(grid_points[i, ]))
}))

bbox <- st_bbox(pts)
# Plot
ggplot() +
  geom_sf(data = na, mapping = aes()) +
  geom_sf(data = grid_cells, fill = NA, color = "grey40") +
  geom_sf(data = pts, color = "red", size = 1) +
  #coord_sf(crs = st_crs(crs_102008)) +
  coord_sf(xlim = c(bbox["xmin"] - 10000, bbox["xmax"] + 10000),
           ylim = c(bbox["ymin"] - 10000, bbox["ymax"] + 10000)) +
  theme_minimal() +
  labs(title = "Presence Points and Reconstructed 50x50 km Grid",
       caption = "Red points = presence records, grey squares = grid cells")


## 