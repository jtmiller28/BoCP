

## Harmonize Nature Serve's Taxonomy 
ns_aligned_names <- fread("/blue/guralnick/millerjared/BoCP/data/processed/nature-serve-names-aligned.csv")
ns_centroids_a <- merge(nature_serve_centroids, ns_aligned_names, by.x = "SNAME", by.y = "user_supplied_name")

## Harmonize Canada Data Taxonomy 
canada1_centroids <- read_sf("/blue/guralnick/millerjared/BoCP/data/raw/grid-data/na-rare-plant-grid-data/NA spatial phylogenetics/canada1.shp")
canada2_centroids <- read_sf("/blue/guralnick/millerjared/BoCP/data/raw/grid-data/na-rare-plant-grid-data/NA spatial phylogenetics/canada2.shp")
canada1_centroids <- canada1_centroids %>% 
  select(NNAME) %>% 
  rename(SNAME = NNAME)
canada2_centroids <- canada2_centroids %>% 
  select(GNAME) %>% 
  rename(SNAME = GNAME)
canada_centroids <- rbind(canada1_centroids, canada2_centroids)
can_aligned_names <- fread("/blue/guralnick/millerjared/BoCP/data/processed/canada-names-aligned.csv")
can_centroids_a <- merge(canada_centroids, can_aligned_names, by.x = "SNAME", by.y = "user_supplied_name" )
canada_grid <- read_sf("./data/raw/grid-data/na-rare-plant-grid-data/NA spatial phylogenetics/NAm_cell_template_10km/NAm_cell_template_10km.shp")
canada_grid <- canada_grid %>% 
  select(-Id, -GRID_ID) %>% 
  mutate(id = 1:n())

can_test <- filter(can_centroids_a, alignedParentName == "Abies balsamea")
joined <- st_join(can_test, canada_grid, join = st_within)

# get p/a if the grid contains a centroid 
presence <- unique(joined$id)
sp_grid <- canada_grid # make a copy for this particular sp
sp_grid <- sp_grid %>% # label presence absence 
  dplyr::mutate(presence = ifelse(id %in% presence, 1, 0))


### Try combining
# Extract a name from both datasets
ns_test <- ns_centroids_a %>%  # 50x50
  filter(alignedParentName == "Abies balsamea") 

can_test <- can_centroids_a %>% #10x10
  filter(alignedParentName == "Abies balsamea")

# Standardize and rasterize
# preform a spatial join to determine presence within a cell
joined_ns <- st_join(ns_test, ns_grid_sf, join = st_intersects)
joined_can <- st_join(can_test, canada_grid, join = st_intersects)
# get p/a if the grid contains a centroid 
presence_ns <- unique(joined_ns$id)
presence_can <- unique(joined_can$id)
sp_grid_ns <- ns_grid_sf # make a copy for this particular sp
sp_grid_can <- canada_grid # make a copy for this particular sp
sp_grid_ns <- sp_grid_ns %>% # label presence absence 
  dplyr::mutate(presence = ifelse(id %in% presence_ns, 1, 0))
sp_grid_can <- sp_grid_can %>% # label presence absence 
  dplyr::mutate(presence = ifelse(id %in% presence_can, 1, 0))
# create an empty raster template, contains extent of original grid at 50x50km, assuure crs also matches
r_template_ns <- rast(ext(sp_grid_ns), resolution = 50000, crs = st_crs(sp_grid_ns)$wkt)
r_template_can <- rast(ext(sp_grid_can), resolution = 10000, crs = st_crs(sp_grid_can)$wkt)
# maps presence to raster cells, choose fun = "max" to take 1 if available
r_presence_ns <- rasterize(vect(sp_grid_ns), r_template_ns, field = "presence", fun = "max")
r_presence_can <- rasterize(vect(sp_grid_can), r_template_can, field = "presence", fun = "max")

r_50_ns <- r_presence_ns
r_10_can <- r_presence_can
crs(r_50_ns ) == crs(r_10_can)
# false, so reproj
r_50_ns_proj <- project(r_50_ns, r_10_can, method = "near") # project r_50 to preserve finer scale
crs(r_50_ns_proj) == crs(r_10_can)
# true, now proceed to make a common extent
combined_extent <- ext(r_50_ns_proj) + ext(r_10_can)
template_50 <- rast(ext = combined_extent, resolution = 50000, crs = crs(r_50_ns_proj))
template_25 <- rast(ext = combined_extent, resolution = 25000, crs = crs(r_10_can))
# resample each to the other grid...
r_10_to_50 <- resample(r_10_can, template_50, method = "near")
r_10_to_25 <- resample(r_10_can, template_25, method = "near")
r_50_to_50 <- resample(r_50_ns_proj, template_50, method = "near") # prob redudant
r_50_to_25 <- resample(r_50_ns_proj, template_25, method = "near")
# create combined P/A grids (summarize the layers into one layer using app() with function max for 1 if present)
combined_50 <- app(c(r_50_to_50, r_10_to_50), fun = max, na.rm = TRUE)
combined_25 <- app(c(r_50_to_25, r_10_to_25), fun = max, na.rm = TRUE)

# ALT 
combined_50 <- r_50_to_50
values(combined_50) <- pmax(values(r_50_to_50), values(r_10_to_50), na.rm = TRUE)

combined_25 <- r_50_to_25
values(combined_25) <- pmax(values(r_50_to_25), values(r_10_to_25), na.rm = TRUE)


# plot to ensure things worked

# Convert rasters to data frames for ggplot
r_50_df <- as.data.frame(combined_50, xy = TRUE, na.rm = FALSE)
r_25_df <- as.data.frame(combined_25, xy = TRUE, na.rm = FALSE)

# Rename raster value column for clarity
colnames(r_50_df)[3] <- "presence"
colnames(r_25_df)[3] <- "presence"

# Transform centroids and grids to same CRS if needed
ns_test <- st_transform(ns_test, crs = st_crs(sp_grid_ns))
can_test <- st_transform(can_test, crs = st_crs(sp_grid_can))

# Plot for 50x50km
bbox <- st_bbox(can_test)
p_50 <- ggplot() +
  geom_raster(data = r_50_df, aes(x = x, y = y, fill = factor(presence))) +
  geom_sf(data = st_geometry(sp_grid_ns), fill = NA, color = "black", size = 0.2) +
  geom_sf(data = ns_test, aes(), color = "blue", size = 1.2, alpha = 0.7) +
  geom_sf(data = can_test, aes(), color = "red", shape = 17, size = 1.2, alpha = 0.7) +
  coord_sf(xlim = c(bbox["xmin"] - 10000, bbox["xmax"] + 10000), # add slight buffers so its viewable
           ylim = c(bbox["ymin"] - 10000, bbox["ymax"] + 10000)) +
  scale_fill_manual(values = c("0" = "white", "1" = "green"), na.value = "grey90", name = "Presence") +
  ggtitle("50x50km Grid Overlay with Combined Presence") +
  theme_minimal()

# Plot for 25x25km
p_25 <- ggplot() +
  geom_raster(data = r_25_df, aes(x = x, y = y, fill = factor(presence))) +
  geom_sf(data = st_geometry(sp_grid_can), fill = NA, color = "black", size = 0.2) +
  geom_sf(data = ns_test, aes(), color = "blue", size = 1.2, alpha = 0.7) +
  geom_sf(data = can_test, aes(), color = "red", shape = 17, size = 1.2, alpha = 0.7) +
  scale_fill_manual(values = c("0" = "white", "1" = "green"), na.value = "grey90", name = "Presence") +
  ggtitle("25x25km Grid Overlay with Combined Presence") +
  theme_minimal() +
  coord_sf()


## Try again...
# set up our raster cells
r_50_ns <- r_presence_ns
r_10_can <- r_presence_can
r_50_ns_proj <- project(r_50_ns, r_10_can, method = "near")
# fix up extents
r_50_ns_proj <- extend(r_50_ns_proj, ext(r_10_can))
r_10_can <- extend(r_10_can, ext(r_50_ns_proj))
ext(r_10_can)
ext(r_50_ns_proj)
# create a template raster to resample to
template_50 <- rast(ext = ext(r_50_ns_proj), resolution = 50000, crs = crs(r_50_ns_proj))
template_25 <- rast(ext = ext(r_50_ns_proj), resolution = 25000, crs = crs(r_50_ns_proj))
r_50_25km <- terra::resample(r_50_ns_proj, template_25)
r_50_50km <- resample(r_50_ns_proj, template_50)
r_10_25km <- resample(r_10_can, template_25)
r_10_50km <- resample(r_10_can, template_50)

# Merge the two rasters at each resolution
r_25km <- (r_50_25km + r_10_25km) > 0
r_50km <- (r_50_50km + r_10_50km) > 0


# Load the ggplot2 package
library(ggplot2)

# Convert the rasters to data frames
r1_df <- as.data.frame(r_50_ns_proj, xy=TRUE)
r2_df <- as.data.frame(r_10_can, xy=TRUE)
r_25km_df <- as.data.frame(r_25km, xy=TRUE)
r_50km_df <- as.data.frame(r_50km, xy=TRUE)

# Rename the column to "presence"
colnames(r1_df)[3] <- "presence"
colnames(r2_df)[3] <- "presence"
colnames(r_25km_df)[3] <- "presence"
colnames(r_50km_df)[3] <- "presence"



r1_df_present <- filter(r1_df, presence == 1)
buffer <- 100000
zoom_df <- subset(r1_df, x >= min(r1_df_present$x) + buffer & x <= max(r1_df_present$x)+ buffer & y >= min(r1_df_present$y) + buffer & y <= max(r1_df_present$y) + buffer)

ggplot(zoom_df, aes(x = x, y = y, fill = factor(presence))) + 
  geom_raster() + 
  scale_fill_manual(values = c("steelblue", "goldenrod"), labels = c("Absence", "Presence")) + 
  theme_void() + 
  labs(title = "Zoomed Raster 1") + 
  coord_fixed(ratio = 1) + 
  theme(legend.position = "bottom")





# Plot the rasters
ggplot(r1_df, aes(x=x, y=y, fill=factor(presence))) + 
  geom_raster() + 
  scale_fill_manual(values=c("blue", "red"), labels=c("Absence", "Presence")) + 
  theme_void() + 
  labs(title="Original Raster 1") + 
  coord_fixed(ratio=1) + 
  theme(legend.position="bottom")

ggplot(r2_df, aes(x=x, y=y, fill=layer)) + 
  geom_raster() + 
  scale_fill_manual(values=c("blue", "red")) + 
  theme_void() + 
  labs(title="Original Raster 2")

ggplot(r_25km_df, aes(x=x, y=y, fill=layer)) + 
  geom_raster() + 
  scale_fill_manual(values=c("blue", "red")) + 
  theme_void() + 
  labs(title="Merged Raster 25km")

ggplot(r_50km_df, aes(x=x, y=y, fill=layer)) + 
  geom_raster() + 
  scale_fill_manual(values=c("blue", "red")) + 
  theme_void() + 
  labs(title="Merged Raster 50km")














# Show plots
print(p_50)
print(p_25)




# export as raster and shapefiles to give flexibility down the line
## raster
# create an empty raster template, contains extent of original grid at 50x50km, assuure crs also matches
r_template <- rast(ext(sp_grid), resolution = 50000, crs = st_crs(sp_grid)$wkt)
# maps presence to raster cells, choose fun = "max" to take 1 if available
r_presence <- rasterize(vect(sp_grid), r_template, field = "presence", fun = "max")



# use aligned parentName to match our BoCP harmonized taxonomy from now on. 
# Rebuild grid for clarity 
# create a 50x50km grid that matches the extent of the NS data, extent was retrieved manually by reviewing the same file mentioned above in QGIS
x_min <- -5105000.0000000000000000
y_min <- -2905000.0000000000000000
x_max<- 3045000.0001220712438226
y_max <- 4645000.0001220703125000
ns_bbox <- st_bbox(c(xmin = x_min, ymin = y_min, xmax = x_max, ymax = y_max), crs = crs_102008)
ns_grid <- st_make_grid(ns_bbox, cellsize = c(50000, 50000), what = "polygons", square = TRUE)
# reformat our grid as an sf object
ns_grid_sf <- st_sf(id = 1:length(ns_grid), geometry = ns_grid)
# create a presence absence grid for our three datasets
