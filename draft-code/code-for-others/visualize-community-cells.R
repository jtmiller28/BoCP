### Title: Visualize Community Cells
### Author: JT Miller
### Date: 05/02/2025

# load libraries
library(data.table)
library(tidyverse)
library(terra)
library(sf)
library(patchwork)

# read in the community P/A 
community_pixels_IN <- fread("/blue/guralnick/millerjared/BoCP/data/processed/grid-project/worldwide/raster-species-PA-visuals/community-richness-50x50/community-sp-PA-makeup-intro-native.csv")

# set up shapefile (delimited to the same pixels as found in clean-rasters-to-geo-scope.R)
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
# set up mollweide proj centered on NA
moll_proj_center <- "+proj=moll +lon_0=-100 +datum=WGS84 +units=m +no_defs"
# reproj our na regions to this projection 
na_regions <- st_transform(na_regions, crs = crs(moll_proj_center))

# covert community pixels to summary form 
community_pixels_summary <- community_pixels_IN %>% 
  mutate(xy = paste(x, y)) %>% 
  group_by(xy) %>% 
  reframe(x,y, n = sum(presence)) %>% 
  ungroup() %>% 
  select(-xy) 




# plot
intro_native_50 <- ggplot() + 
  geom_tile(community_pixels_summary, mapping = aes(x = x, y = y, fill = n)) + 
  scale_fill_viridis_c(na.value = "transparent") + 
  geom_sf(na_regions, fill = NA, color = "black", size = 0.3, mapping = aes()) +
  theme_minimal() + 
  ggtitle("Introduced & Native \n 50x50km Presence Cells for Plant Community",) + 
  xlab("Longitude") + 
  ylab("Latitude") +
  theme(plot.title = element_text(hjust = 0.5))
intro_native_50
ggsave("/blue/guralnick/millerjared/BoCP/data/processed/grid-project/worldwide/raster-species-PA-visuals/community-richness-50x50/community-richness-intro-native-50km.png", width = 10, height = 10)

log_intro_native_50 <- ggplot() + 
  geom_tile(community_pixels_summary, mapping = aes(x = x, y = y, fill = log(n+1))) + 
  scale_fill_viridis_c(na.value = "transparent") + 
  geom_sf(na_regions, fill = NA, color = "black", size = 0.3, mapping = aes()) +
  theme_minimal() + 
  ggtitle("Introduced & Native \n 50x50km Log Transformed Presence Cells \n for Plant Community",) + 
  xlab("Longitude") + 
  ylab("Latitude") + 
  theme(plot.title = element_text(hjust = 0.5))
log_intro_native_50
ggsave("/blue/guralnick/millerjared/BoCP/data/processed/grid-project/worldwide/raster-species-PA-visuals/community-richness-50x50/community-richness-log-intro-native-50km.png", width = 10, height = 10)


## Same thing for 25x25
community_pixels_25 <- fread("/blue/guralnick/millerjared/BoCP/data/processed/grid-project/worldwide/raster-species-PA-visuals/community-richness-25x25/community-sp-PA-makeup-intro-native.csv")

community_pixels_25summary <- community_pixels_25 %>% 
  mutate(xy = paste(x, y)) %>% 
  group_by(xy) %>% 
  reframe(x,y, n = sum(presence)) %>% 
  ungroup() %>% 
  select(-xy) 

intro_native_25 <- ggplot() + 
  geom_tile(community_pixels_25summary, mapping = aes(x = x, y = y, fill = n)) + 
  scale_fill_viridis_c(na.value = "transparent") + 
  geom_sf(na_regions, fill = NA, color = "black", size = 0.3, mapping = aes()) +
  theme_minimal() + 
  ggtitle("Introduced & Native \n 25x25km Presence Cells for Plant Community",) + 
  xlab("Longitude") + 
  ylab("Latitude") + 
  theme(plot.title = element_text(hjust = 0.5))
intro_native_25

ggsave("/blue/guralnick/millerjared/BoCP/data/processed/grid-project/worldwide/raster-species-PA-visuals/community-richness-25x25/community-richness-intro-native-25km.png", width = 10, height = 10)

log_intro_native_25 <- ggplot() + 
  geom_tile(community_pixels_25summary, mapping = aes(x = x, y = y, fill = log(n+1))) + 
  scale_fill_viridis_c(na.value = "transparent") + 
  geom_sf(na_regions, fill = NA, color = "black", size = 0.3, mapping = aes()) +
  theme_minimal() + 
  ggtitle("Introduced & Native \n 25x25km Log Transformed Presence Cells \n for Plant Community",) + 
  xlab("Longitude") + 
  ylab("Latitude") +
  theme(plot.title = element_text(hjust = 0.5))
log_intro_native_25 
ggsave("/blue/guralnick/millerjared/BoCP/data/processed/grid-project/worldwide/raster-species-PA-visuals/community-richness-25x25/community-richness-log-intro-native-25km.png", width = 10, height = 10)

### Do the same for just the native data
# read in the community P/A 
community_pixels_N <- fread("/blue/guralnick/millerjared/BoCP/data/processed/grid-project/worldwide/raster-species-PA-visuals/community-richness-50x50/community-sp-PA-makeup.csv")

# set up shapefile (delimited to the same pixels as found in clean-rasters-to-geo-scope.R)
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
# set up mollweide proj centered on NA
moll_proj_center <- "+proj=moll +lon_0=-100 +datum=WGS84 +units=m +no_defs"
# reproj our na regions to this projection 
na_regions <- st_transform(na_regions, crs = crs(moll_proj_center))

# covert community pixels to summary form 
community_pixels_summary <- community_pixels_N %>% 
  mutate(xy = paste(x, y)) %>% 
  group_by(xy) %>% 
  reframe(x,y, n = sum(presence)) %>% 
  ungroup() %>% 
  select(-xy) 




# plot
native_50 <- ggplot() + 
  geom_tile(community_pixels_summary, mapping = aes(x = x, y = y, fill = n)) + 
  scale_fill_viridis_c(na.value = "transparent") + 
  geom_sf(na_regions, fill = NA, color = "black", size = 0.3, mapping = aes()) +
  theme_minimal() + 
  ggtitle("Native Only \n 50x50km Presence Cells for Plant Community",) + 
  xlab("Longitude") + 
  ylab("Latitude") +
  theme(plot.title = element_text(hjust = 0.5))
native_50
ggsave("/blue/guralnick/millerjared/BoCP/data/processed/grid-project/worldwide/raster-species-PA-visuals/community-richness-50x50/community-richness-native-50km.png", width = 10, height = 10)

log_native_50 <- ggplot() + 
  geom_tile(community_pixels_summary, mapping = aes(x = x, y = y, fill = log(n+1))) + 
  scale_fill_viridis_c(na.value = "transparent") + 
  geom_sf(na_regions, fill = NA, color = "black", size = 0.3, mapping = aes()) +
  theme_minimal() + 
  ggtitle("Native Only \n 50x50km Log Transformed Presence Cells \n for Plant Community",) + 
  xlab("Longitude") + 
  ylab("Latitude") +
  theme(plot.title = element_text(hjust = 0.5))
log_native_50
ggsave("/blue/guralnick/millerjared/BoCP/data/processed/grid-project/worldwide/raster-species-PA-visuals/community-richness-50x50/community-richness-log-native-50km.png", width = 10, height = 10)


## Same thing for 25x25
community_pixels_25_N <- fread("/blue/guralnick/millerjared/BoCP/data/processed/grid-project/worldwide/raster-species-PA-visuals/community-richness-25x25/community-sp-PA-makeup.csv")

community_pixels_25summary <- community_pixels_25_N %>% 
  mutate(xy = paste(x, y)) %>% 
  group_by(xy) %>% 
  reframe(x,y, n = sum(presence)) %>% 
  ungroup() %>% 
  select(-xy) 

native_25 <- ggplot() + 
  geom_tile(community_pixels_25summary, mapping = aes(x = x, y = y, fill = n)) + 
  scale_fill_viridis_c(na.value = "transparent") + 
  geom_sf(na_regions, fill = NA, color = "black", size = 0.3, mapping = aes()) +
  theme_minimal() + 
  ggtitle("Native Only \n 25x25km Presence Cells for Plant Community",) + 
  xlab("Longitude") + 
  ylab("Latitude") +
  theme(plot.title = element_text(hjust = 0.5))
native_25

ggsave("/blue/guralnick/millerjared/BoCP/data/processed/grid-project/worldwide/raster-species-PA-visuals/community-richness-25x25/community-richness-native-25km.png", width = 10, height = 10)

log_native_25 <- ggplot() + 
  geom_tile(community_pixels_25summary, mapping = aes(x = x, y = y, fill = log(n+1))) + 
  scale_fill_viridis_c(na.value = "transparent") + 
  geom_sf(na_regions, fill = NA, color = "black", size = 0.3, mapping = aes()) +
  theme_minimal() + 
  ggtitle("Native Only \n 25x25km Log Transformed Presence Cells \n for Plant Community",) + 
  xlab("Longitude") + 
  ylab("Latitude") +
  theme(plot.title = element_text(hjust = 0.5))
log_native_25 
ggsave("/blue/guralnick/millerjared/BoCP/data/processed/grid-project/worldwide/raster-species-PA-visuals/community-richness-25x25/community-richness-log-native-25km.png", width = 10, height = 10)

#########################################################################################################################################################################################################################

# Understand how outliers influence these data
# outlier_df <- fread("./data/processed/grid-project/worldwide/sp-occs-outliers.csv")
# # only consider data that has some presences 
# outlier_df <- outlier_df %>% 
#   filter(num_all_occs_50 > 0) # the coarsest grain 
# 
# ggplot(outlier_df, mapping = aes(x = num_outliers_50)) + 
#   geom_histogram()
# 
# ggplot(outlier_df, mapping = aes(x = num_outliers_25)) + 
#   geom_histogram()
# 

### Additionally, visualize how what net change in pixels there is given introduced + native vs native data pulls
community_pixels_NI <- fread("/blue/guralnick/millerjared/BoCP/data/processed/grid-project/worldwide/raster-species-PA-visuals/community-richness-50x50/community-sp-PA-makeup-intro-native.csv")
community_pixels_N <- fread("/blue/guralnick/millerjared/BoCP/data/processed/grid-project/worldwide/raster-species-PA-visuals/community-richness-50x50/community-sp-PA-makeup.csv")

# retain background to append later
background_rast <- community_pixels_NI %>% 
  filter(species == "background")

# Remove redundant cells 
native_filtered <- community_pixels_N %>% 
  filter(presence == 1) %>% 
  select(x,y, species)

native_introduced_filtered <- community_pixels_NI %>% 
  filter(presence == 1)

# antijoin
introduced_only <- native_introduced_filtered %>% 
  anti_join(native_filtered, by = c("x", "y", "species")) %>% 
  #filter(species == "Lonicera japonica") %>%  # testing, remove later
  rbind(background_rast)

# summary 

introduced_only_summary <- introduced_only %>% 
  mutate(xy = paste(x, y)) %>% 
  group_by(xy) %>% 
  reframe(x,y, n = sum(presence)) %>% 
  ungroup() %>% 
  select(-xy) 

I_50 <- ggplot() + 
  geom_tile(introduced_only_summary, mapping = aes(x = x, y = y, fill = n)) + 
  scale_fill_viridis_c(na.value = "transparent") + 
  geom_sf(na_regions, fill = NA, color = "black", size = 0.3, mapping = aes()) +
  theme_minimal() + 
  ggtitle("Introduced Effect \n 50x50km Presence Cells for Plant Community",) + 
  xlab("Longitude") + 
  ylab("Latitude") + 
  theme(plot.title = element_text(hjust = 0.5))
I_50
ggsave("/blue/guralnick/millerjared/BoCP/data/processed/grid-project/worldwide/raster-species-PA-visuals/community-richness-50x50/community-richness-intro-effect-50km.png", width = 10, height = 10)
log_I_50 <- ggplot() + 
  geom_tile(introduced_only_summary, mapping = aes(x = x, y = y, fill = log(n+1))) + 
  scale_fill_viridis_c(na.value = "transparent") + 
  geom_sf(na_regions, fill = NA, color = "black", size = 0.3, mapping = aes()) +
  theme_minimal() + 
  ggtitle("Introduced Effect \n 50x50km Log Transformed Presence Cells \n for Plant Community",) + 
  xlab("Longitude") + 
  ylab("Latitude") +
  theme(plot.title = element_text(hjust = 0.5))
log_I_50
ggsave("/blue/guralnick/millerjared/BoCP/data/processed/grid-project/worldwide/raster-species-PA-visuals/community-richness-50x50/community-richness-intro-effect-log-50km.png", width = 10, height = 10)

# repeat for 25x25
community_pixels_NI <- fread("/blue/guralnick/millerjared/BoCP/data/processed/grid-project/worldwide/raster-species-PA-visuals/community-richness-25x25/community-sp-PA-makeup-intro-native.csv")
community_pixels_N <- fread("/blue/guralnick/millerjared/BoCP/data/processed/grid-project/worldwide/raster-species-PA-visuals/community-richness-25x25/community-sp-PA-makeup.csv")

# retain background to append later
background_rast <- community_pixels_NI %>% 
  filter(species == "background")

# Remove redundant cells 
native_filtered <- community_pixels_N %>% 
  filter(presence == 1) %>% 
  select(x,y, species)

native_introduced_filtered <- community_pixels_NI %>% 
  filter(presence == 1)

# antijoin
introduced_only <- native_introduced_filtered %>% 
  anti_join(native_filtered, by = c("x", "y", "species")) %>% 
  #filter(species == "Lonicera japonica") %>%  # testing, remove later
  rbind(background_rast)

# summary 

introduced_only_summary <- introduced_only %>% 
  mutate(xy = paste(x, y)) %>% 
  group_by(xy) %>% 
  reframe(x,y, n = sum(presence)) %>% 
  ungroup() %>% 
  select(-xy) 

I_25 <- ggplot() + 
  geom_tile(introduced_only_summary, mapping = aes(x = x, y = y, fill = n)) + 
  scale_fill_viridis_c(na.value = "transparent") + 
  geom_sf(na_regions, fill = NA, color = "black", size = 0.3, mapping = aes()) +
  theme_minimal() + 
  ggtitle("Introduced Effect \n 25x25km Presence Cells for Plant Community",) + 
  xlab("Longitude") + 
  ylab("Latitude") + 
  theme(plot.title = element_text(hjust = 0.5))
I_25
ggsave("/blue/guralnick/millerjared/BoCP/data/processed/grid-project/worldwide/raster-species-PA-visuals/community-richness-25x25/community-richness-intro-effect-25km.png", width = 10, height = 10)
log_I_25 <- ggplot() + 
  geom_tile(introduced_only_summary, mapping = aes(x = x, y = y, fill = log(n+1))) + 
  scale_fill_viridis_c(na.value = "transparent") + 
  geom_sf(na_regions, fill = NA, color = "black", size = 0.3, mapping = aes()) +
  theme_minimal() + 
  ggtitle("Introduced Effect \n 25x25km Log Transformed Presence Cells \n for Plant Community",) + 
  xlab("Longitude") + 
  ylab("Latitude") +
  theme(plot.title = element_text(hjust = 0.5))
log_I_25
ggsave("/blue/guralnick/millerjared/BoCP/data/processed/grid-project/worldwide/raster-species-PA-visuals/community-richness-50x50/community-richness-intro-effect-log-25km.png", width = 10, height = 10)


# patchwork these plots into a pdf
combined_plot <- (native_25 | native_50) /
                 (intro_native_25 | intro_native_50) /
                 (I_25 | I_50) 
ggsave("combined_plots.pdf", combined_plot, width = 12, height = 8)

combined_plot <- (log_native_25 | log_native_50) /
                 (log_intro_native_25 | log_intro_native_50) /
                 (log_I_25 | log_I_50)
ggsave("combined_log_plots.pdf", combined_plot, width = 12, height = 8)


# attempt a faster viewing pdf
library(grid)
library(gridExtra)
pdf("plots_rasterized.pdf", width = 12, height = 8)
grid.arrange(
  rasterGrob(png::readPNG("./data/processed/grid-project/worldwide/raster-species-PA-visuals/community-richness-25x25/community-richness-native-25km.png")),
  rasterGrob(png::readPNG("./data/processed/grid-project/worldwide/raster-species-PA-visuals/community-richness-50x50/community-richness-native-50km.png")),
  rasterGrob(png::readPNG("./data/processed/grid-project/worldwide/raster-species-PA-visuals/community-richness-25x25/community-richness-intro-native-25km.png")),
  rasterGrob(png::readPNG("./data/processed/grid-project/worldwide/raster-species-PA-visuals/community-richness-50x50/community-richness-intro-native-50km.png")),
  rasterGrob(png::readPNG("./data/processed/grid-project/worldwide/raster-species-PA-visuals/community-richness-25x25/community-richness-intro-effect-25km.png")),
  rasterGrob(png::readPNG("./data/processed/grid-project/worldwide/raster-species-PA-visuals/community-richness-50x50/community-richness-intro-effect-50km.png")),
  ncol = 2
)
dev.off()

pdf("log_plots_rasterized.pdf", width = 8, height = 8)
grid.arrange(
  rasterGrob(png::readPNG("./data/processed/grid-project/worldwide/raster-species-PA-visuals/community-richness-25x25/community-richness-log-native-25km.png")),
  rasterGrob(png::readPNG("./data/processed/grid-project/worldwide/raster-species-PA-visuals/community-richness-50x50/community-richness-log-native-50km.png")),
  rasterGrob(png::readPNG("./data/processed/grid-project/worldwide/raster-species-PA-visuals/community-richness-25x25/community-richness-log-intro-native-25km.png")),
  rasterGrob(png::readPNG("./data/processed/grid-project/worldwide/raster-species-PA-visuals/community-richness-50x50/community-richness-log-intro-native-50km.png")),
  rasterGrob(png::readPNG("./data/processed/grid-project/worldwide/raster-species-PA-visuals/community-richness-25x25/community-richness-log-intro-effect-25km.png")),
  rasterGrob(png::readPNG("./data/processed/grid-project/worldwide/raster-species-PA-visuals/community-richness-50x50/community-richness-log-intro-effect-50km.png")),
  ncol = 2
)
dev.off()