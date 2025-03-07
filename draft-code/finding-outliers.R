### A script for explaining how we are identifying outliers
# Author: JT Miller 
# Date: 09-26-2024
# Project: BoCP 

library(data.table)
library(sf)
library(dplyr)
library(ggplot2)
library(CoordinateCleaner)
library(geosphere)
library(cowplot)

# Create a reproducible dataset of clustered points: 
## Set seed for reproducibility
set.seed(123)

## Load shp file
level3_regions <- read_sf("./data/raw/level3-wgsrpd/level3.shp")
## Define a function to generate clustered points with varying points per cluster
generate_clustered_points <- function(n_clusters, points_range, lon_range, lat_range, spread) {
  points <- data.frame(lon = numeric(0), lat = numeric(0))
  
  # Generate cluster centers
  cluster_centers <- data.frame(
    lon = runif(n_clusters, lon_range[1], lon_range[2]),
    lat = runif(n_clusters, lat_range[1], lat_range[2])
  )
  
  # Generate points around each cluster center
  for (i in 1:n_clusters) {
    # Randomly select the number of points for this cluster
    points_per_cluster <- sample(points_range[1]:points_range[2], 1)
    
    # Generate points for the cluster
    lon_vals <- rnorm(points_per_cluster, mean = cluster_centers$lon[i], sd = spread)
    lat_vals <- rnorm(points_per_cluster, mean = cluster_centers$lat[i], sd = spread)
    cluster_points <- data.frame(lon = lon_vals, lat = lat_vals)
    
    # Add points to the dataframe
    points <- rbind(points, cluster_points)
  }
  
  return(points)
}

# Parameters for clustering
n_clusters <- 5 # Number of clusters
points_range <- c(3, 25) # Range of points per cluster (min, max)
lon_range <- c(-120, -80) # Longitude range for North America (approx)
lat_range <- c(35, 50)    # Latitude range for North America (approx)
spread <- 5 # Standard deviation of the cluster distance

# Generate clustered points with varying points per cluster
points_df <- generate_clustered_points(n_clusters, points_range, lon_range, lat_range, spread)

# Add species column
points_df <- points_df %>%
  mutate(species = "species x") %>% 
  mutate(id = 1:n())
# Convert to an sf object with WGS84 projection
points_sf <- st_as_sf(points_df, coords = c("lon", "lat"), crs = 4326, remove = FALSE)

# View first few rows
print(head(points_sf))

# Plot to visualize the clustered points
ggplot() +
  geom_sf(data = points_sf, aes(), color = "black") +
  theme_minimal() +
  labs(title = "Sample Distribution for 'species x'",
       x = "Longitude", y = "Latitude")

## Demonstrate cc_outlier functionality 

# First demonstrating the the the distance method: The minimum distance is set by tdi in km, if another point does not fall within this distance, drop it. 
id_dist <- CoordinateCleaner::cc_outl(points_df, # dont use the sf version for this, will throw error due to the geom col
                                      lon = "lon",
                                      lat = "lat", 
                                      species = "species", 
                                      method = "distance", # options are "distance", "quantile", and "mad"
                                      mltpl = 5, # The multiplier of the interquantile range (for method = "quantile" or median absolute deviation (for method = "mad)
                                      tdi = 100, # the minimum absolute distance of a record to all other records of a species to be IDed as an outlier in km
                                      value = "flagged", # determines type of output being cleaned dataset ("clean") or flagged "flagged"
                                      sampling_thresh = 0, # the cut off threshold for the sampling correction. Indicates the quantile of sampling in which outliers should be ignored. 
                                      verbose = TRUE, 
                                      min_occs = 7, # minimum number of geographically unique datapoints needed for a species to be tested. 
                                      thinning = FALSE, # forces a raster approximation for the distance calc. Useful for species with more than 10,000 records to speed up computation or for small datasets that have very uneven sampling in their datasets. 
                                      thinning_res = 0.5 # the resolution for the spatial thinning in decimal degrees, default = 0.5
)

id_dist <- !id_dist # because its annoying to have them ordered opposite...
points_sf$distOutlierFlagged <- id_dist # append the outleirs as a boolean field 

# Create convex hull and bounding box
hull <- st_convex_hull(st_union(points_sf))
hull_bbox <- st_bbox(hull)

quantile_id <- unique(points_sf$quantile)

# create some buffers to illustrate distance overlaps. 
points_sf$buffer <- st_buffer(points_sf, dist = 100000) # 100km to match the tdi

buffer_fills <- c("green", "red")

# Corrected Plot
ggplot() +
  geom_sf(data = level3_regions, color = "black") +  # Load the level3 regions correctly
  geom_sf(data = points_sf$buffer, aes(fill = distOutlierFlagged), alpha = 0.2) +  # Use aes for mapping colors
  geom_sf(data = points_sf, aes(color = distOutlierFlagged)) +  # Map color for points flagged for removal
  scale_fill_manual(name = "Points Removal", values = c("FALSE" = "green", "TRUE" = "red")) +  # Fill for buffer circles
  scale_color_manual(name = "Points Removal", values = c("FALSE" = "black", "TRUE" = "red")) +  # Color for points
  coord_sf(xlim = c(hull_bbox[1], hull_bbox[3]), ylim = c(hull_bbox[2], hull_bbox[4])) +  # Zoom to convex hull
  theme_minimal() +
  ggtitle(paste0("Native Distribution for species x\nUsing Distance Method with 100km Outlier Threshold"))

# generally from messing around with the parameters it seems that the distance method performs well when there is a large amount of data, or the data is well clustered. If there a low data density/the data is in sporadic large clusters this method fails. 

## Second demonstrating the the the mad method: "Median Absolute Deviation", a record is flagged as an outlier if the mean distance to all other records is larger than the median of the mean distance of all points plus/minus the mad of the mean distances of all records of the species * mltpl
# distance calc
dist <- geosphere::distm(points_df[, c("lon", "lat")], 
                         fun = geosphere::distHaversine) / 1000 

# set diagonale to NA, so that it does not influence the mean
dist[dist==0] <- NA

# calc minimum distance to any next pt 
mins <- apply(dist, 1, min, na.rm = TRUE)

# obtain the median 
quo <- stats::median(mins, na.rm = TRUE)
# get the mad stat 
tester <- stats::mad(mins, na.rm = TRUE)

# define mltpl
mltpl <- 5

# identify outliers 
outlier_indices <- which(mins > (quo + tester * mltpl))

# flag points as outliers 
points_df <- points_df %>% 
  mutate(is_outlier = id %in% outlier_indices)
 
cutoff_upper <- quo+tester*mltpl
cutoff_lower <- quo-tester*mltpl

# Plot histogram of minimum distances
mad_hist_plot <- ggplot(data = data.frame(min_distance = mins), aes(x = min_distance)) +
  geom_histogram(binwidth = 10, color = "black", fill = "lightblue") +
  labs(title = "Histogram of Minimum Distances",
       x = "Minimum Distance (km)",
       y = "Frequency") +
  geom_vline(xintercept = cutoff_upper, linetype = "dashed", color = "red") +
  geom_vline(xintercept = cutoff_lower, linetype = "dashed", color = "blue") +
  annotate("text", x = cutoff_upper, y = 1, label = paste("Upper Cutoff: ", round(cutoff_upper, 2)), 
           vjust = -1, color = "red") +
  annotate("text", x = cutoff_lower, y = 1, label = paste("Lower Cutoff: ", round(cutoff_lower, 2)), 
           vjust = -1, color = "blue") +
  theme_minimal()

id_mad <- CoordinateCleaner::cc_outl(points_df, # dont use the sf version for this, will throw error due to the geom col
                                      lon = "lon",
                                      lat = "lat", 
                                      species = "species", 
                                      method = "mad", # options are "distance", "quantile", and "mad"
                                      mltpl = 5, # The multiplier of the interquantile range (for method = "quantile" or median absolute deviation (for method = "mad)
                                      tdi = 100, # the minimum absolute distance of a record to all other records of a species to be IDed as an outlier in km
                                      value = "flagged", # determines type of output being cleaned dataset ("clean") or flagged "flagged"
                                      sampling_thresh = 0, # the cut off threshold for the sampling correction. Indicates the quantile of sampling in which outliers should be ignored. 
                                      verbose = TRUE, 
                                      min_occs = 7, # minimum number of geographically unique datapoints needed for a species to be tested. 
                                      thinning = FALSE, # forces a raster approximation for the distance calc. Useful for species with more than 10,000 records to speed up computation or for small datasets that have very uneven sampling in their datasets. 
                                      thinning_res = 0.5 # the resolution for the spatial thinning in decimal degrees, default = 0.5
)

id_mad <- !id_mad # because its annoying to have them ordered opposite...
points_sf$madOutlierFlagged <- id_mad # append the outleirs as a boolean field 

# Create convex hull and bounding box
hull <- st_convex_hull(st_union(points_sf))
hull_bbox <- st_bbox(hull)


# Plot the MAD method results
mad_spatial_plot <- ggplot() +
  geom_sf(data = level3_regions, color = "black") +  # Level3 regions background
  geom_sf(data = points_sf, aes(color = madOutlierFlagged)) +  # Points flagged as outliers
  scale_fill_manual(name = "Outlier Buffer", values = c("FALSE" = "green", "TRUE" = "red")) +  # Buffers fill
  scale_color_manual(name = "Points Removal", values = c("FALSE" = "black", "TRUE" = "red")) +  # Point color
  coord_sf(xlim = c(hull_bbox[1], hull_bbox[3]), ylim = c(hull_bbox[2], hull_bbox[4])) +  # Zoom to convex hull
  theme_minimal() 

combined_plot <- plot_grid(mad_spatial_plot, mad_hist_plot, labels = c("A", "B"), ncol = 2)

combined_plot

## Third demonstrating the the quantile: a boxplot method is used and records are flagged as outliers if their mean distance to all other records of thespecies is larger than mltpl * the interquantile range of the mean distance of all records for the species. 
dist <- geosphere::distm(points_df[, c("lon", "lat")], 
                         fun = geosphere::distHaversine) / 1000 

# set diagonale to NA, so that it does not influence the mean
dist[dist==0] <- NA

# calc minimum distance to any next pt 
mins <- apply(dist, 1, min, na.rm = TRUE)

quo <- quantile(mins, c(0.25, 0.75), na.rm = TRUE)

# flag outliers 
iqr <- IQR(mins, na.rm = TRUE)

# flag outliers using the upper bound criterion
out <- which(mins > (quo[2] + iqr*mltpl))

id_quantile <- !seq_along(mins) %in% out
points_sf$quantileOutlierFlagged <- id_quantile

# Create the boxplot of minimum distances
boxplot_quantile <- ggplot(data = data.frame(min_distance = mins), aes(y = min_distance)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 2) +
  labs(title = "Boxplot of Minimum Distances",
       y = "Minimum Distance (km)") +
  theme_minimal() +
  geom_hline(yintercept = quo[2] + iqr * mltpl, linetype = "dashed", color = "red") +
  annotate("text", x = 1, y = quo[2] + iqr * mltpl + 1, 
           label = paste("Upper Cutoff: ", round(quo[2] + iqr * mltpl, 2)), 
           vjust = -0.5, color = "red")

id_quantile <- CoordinateCleaner::cc_outl(points_df, # dont use the sf version for this, will throw error due to the geom col
                                     lon = "lon",
                                     lat = "lat", 
                                     species = "species", 
                                     method = "quantile", # options are "distance", "quantile", and "mad"
                                     mltpl = 5, # The multiplier of the interquantile range (for method = "quantile" or median absolute deviation (for method = "mad)
                                     tdi = 100, # the minimum absolute distance of a record to all other records of a species to be IDed as an outlier in km
                                     value = "flagged", # determines type of output being cleaned dataset ("clean") or flagged "flagged"
                                     sampling_thresh = 0, # the cut off threshold for the sampling correction. Indicates the quantile of sampling in which outliers should be ignored. 
                                     verbose = TRUE, 
                                     min_occs = 7, # minimum number of geographically unique datapoints needed for a species to be tested. 
                                     thinning = FALSE, # forces a raster approximation for the distance calc. Useful for species with more than 10,000 records to speed up computation or for small datasets that have very uneven sampling in their datasets. 
                                     thinning_res = 0.5 # the resolution for the spatial thinning in decimal degrees, default = 0.5
)

id_quantile <- !id_quantile # because its annoying to have them ordered opposite...
points_sf$quantileOutlierFlagged <- id_quantile # append the outleirs as a boolean field 

# Create convex hull and bounding box
hull <- st_convex_hull(st_union(points_sf))
hull_bbox <- st_bbox(hull)


# Plot the MAD method results
quantile_spatial_plot <- ggplot() +
  geom_sf(data = level3_regions, color = "black") +  # Level3 regions background
  geom_sf(data = points_sf, aes(color = quantileOutlierFlagged)) +  # Points flagged as outliers
  scale_fill_manual(name = "Outlier Buffer", values = c("FALSE" = "green", "TRUE" = "red")) +  # Buffers fill
  scale_color_manual(name = "Points Removal", values = c("FALSE" = "black", "TRUE" = "red")) +  # Point color
  coord_sf(xlim = c(hull_bbox[1], hull_bbox[3]), ylim = c(hull_bbox[2], hull_bbox[4])) +  # Zoom to convex hull
  theme_minimal() 

combined_plot_quantile <- plot_grid(boxplot_quantile, quantile_spatial_plot)

combined_plot_quantile
