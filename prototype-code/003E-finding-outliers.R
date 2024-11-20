### A script for identifying outliers
# Author: JT Miller 
# Date: 09-26-2024
# Project: BoCP 

# Load packages
library(data.table)
library(sf)
library(dplyr)
library(ggplot2)
library(geosphere)

# set dir
setwd("/home/millerjared/blue_guralnick/millerjared/BoCP/")


# use name alignment to read in data
na_plants_alignment <- fread("./data/processed/finalized_name_alignment_wcvp.csv")

# load shapefile
level3_regions <- read_sf("./data/raw/level3-wgsrpd/level3.shp")

# build a vector of names to index by, then index according to array tasks
accepted_name_v <- unique(na_plants_alignment$acceptedNameParent)
accepted_name <- accepted_name_v[task_id] # index by task ID
accepted_name_filestyle <- gsub(" ", "-", accepted_name)

# read in these names 
  if(file.exists(paste0("./data/processed/regioned-data/native/", accepted_name_filestyle, ".csv"))){
  
  # load the data
  occur_data <- fread(paste0("./data/processed/regioned-data/native/", accepted_name_filestyle, ".csv"))
} else { # If the name does not exist at this step, include this in the final output table. 
  print(paste("No occur data for", accepted_name))
} # end of if statement


# Calculate the pairwise dists
dist_matrix <- distm(as.matrix(occur_data[, .(roundedLongitude, roundedLatitude)]), fun = distHaversine)

# use only the upper triangle of dists
upper_triangle <- upper.tri(dist_matrix, diag = FALSE)

# flatten the upper triangle to obtain the pairwise dists
dist_vector <- dist_matrix[upper_triangle]
# associate pairwise dists with uuids
uuid_combinations <- expand.grid(uuid1 = occur_data$uuid, uuid2 = occur_data$uuid)
uuid_combinations <- uuid_combinations[upper_triangle, ] # create all combinations of points except for identical

# convert to dataframe for ggploting
dist_df <- data.frame(uuid_combinations, distance = dist_vector) 
dist_df <- dist_df %>% mutate(distance_km = distance/1000)

# calc the mean, standard deviation, and the upper threshold (lower is not needed as I dont really care if things are close)
# mean_dist <- mean(dist_df$distance_km, na.rm = TRUE)
# sd_dist <- sd(dist_df$distance_km, na.rm = TRUE)
# upper_threshold <- mean_dist + 2 * sd_dist

# calculate Q1, Q3, and IQR
Q1 <- quantile(dist_df$distance_km, 0.25, na.rm = TRUE)
Q3 <- quantile(dist_df$distance_km, 0.75, na.rm = TRUE)
IQR_value <- IQR(dist_df$distance_km, na.rm = TRUE)

# calculate the upper threshold for outliers
upper_threshold <- Q3 + 1.5 * IQR_value

# identify points that increase the distance past the upper threshold
high_dist_df <- dist_df %>% filter(distance_km > upper_threshold)

# find uuids that cause this
high_dist_points <- unique(c(high_dist_df$uuid1, high_dist_df$uuid2))

# demonstrate what removal will look like
ggplot(dist_df, aes(x = distance_km)) + 
  geom_histogram(binwidth = 100, fill = "lightblue", color = "black") + # Adjust binwidth
  geom_vline(aes(xintercept = upper_threshold, color = "Cut-off threshold"), linetype = "dashed", size = 1, show.legend = TRUE) + 
  labs(title = "Distribution of Pairwise Distances", x = "Distance (km)", y = "Count") + 
  scale_color_manual(name = "Legend", values = c("Cut-off threshold" = "red")) + 
  theme_minimal()

# boxplot with horizontal orientation
boxplot(dist_df$distance_km, horizontal = TRUE, main = "Boxplot of Pairwise Distances", xlab = "Distance (km)")

# sdd red points for those above the upper threshold
outliers <- dist_df$distance_km[dist_df$distance_km > upper_threshold]
outlier_positions <- which(dist_df$distance_km > upper_threshold)

# sdd red points for the outliers
points(outliers, rep(1, length(outliers)), col = "red", pch = 16)

# sdd a dashed vertical line indicating the threshold
abline(v = upper_threshold, col = "red", lty = 2)

occur_data_sf <- st_as_sf(occur_data, coords = c("roundedLongitude", "roundedLatitude"), crs = 4326)
occur_data_sf <- mutate(occur_data_sf, high_dist_pts = uuid %in% high_dist_points)

# Create convex hull and bounding box
hull <- st_convex_hull(st_union(occur_data_sf))
hull_bbox <- st_bbox(hull)

# Plot spatial data with points marked for removal in red
ggplot() +
  geom_sf(data = level3_regions, mapping = aes()) +  # Add the background level3 regions
  geom_sf(data = occur_data_sf, mapping = aes(color = high_dist_pts)) +  # Map fill based on whether the point is to be removed
  scale_color_manual(name = "Points Removal", values = c("TRUE" = "red", "FALSE" = "black")) +  # Use red for points flagged for removal
  coord_sf(xlim = c(hull_bbox[1], hull_bbox[3]), ylim = c(hull_bbox[2], hull_bbox[4])) +  # Zoom to convex hull
  theme_minimal()



  