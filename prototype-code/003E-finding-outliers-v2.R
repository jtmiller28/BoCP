### A script for identifying outliers v2
# Author: JT Miller 
# Date: 09-26-2024
# Project: BoCP 

# Load packages
library(data.table)
library(sf)
library(dplyr)
library(ggplot2) 
library(geosphere)
library(dbscan)

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

# DBSCAN clustering method. 

# convert occurrence data to matrix for clustering

coords_matrix <- as.matrix(occur_data[, .(roundedLongitude, roundedLatitude)])

# find suitable DBSCAN params
kNNdistplot(coords_matrix, minPts = 5)
# run DBSCAN using the Haversine distance matrix
clustering <- dbscan::dbscan(coords_matrix, eps = 0.1, minPts = 5) # eps is the neighborhood radius, minpts is the min required to form a cluster. 

# add cluster labels
occur_data$cluster <- as.factor(clustering$cluster)


