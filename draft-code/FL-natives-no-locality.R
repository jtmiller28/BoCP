### A script for determining florida plants from our datasets that lack occurrence data 
# Author: JT Miller 
# Date: 10-31-2024
# Project: BoCP 

# Load Libraries
library(data.table)
library(dplyr)
library(sf)

# set dir 
setwd("/home/millerjared/blue_guralnick/millerjared/BoCP/")

# Set up array logic
start_num <- as.numeric(Sys.getenv("START_NUM"))
task_id <- as.numeric(start_num)
part <- paste0("part", task_id)

# load data 
wcvp_distribution_tb <- fread("./data/raw/wcvp_distribution.csv") # distributional data
wcvp_name_tb <- fread("./data/raw/wcvp_names.csv")
na_plants_alignment <- fread("./data/processed/finalized_name_alignment_wcvp.csv")
level3_regions <- read_sf("./data/raw/level3-wgsrpd/level3.shp")

# select for relevant info 
wcvp_names <- wcvp_name_tb %>% 
  select(plant_name_id, taxon_name, taxon_authors, taxon_status)

# merge names with distribution info.
wcvp_names_w_dist <- merge(wcvp_names, wcvp_distribution_tb, by = "plant_name_id")

na_wcvp_names_w_dist <- wcvp_names_w_dist %>% 
  filter(taxon_name %in% na_plants_alignment$acceptedName) # use acceptedName as this will create the full range for the parentTaxon

# select for relevant info
na_wcvp_names_w_dist <- select(na_wcvp_names_w_dist, plant_name_id, taxon_name, continent, region_code_l2, region, area_code_l3,
                               area, introduced, extinct, location_doubtful)

# determine parentTaxon using the name alignment's criteria 
na_wcvp_parentTaxon_w_dist <- na_wcvp_names_w_dist %>% 
  left_join(na_plants_alignment, by = join_by(taxon_name == acceptedName)) %>% 
  rename(parentTaxon = acceptedNameParent) %>% 
  select(names(na_wcvp_names_w_dist), parentTaxon)

# Limit to florida taxon
fl_wcvp_parentTaxon_w_dist <- filter(na_wcvp_parentTaxon_w_dist, area =="Florida")

# build a vector of names to index by, then index according to array tasks
accepted_name_v <- unique(na_plants_alignment$acceptedNameParent)
accepted_name_fl_v <- intersect(accepted_name_v, fl_wcvp_parentTaxon_w_dist$parentTaxon)
accepted_name <- accepted_name_v[task_id] # index by task ID
accepted_name_filestyle <- gsub(" ", "-", accepted_name)

# read in the data of the names we know we have data for
fl_names_w_data <- fread("./data/processed/fl-plants-project/name-summary/fl-plant-summary.csv")


fl_names_w_out_occur_data <- setdiff(accepted_name_fl_v, fl_names_w_data$parentTaxon)

fl_names_w_out_occur_data_dt <- as.data.table(fl_names_w_out_occur_data)
fl_names_w_out_occur_data_dt <- rename(fl_names_w_out_occur_data_dt, parentTaxon = fl_names_w_out_occur_data)

# check to see if they are endemic or not



fwrite(fl_names_w_out_occur_data_dt, "./data/processed/fl-plants-project/name-summary/fl-plants-w-out-occur-data.csv")


