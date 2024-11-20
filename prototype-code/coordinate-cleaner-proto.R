# 003C-coordinate-clean
# Author: JT Miller 
# Date: 08-02-2024
# Project: BoCP 

## Cleaning coordinates for flagged errorneous records + botanical gardens/institutions. 

# load libraries 
library(dplyr)
library(data.table)
library(CoordinateCleaner)

# Load modified gatoRs functions
source("finished-code/r/gatoRs-fxns-edited.R") # edited gatoRs fxns to allow for only pulling idigbio

# Read in the name alignment
name_alignment <- fread("./data/processed/finalized_name_alignment_wcvp.csv")
# Create an accepted name vector & filestyle version.
accepted_name_v <- unique(name_alignment$acceptedNameParent)
# index when necessary
accepted_name_v <- accepted_name_v
accepted_name_filestyle_v <- gsub(" ", "-", accepted_name_v)

# load names and take distinct records for relevant fields.
for(i in 1:length(accepted_name_v)){
  if(file.exists(paste0("./data/processed/distinct-clean-data/", accepted_name_filestyle_v[i], ".csv"))){
    occur_data <- fread(paste0("./data/processed/distinct-clean-data/", accepted_name_filestyle_v[i], ".csv"))
  } else{
    print(paste("No occur data for", accepted_name_v[[i]]))
  } 
  
  records_prior_coord_cleaning <- nrow(occur_data)
  
  # Run coord clean using all of the defaults besides outliers (I'd prefer to handle these in-house.)
  occur_data <- CoordinateCleaner::clean_coordinates(occur_data, 
                                                     lon = "longitude", 
                                                     lat = "latitude", 
                                                     species = "scientificNameParsed", 
                                                     tests = c("capitals", "centroids", "equal", "gbif", "institutions", "seas",
                                                               "zeros")
                                                     )
  # remove records that do not pass coordinate clean checks                                                   
  occur_data <- occur_data %>% 
    filter(.summary == TRUE) %>% 
    select(-.val, -.equ, -.zer, -.cap, -.cen, -.sea, -.gbf, -.inst, -.summary)
  records_post_coord_cleaning <- nrow(occur_data)

  # create table 
  coord_clean_info_tb <- data.frame(
    name = accepted_name_v[i],
    priorCoordClean = records_prior_coord_cleaning,
    postCoordClean = records_post_coord_cleaning
  )
  fwrite(occur_data, paste0("./data/processed/coord-clean-data/", accepted_name_filestyle_v[i], ".csv"))
  fwrite(coord_clean_info_tb , file = "./data/processed/cleaning-summaries/coord_clean_info_tb.csv", append = TRUE, col.names = !file.exists("./data/processed/cleaning-summaries/distinct_clean_info_tb.csv"))
  print(i) # print where we're at in the loop...
}