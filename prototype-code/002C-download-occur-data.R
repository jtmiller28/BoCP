# 002-download Symbiota Data
# Author: JT Miller 
# Date: 07/08/2024
# Project: BoCP 

# Load packages
library(data.table)
library(dplyr)

# read in data on vascular plants obtained from Ed Gilbert at Symbiota. Append colnames as those didnt quite make it on the data
colnames_df <- fread("/blue/guralnick/millerjared/BoCP/data/symbiota-dws/symbiota/field_heading.csv")
new_cols <- names(colnames_df)
cch2 <- fread("/blue/guralnick/millerjared/BoCP/data/symbiota-dws/symbiota/cch2_export.csv")
colnames(cch2) <- new_cols
cnh <- fread("/blue/guralnick/millerjared/BoCP/data/symbiota-dws/symbiota/cnh_export.csv")
colnames(cnh) <- new_cols
sienet <- fread("/blue/guralnick/millerjared/BoCP/data/symbiota-dws/symbiota/seinet_export.csv")
colnames(sienet) <- new_cols

# Combine data
symbiota_data <- rbind(cch2, cnh, sienet)

# Restructure NAs (currently are \Ns)
symbiota_data[symbiota_data == "\\N"] <- NA

# The sciname and authorship fields appear to be preprocessed, so we'll move on. 
# Index names to pull from the data. 
# Read in the name alignment
name_alignment <- fread("./data/processed/finalized_name_alignment_wcvp.csv")
accepted_name_v <- unique(name_alignment$acceptedNameParent)
accepted_name_filestyle_v <- gsub(" ", "-", accepted_name_v)
name_list <- list() # initialize empty list
## Reorganize names into a nested list 
bench::bench_time({
  for(i in 1:length(accepted_name_v)){
    p <- name_alignment[acceptedNameParent == accepted_name_v[[i]]] # create a for-loop that goes through by the accepted name and grabs synonyms, storing them both as a vector in a list. 
    s <- p[!is.na(synonym)]
    x <- print(s$synonym)
    a <- print(unique(p$acceptedNameParent))
    name_list[[i]] <- unique(c(a,x))
  }
})

# Function to print the loading bar
print_loading_bar <- function(current, total) {
  width <- 50
  progress <- current / total
  bar <- paste(rep("=", floor(progress * width)), collapse = "")
  space <- paste(rep(" ", width - floor(progress * width)), collapse = "")
  percentage <- sprintf("%3.0f%%", progress * 100)
  cat(sprintf("\r|%s%s| %s", bar, space, percentage))
}

# init name length
total_names <- length(accepted_name_v)


failed_names_holder <- list()
for(i in 1:length(accepted_name_v)){
  occurrence_raw_dat <- symbiota_data[sciname %in% name_list[[i]]] 
  if(nrow(occurrence_raw_dat) > 0){
    info_df_temp <- data.frame(
      acceptedParentName = accepted_name_v[[i]], # Store the acceptedNameParent 
      numOccurrences = nrow(occurrence_raw_dat)) # Store information on the number of occurrence records in raw download
    fwrite(info_df_temp, file = "./data/symbiota-dws/raw-info-tbls.csv", append = TRUE, col.names = !file.exists("./data/symbiota-dws/raw-info-tbls.csv"))
    fwrite(occurrence_raw_dat, file = paste0("./data/symbiota-dws/raw/", accepted_name_filestyle_v[[i]], ".csv"))
  }
  else{
  }
  print_loading_bar(i, total_names)  # Update the loading bar
}
cat("\n")  # Move to a new line after the loading bar is complete

