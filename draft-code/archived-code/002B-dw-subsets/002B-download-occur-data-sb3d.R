# 002-download-occur-data Rewrite
# Author: JT Miller 
# Date: 07/08/2024
# Project: BoCP 

# Load packages
library(gatoRs)
library(data.table)
library(dplyr)

# set working dir
setwd("/blue/guralnick/millerjared/BoCP")
# call edited gatoRs fxns
source("finished-code/r/gatoRs-fxns-edited.R") # edited gatoRs fxns to allow for only pulling gbif

# Read in the name alignment
name_alignment <- fread("./data/processed/finalized_name_alignment_wcvp.csv")

# Prep name list
## Note that we use acceptedNameParents because we do not plan to model at the infraspecific level. 
## Structure these names in a list with acceptedNameParent + synonym(s)
accepted_name_v <- unique(name_alignment$acceptedNameParent)
### Subset to do only a portion of names for this script. 
accepted_name_v <- accepted_name_v[4501:6500]
accepted_name_v <- accepted_name_v[1751:2000]
## Make a file style version of this vector for storage purpsoses
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

## Download with retry attempts
# Initialize a list in case failures occur
failed_names_holder <- list()

# Recursive function for retrying download
retry_download <- function(i, retry_count) {
  if (retry_count >= 20) {
    # Customize the error handler for max retries exceeded
    failed_names_holder[[length(failed_names_holder) + 1]] <<- c(accepted_name_v[[i]], "Maximum retries exceeded", format(Sys.time(), "%a %b %d %X %Y"))
    failed_names_df <- as.data.frame(do.call(rbind, failed_names_holder))
    colnames(failed_names_df) <- c("acceptedParentName", "errorMsg", "timeStamp")
    fwrite(failed_names_df, "./data/gbif-dws/raw-dw-info-tbls/failed_names3.csv", append = TRUE, col.names = FALSE)
    return(NULL)
  }
  
  # Try and Catch errors
  result <- tryCatch({
    occurrence_raw_dw <- gators_download_edited(name_list[[i]], call.idigbio = FALSE, call.gbif = TRUE)
    return(occurrence_raw_dw)
  }, error = function(e) {
    if (e$message != "No Records Found.") {
      Sys.sleep(retry_count * 15)
      print(paste("Download attempt", retry_count + 1, "for", accepted_name_v[[i]], "failed. Retrying with delay", retry_count * 30, "second delay"))
      return(retry_download(i, retry_count + 1))
    } else {
      # Only record the failure once per `i`
      if (!any(sapply(failed_names_holder, function(x) x[1] == accepted_name_v[[i]]))) {
        failed_names_holder[[length(failed_names_holder) + 1]] <<- c(accepted_name_v[[i]], e$message, format(Sys.time(), "%a %b %d %X %Y"))
        failed_names_df <- as.data.frame(do.call(rbind, failed_names_holder))
        colnames(failed_names_df) <- c("acceptedParentName", "errorMsg", "timeStamp")
        fwrite(failed_names_df, "./data/gbif-dws/raw-dw-info-tbls/failed_names3.csv", append = TRUE, col.names = FALSE)
      }
      return(NULL)
    }
  })
  
  return(result)
}

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

# Proceed to download occurrences 
for (j in 1:length(accepted_name_v)) {
  time_taken <- bench::bench_time({
    occurrence_data <- retry_download(j, 0)  # Initial call with retry_count set to 0
  }) # end of bench time
  if (!is.null(occurrence_data)) {
    info_df_temp <- data.frame(
      acceptedParentName = accepted_name_v[[j]], # Store the acceptedNameParent 
      numOccurrences = nrow(occurrence_data), # Store information on the number of occurrence records in raw download
      timeTaken = time_taken[["real"]] # Store the time (wall time) taken
    )
    fwrite(info_df_temp, file = "./data/gbif-dws/raw-dw-info-tbls/raw-info-tbl3.csv", append = TRUE, col.names = !file.exists("./data/gbif-dws/raw-dw-info-tbls/raw-info-tbl3.csv"))
    # Process the occurrence data if not NULL
    print(paste("Calling data for", accepted_name_v[[j]]))
    fwrite(occurrence_data, file = paste0("./data/gbif-dws/raw/", accepted_name_filestyle_v[[j]], ".csv"))
  } else {
    print(paste("No data found for", accepted_name_v[[j]]))
  }
  print_loading_bar(j, total_names)  # Update the loading bar
}
cat("\n")  # Move to a new line after the loading bar is complete





