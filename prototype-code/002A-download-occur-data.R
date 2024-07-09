# 002-download-occur-data Rewrite
# Author: JT Miller 
# Date: 07/08/2024
# Project: BoCP 

# Load packages
library(gatoRs)
library(data.table)
library(dplyr)
source("prototype-code/gatoRs-fxns-edited.R") # edited gatoRs fxns to allow for only pulling idigbio

# Read in the name alignment
name_alignment <- fread("/home/jtmiller/my_elements/jtmiller/BoCP/data/processed/finalized_name_alignment_wcvp.csv")

# Prep name list
## Note that we use acceptedNameParents because we do not plan to model at the infraspecific level. 
## Structure these names in a list with acceptedNameParent + synonym(s)
accepted_name_v <- unique(name_alignment$acceptedNameParent)
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

# Build a retry fxn in case we overping the API 
failed_names_holder <- list() # Init a list in case failures occur
retry_download <- function(i, retry_count) { # Build a retry_downloader fxn that proceeds to try if an API error occurs
  if (retry_count >= 20) { # Allow up to 20 attempts (maximum of 10 minute cooling period)
    failed_names_holder[[i]] <<- c(accepted_name_v[[i]], "Maximum retries exceeded", format(Sys.time(), "%a %b %d %X %Y")) # Customize the error handler 
    failed_names_df <- as.data.frame(do.call(rbind, failed_names_holder))
    colnames(failed_names_df) <- c("acceptedParentName", "errorMsg", "timeStamp")
    fwrite(failed_names_df, "/home/jtmiller/my_elements/jtmiller/BoCP/data/idigbio-dws/raw-dw-info-tbls/failed_names.csv", append = TRUE, col.names = FALSE)
    return(NULL) # Null otherwise 
  } # end of if condiitonal
  
  tryCatch({ # Try and Catch errors
    occurrence_raw_dw <- gators_download_edited(name_list[[i]], 
                                                call.idigbio = TRUE, 
                                                call.gbif = FALSE)
    
    return(occurrence_raw_dw) # Return the raw dw
    print("No Null Found")
  }, error = function(e) { # create error report
    if (e$message != "No records found.") { # If errors besides the following occurs, retry
      Sys.sleep(retry_count*30)  # Apply time-cooling period: half a minute between API tries...
      print(paste("Download attempt", retry_count + 1, "for", accepted_name_v[[i]], "failed. Retrying with delay", print(retry_count*30), "second delay"))
      return(retry_download(i, retry_count + 1))  # Retry the download
    } else {
      failed_names_holder[[i]] <<- c(accepted_name_v[[i]], e, format(Sys.time(), "%a %b %d %X %Y")) # Store other errors in our error handler as well. 
      failed_names_df <- as.data.frame(do.call(rbind, failed_names_holder))
      colnames(failed_names_df) <- c("acceptedParentName", "errorMsg", "timeStamp")
      fwrite(failed_names_df, "/home/jtmiller/my_elements/jtmiller/BoCP/data/idigbio-dws/raw-dw-info-tbls/failed_names.csv", append = TRUE, col.names = FALSE)
      return(NULL)
      print("NULL found")
    } # end of else of eror reporting
  }) # end of try catch
} # End of retry download fxn

# Loop through Raw Data Downloads from iDigBio by Names

for(i in 1:length(name_list)){
  time_taken <- bench::bench_time({
  occurrence_raw_dw <- retry_download(i,0) # see fxn above, retries in cases where API call fails
  })
  print(paste("trying", accepted_name_v[[i]]))
  if(!is.null(occurrence_raw_dw)){
    info_df_temp <- data.frame(
      acceptedParentName = accepted_name_v[[i]], # Store the acceptedNameParent 
      numOccurrences = nrow(occurrence_raw_dw), # Store information on the number of occurrence records in raw download
      timeTaken = time_taken[["real"]] # Store the time (wall time) taken
    )
    fwrite(info_df_temp, file = "/home/jtmiller/my_elements/jtmiller/BoCP/data/idigbio-dws/raw-dw-info-tbls/raw-info-tbl.csv", append = TRUE, col.names = !file.exists("/home/jtmiller/my_elements/jtmiller/BoCP/data/idigbio-dws/raw-dw-info-tbls/raw-info-tbl.csv"))
  
  }
  if(is.null(occurrence_raw_dw)){
    next # skip to the next iteration if the download failed or has 0 data. 
  }
  
  fwrite(occurrence_raw_dw, file = paste0("/home/jtmiller/my_elements/jtmiller/BoCP/data/idigbio-dws/raw/", accepted_name_filestyle_v[[i]], "-raw.csv"))
}


