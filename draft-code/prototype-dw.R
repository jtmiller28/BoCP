# Prototype Download Script for iDigBio occurrence records
## Goals:
# 1. Pull idigbio occurrence data for taxa of interest
# 2. Filter down data to ensure taxonomy is congruent 
# 3. Create appropriate error handlers to deal with automation

# Load packages
library(gatoRs)
library(data.table)
library(dplyr)
library(rgnparser)
source("prototype-code/gatoRs-fxns-edited.R") # edited gatoRs fxns to allow for only pulling idigbio
# Some path shenanigans, as gnparser and rstudioserver do not communicate well. Basically grab our default path and append the homebrew path I made, and read that in as the System Defualt variable. 
my_path <- Sys.getenv("PATH") # grab our path
Sys.setenv(PATH = paste0(my_path, ":/home/linuxbrew/.linuxbrew/Cellar/gnparser/1.10.1/bin"))

# Read in the name alignment
name_alignment <- fread("/home/jtmiller/my_elements/jtmiller/BoCP/data/processed/finalized_name_alignment_wcvp.csv")

# Structure these names in a list with acceptedName + synonym(s)
accepted_name_v <- unique(name_alignment$acceptedName)
# Make a file style version of this vector for storage purpsoses
accepted_name_filestyle_v <- gsub(" ", "-", accepted_name_v)
name_list <- list() # initialize empty list

bench::bench_time({
for(i in 1:length(accepted_name_v)){
  p <- name_alignment[acceptedName == accepted_name_v[[i]]] # create a for-loop that goes through by the accepted name and grabs synonyms, storing them both as a vector in a list. 
  s <- p[!is.na(synonyms)]
  x <- print(s$synonyms)
  a <- print(unique(p$acceptedName))
  name_list[[i]] <- unique(c(a,x))
}
})

for(i in 1:length(name_list)){
# Proceed to download
bench::bench_time({
occurrence_raw_dw <- gators_download_edited(name_list[[648]], call.idigbio = TRUE, call.gbif = FALSE)
# Name Clean 
## Filter for name alignment
occurrence_raw_dw <- as.data.table(occurrence_raw_dw)
occurrence_parsed <- gn_parse_tidy(occurrence_raw_dw$scientificName)  # break up names
canonicalFull <- occurrence_parsed[,6] # the canonicalfull name broken up during parsing
authorship <- occurrence_parsed[,7] # the authorship name broken up during parsing 
occurrence_df <- cbind(occurrence_raw_dw, canonicalFull, authorship)
occurrence_tax1_cl <- occurrence_df[canonicalfull %in% name_list[[648]]] # assure that what we are downloading for is the same
## Check multiple mapping status, if TRUE construct require aligned authorship for cleaning 
name_alignment_snip <- name_alignment[acceptedName == accepted_name_v[648]]
name_alignment_accepted_name_info <- name_alignment_snip[1, ] # grab any case of the accepted name
  if(name_alignment_accepted_name_info$acceptedNameMultMapsPossible == TRUE){ # if accepted name does have multiple mappings based on authorship possible
    if(name_alignment_accepted_name_info$acceptedNameAuthorshipMatch == TRUE & name_alignment_accepted_name_info$acceptedNameMultMapResolutionPossible == TRUE){ # and if resolution of these was possible during name alignment
        occurrence_tax1_cl <- occurrence_tax1_cl[, spacelessAuthorship := gsub(" ", "", authorship)]
        accepted_name_occ_df_subset <- occurrence_tax1_cl[canonicalfull == accepted_name_v[648]]
        synonym_names_occ_df_subset <- occurrence_tax1_cl[canonicalfull != accepted_name_v[648]]
        # require that the authorship be present and exactly match. 
        accepted_name_occ_df_subset_c <- accepted_name_occ_df_subset[tolower(spacelessAuthorship) == tolower(name_alignment_accepted_name_info$acceptedNameSpacelessWCVPAuthorship)]
        # bring df back together 
        occurrence_tax2_cl <- rbind(accepted_name_occ_df_subset_c, synonym_names_occ_df_subset)
        accepted_name_mult_map_loss <- nrow(accepted_name_occ_df_subset) - nrow(accepted_name_occ_df_subset_c) # for storing info
    } else {
      # dead case, alignment would not proceed if authorship match was not TRUE when multmaps is TRUE
      accepted_name_mult_map_loss <- 0 
    }
    
  } else { 
    occurrence_tax2_cl <- occurrence_tax1_cl # if acceptedNameMultMapsPossible is FALSE, then we can proceed as normal 
    
  }
  if(any(name_alignment_snip$synonymMultMapPossible == TRUE)){ # if multiple mappings is ever possible in synonyms, proceed to check authorship for those names
    # Extract cases where this is true
    syn_mult_map_cases <- name_alignment_snip[synonymMultMapPossible == TRUE]
    synonym_names_occ_df_subset <- occurrence_tax2_cl[canonicalfull %in% syn_mult_map_cases$synonyms]
    accepted_names_occ_df_subset <- occurrence_tax2_cl[canonicalfull == accepted_name_v[648]]
    num_of_checks <- nrow(syn_mult_map_cases)
    holder <- data.frame()
    for(i in seq(num_of_checks)){
      synonym_names_occ_df_subset_i <- synonym_names_occ_df_subset[canonicalfull %in% syn_mult_map_cases$synonyms[2]]
      synonym_alignment_info <- syn_mult_map_cases[648] # extract the case for authorship
      if(nrow(synonym_names_occ_df_subset_i) > 0){ # if anything even exists for said synonym
        synonym_names_occ_df_subset_ic <- synonym_names_occ_df_subset_i[tolower(spacelessAuthorship) == tolower(synonym_alignment_info$synonymNameWCVPAuthorship)]
        holder <- rbind(synonym_names_occ_df_subset_ic) # bind the data. 
      } else {
       # if there is no data for a synonym then do nothing
      }
      
    }
    occurrence_tax3_cl <- rbind(accepted_names_occ_df_subset, holder) # rebuild
    synonym_name_mult_map_loss_info <- nrow(occurrence_tax2_cl) - nrow(occurrence_tax3_cl) # for storing info 
  } else {
    occurrence_tax3_cl <- occurrence_tax2_cl # if no synonyms are needed to be remapped, proceed without 
    synonym_name_mult_map_loss_info <- 0 # store info 
  }
  
# Initiate cleaning for records without georeference (lat lon)

# Build a dataframe noting the total raw records, the cleaned taxonomic records (change to accomodate rbind structure...)
info_df <- data.frame(acceptedName = as.character(accepted_name_v[[648]]), 
                      idigbioRawRecordCount = as.integer(nrow(occurrence_raw_dw)),
                      taxonomicAligmentClean = as.integer(nrow(occurrence_tax1_cl)),
                      taxonomicAcceptedNameMultClean = as.integer(nrow(occurrence_tax2_cl)),
                      taxonomicSynonymNameMultClean = as.integer(nrow(occurrence_tax3_cl)))


# Write out data into CSVs 
fwrite(occurrence_raw_dw, paste("/home/jtmiller/my_elements/jtmiller/BoCP/data/idigbio-dws/raw/", accepted_name_filestyle_v[[1]]))

})
}


