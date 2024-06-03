# The most basic way to grab counts for each name
# note that these will grab more than we desire but 
# will reflect the init data download on the gbif side

### Load Libraries 
library(rgbif)
library(data.table)
library(tidyverse)
library(taxadb)
### Load names  to count 
names <- fread("/home/millerjared/blue_guralnick/millerjared/BoCP/data/processed/finalized_name_alignment.csv")
gbif_name_match <- names %>% 
  mutate(gbif_id = get_ids(acceptedName, "gbif")) %>% 
  mutate(gbif_harmonizedName = get_names(gbif_id, "gbif"))

gbif_match_fails <- gbif_name_match %>% 
  filter(is.na(gbif_harmonizedName))

gbif_match_s_retry <- gbif_match_fails %>% 
  mutate(gbif_id2 = get_ids(synonym, "gbif")) %>% 
  mutate(gbif_harmonizedName2 = get_names(gbif_id2, "gbif"))

gbif_match_s_retry_f <- gbif_match_s_retry %>% 
  group_by(acceptedName) %>% 
  mutate(foundMatch = case_when(
    any(!is.na(gbif_harmonizedName2)) ~ TRUE,
    TRUE ~ FALSE
  ))

gbif_match_s_retry_d <- gbif_match_s_retry_f %>% 
  filter(!is.na(gbif_harmonizedName2) & foundMatch == TRUE) %>% 
  distinct(acceptedName, .keep_all = TRUE)

unresolvables <- gbif_match_s_retry_f %>%  
  filter(is.na(gbif_harmonizedName2) & foundMatch == FALSE) %>% 
  distinct(acceptedName, .keep_all = TRUE)

# Grab the successful matches 
gbif_match_successes <- gbif_name_match %>% 
  filter(!is.na(gbif_harmonizedName))
gbif_match_retry_successes <- gbif_match_s_retry_d

harmonized_names <- unique(c(gbif_match_successes$gbif_harmonizedName, gbif_match_retry_successes$gbif_harmonizedName2))
### Run a simple count estimate 
result_table <- data.frame(harmonizedName = character(), totalCount = numeric(), stringsAsFactors = FALSE) # create a storage df
harmonized_names_filestyle <- gsub(" ", "-", harmonized_names)
for(i in 24700:length(harmonized_names)){ # change indexing to 1:length(harmonized_names)
  name_to_count <- harmonized_names[[i]] # pull a name 
  
  name_suggest_df <- rgbif::name_suggest(name_to_count) # call the api for suggested names associated with this parent name
  count_storage <- list() # create a storage list 
  
  for(j in 1:length(name_suggest_df$data$key)){ # loop through each associated name within name_suggest_df 
    name_count <- rgbif::occ_count(taxonKey = name_suggest_df$data$key[j])
    provided_name <- name_to_count # to associate the name we'll call this by 
    name_suggest <- name_suggest_df$data$canonicalName[1]
    canonical_name <- name_suggest_df$data$canonicalName[j]
    if(!is.null(name_suggest)){
    count_storage[[j]] <- list(providedName = provided_name,
                                  canonicalName = canonical_name,
                                  gbifNameSuggest = name_suggest, 
                                  nameCount = name_count)
    
    }
  } # end of j loop
  
  # Do some backwards manipulation to find the number of records associated with the present name 
  count_df <- rbindlist(count_storage)
  if(nrow(count_df) > 0){ # conditional to ensure that there is information. 
    count_df <- count_df %>% # create two new fields, one that notes the total of occurrences and one that notes the adjusted counts 
      mutate(adjustNameCount = case_when(
        row_number() == 1 & length(unique(count_df$canonicalName)) == 1 ~ nameCount,
        row_number() == 1 ~ nameCount - sum(count_df$nameCount[2:nrow(count_df)]),
        row_number() != 1 ~ nameCount
      )) %>% 
      mutate(totalNameCount = max(count_df$nameCount)) # give the max count assoicated with each name 
      
    # write out in a directory for storage, 
    fwrite(count_df, paste("/home/millerjared/blue_guralnick/millerjared/BoCP/data/processed/gbif-count-tables-full/",harmonized_names_filestyle[[i]], "counts.csv" ))
    # write out a df for totals per representative name
    
  }
result_table <- rbind(result_table, data.frame(harmonizedName = harmonized_names[[i]], totalCount = max(count_df$nameCount)))
print(paste("finished name ", i, "out of ", max(length(harmonized_names))))
} # end of i loop

fwrite(result_table, "/home/millerjared/blue_guralnick/millerjared/BoCP/data/processed/final_counts_gbif5.csv")


# Append the several sessions worth of data downloads together 
final_counts1 <- fread("/blue/guralnick/millerjared/BoCP/data/processed/final_counts_gbif1.csv")
final_counts2 <- fread("/blue/guralnick/millerjared/BoCP/data/processed/final_counts_gbif2.csv")
final_counts3 <- fread("/blue/guralnick/millerjared/BoCP/data/processed/final_counts_gbif3.csv")
final_counts4 <- fread("/blue/guralnick/millerjared/BoCP/data/processed/final_counts_gbif4.csv")
final_counts5 <- fread("/blue/guralnick/millerjared/BoCP/data/processed/final_counts_gbif5.csv")
final_counts6 <- fread("/blue/guralnick/millerjared/BoCP/data/processed/final_counts_gbif6.csv")

final_counts_full <- rbind(final_counts1, final_counts2, final_counts3, final_counts4, final_counts5, final_counts6)
final_counts_full <- distinct(final_counts_full, harmonizedName, totalCount, .keep_all = TRUE)
final_counts_full$totalCount <- ifelse(final_counts_full$totalCount == "-Inf", 0, final_counts_full$totalCount)
fwrite(final_counts_full, "/blue/guralnick/millerjared/BoCP/data/processed/final-rough-gbif-counts.csv")
