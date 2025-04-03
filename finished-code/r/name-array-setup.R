### A script for creating a temporary object noting what names need to be finished in a directory

# read the taxonomic harmonized list of names + synonyms
name_list <- readRDS("/blue/guralnick/millerjared/BoCP/data/processed/name_list.rds")

# for the last flagging step. 
# check for names that have occurrence data
list_found_names <- list.files("/blue/guralnick/millerjared/BoCP/data/processed/dupe-flagged/")
list_found_names <- gsub("-", " ", list_found_names)
list_found_names <- gsub(".csv", "", list_found_names)

# check current contents 
list_done_names <- list.files("/blue/guralnick/millerjared/BoCP/data/processed/fully-flagged-data/")
list_done_names <- gsub("-", " ", list_done_names)
list_done_names <- gsub(".csv", "", list_done_names)

# find the difference
unfinished_names <- setdiff(list_found_names, list_done_names)
write_rds(unfinished_names, "/blue/guralnick/millerjared/BoCP/data/processed/unfinished_flag_names_temp.rds")
