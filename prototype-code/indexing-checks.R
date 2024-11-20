## Testing script to see where we're at for indexing names on jobs 
library(data.table)
library(dplyr)
names_completed <- fread("/home/millerjared/blue_guralnick/millerjared/BoCP/data/gbif-dws/raw-dw-info-tbls/raw-info-tbl16.csv")
failed_names_checked <- fread("/home/millerjared/blue_guralnick/millerjared/BoCP/data/gbif-dws/raw-dw-info-tbls/failed_names16.csv", header = FALSE)
# Read in the name alignment
name_alignment <- fread("./data/processed/finalized_name_alignment_wcvp.csv")
accepted_name_v <- unique(name_alignment$acceptedNameParent)
### Subset to do only a portion of names for this script. 

accepted_name_v <- accepted_name_v[30501:32867]

# append failed names together with finished
finished_names <- c(names_completed$acceptedParentName, unique(failed_names_checked$V1))

check <- setdiff(accepted_name_v, finished_names)

grep("Lemna minuta", accepted_name_v)
