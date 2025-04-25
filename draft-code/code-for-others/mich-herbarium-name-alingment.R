# Michigian Herbarium Name Alignment 
# Author: JT Miller 
# Date: 03-14-2025
# Project: BoCP 

# Load libraries
library(data.table)
library(rgnparser)
library(tidyverse)
library(taxadb)

#### Set up pathing to use gnparser
my_path <- Sys.getenv("PATH") # grab our path
Sys.setenv(PATH = paste0(my_path, "/home/millerjared/gnparser"))

# Load Michigian Names
mich_names <- fread("/blue/guralnick/millerjared/BoCP/data/raw/mich-raw-names.csv")

# Select to simplify 
mich_names <- mich_names %>% select(name = sciname) %>% distinct()

# Load name alignment 
name_alignment <- fread("/blue/guralnick/millerjared/BoCP/data/processed/wcvp-ncbi-alignment-na.csv")

# Align Names 
mich_names_aligned <- mich_names %>% 
  left_join(name_alignment, by = "name") %>% 
  distinct(name, ncbiPossibleIDs, .keep_all = TRUE) %>%  # we gain a few duplicates during the join due to there being multiple values of ncbiAcceptedNameUsageID
  mutate(alignmentFound = ifelse(is.na(alignedName), FALSE, TRUE))

mich_names_aligned %>% group_by(alignmentFound) %>% summarize(n = n())

fwrite(mich_names_aligned, "/blue/guralnick/millerjared/BoCP/outputs/mich-aligned-names.csv")
