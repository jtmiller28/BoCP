## just some scratch code to produce tables of sp summaries
library(data.table)
library(tidyverse)

species_files <- list.files("/blue/guralnick/millerjared/BoCP/data/processed/species-occs/")
accepted_names <- gsub("-", " ", species_files)
accepted_names <- gsub(".csv", "", accepted_names)
for(i in 1:length(species_files)){
  occ_data <- fread(paste0("/blue/guralnick/millerjared/BoCP/data/processed/species-occs/", species_files[i]))
  occ_data <- occ_data %>% 
    summarize(accepted_names[i], n = n())
}