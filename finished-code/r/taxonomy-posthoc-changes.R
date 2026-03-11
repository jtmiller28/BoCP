### Title: Taxonomy Post-Hoc Changs
### Author: JT Miller
### Date: 02-03/2026

## Load Packages
library(data.table)
library(dplyr)
library(readxl)

## Load in full taxonomy 
wcvp_ncbi_taxonomy <- fread("/blue/guralnick/millerjared/PlantSweepeR/data/processed/wcvp-ncbi-alignment-2025.csv")

## Load in Joey's adjusted taxonomic list
mich_adjusted_names <- read_excel("/blue/guralnick/millerjared/BoCP/data/processed/bocp_na_new_names_to_master_taxonomy.xlsx")

## Take name's from each, look at where they go, adjust for conflict 
names <- mich_adjusted_names$Name
# find names in wcvp_ncbi_taxonomy
wcvp_ncbi_taxonomy2 <- wcvp_ncbi_taxonomy %>% 
  filter(name %in% names)
