### Stage Occurrence Data for Modeling
### Author: JT Miller
### Date: 1-29-2026

## libraries
library(arrow)
library(dplyr)
library(data.table)
library(sf)
library(tidyverse)


## Open dataset
ds_parq <- open_dataset("/blue/guralnick/millerjared/PlantSweepeR/data/processed/sp-partitioned-occs-flagged_parquet/")

## read in level3 botanical regions
#level3_regions <- sf::read_sf("/blue/guralnick/millerjared/BoCP/data/raw/level3-wgsrpd/level3.shp")

## Retrieve Sp that have any occurrences in FL
FL_all_possible_taxa <- ds_parq |>
  filter(area_code_l3 == "FLA") |> 
  distinct(species, wcvpRangeStatus) |>
  group_by(species) |>
  summarize(FL_status = ifelse(any(wcvpRangeStatus == "native"), "native", "introduced or undocumented")) |>
  collect()

## Filter to only include native taxa
FL_natives <- FL_all_possible_taxa |>
  filter(FL_status == "native")

## Bring in species list for FL natives according to POWO
FL_2025_taxa <- fread("/blue/guralnick/millerjared/BoCP/data/processed/fl-names-2025-update.csv")
FL_2025_natives <- FL_2025_taxa %>% filter(FL_native == TRUE)

## Grab these species, compile datasets for models 
for(i in 1:length(unique(FL_natives$species))){
  sp_name <- unique(FL_natives$species)[i]
  sp_occs <- ds_parq |> 
    filter(species == sp_name) |> 
    collect()
  fwrite(sp_occs, paste0("/blue/guralnick/millerjared/BoCP/data/processed/FL-occs/", gsub(" ", "_", sp_name), ".csv"))
  print(paste("finished", sp_name, "taxa", i, "out of", length(unique(FL_natives$species))))
}

## Add in taxa suggesteed to be native to FL but dont have occs
FL_misses <- setdiff(FL_2025_natives$alignedParentName, FL_natives$species)


## Grab these species if they exist, compile datasets for models 
for(i in 1:length(unique(FL_misses))){
  sp_name <- unique(FL_misses)[i]
  sp_occs <- ds_parq |> 
    filter(species == sp_name) |> 
    collect()
  if(nrow(sp_occs) > 0){
  fwrite(sp_occs, paste0("/blue/guralnick/millerjared/BoCP/data/processed/FL-occs/", gsub(" ", "_", sp_name), ".csv"))
  } else{
    print(paste("no data for", sp_name))
  }
    print(paste("finished", sp_name, "taxa", i, "out of", length(unique(FL_natives$species))))
}

## Due to some nativity issues with my initial run, we're going to rerun a few species
fl_finished <- list.files("/blue/guralnick/millerjared/BoCP/data/processed/FL-occs/")
range_check_status <- fread("/blue/guralnick/millerjared/PlantSweepeR/data/processed/names-that-changed-range-status.csv")
names_to_redo <- intersect(fl_finished, range_check_status$species)
names_to_redo <- gsub(".csv", "", names_to_redo)
names_to_redo <- gsub("_", " ", names_to_redo)
for(i in 1:length(names_to_redo)){
  sp_name <- unique(names_to_redo)[i]
  sp_occs <- ds_parq |> 
    filter(species == sp_name) |> 
    collect()
  if(nrow(sp_occs) > 0){
    fwrite(sp_occs, paste0("/blue/guralnick/millerjared/BoCP/data/processed/FL-data-redo/", gsub(" ", "_", sp_name), ".csv"))
  } else{
    print(paste("no data for", sp_name))
  }
  print(paste("finished", sp_name, "taxa", i, "out of", length(names_to_redo)))
}

## Pacific NorthWest Taxa
pnw_genera <- fread("/blue/guralnick/millerjared/BoCP/data/processed/fam_gen_pnw.csv")
# find species in our tables
wcvp_na_table <- fread("/blue/guralnick/millerjared/PlantSweepeR/data/processed/wcvp-ncbi-alignment-na.csv")
wcvp_na_genera <- wcvp_na_table %>% 
  mutate(alignedGenus = word(alignedParentName, 1)) %>% 
  filter(alignedGenus %in% pnw_genera$Genus)

pnw_sp <- unique(wcvp_na_genera$alignedParentName)
range_check_status <- fread("/blue/guralnick/millerjared/PlantSweepeR/data/processed/names-that-changed-range-status.csv")
range_check_status <- range_check_status %>% 
  mutate(species = gsub(".csv", "", species)) %>% 
  mutate(species = gsub("_", " ", species))
names_to_redo <- intersect(pnw_sp, range_check_status$species)

## Grab these species if they exist, compile datasets for models 
for(i in 1:length(unique(pnw_sp))){
  sp_name <- unique(pnw_sp)[i]
  sp_occs <- ds_parq |> 
    filter(species == sp_name) |> 
    collect()
  if(nrow(sp_occs) > 0){
    fwrite(sp_occs, paste0("/blue/guralnick/millerjared/BoCP/data/processed/pnw-occs/", gsub(" ", "_", sp_name), ".csv"))
  } else{
    print(paste("no data for", sp_name))
  }
  print(paste("finished", sp_name, "taxa", i, "out of", length(unique(pnw_sp))))
}

## And redo names 
## Grab these species if they exist, compile datasets for models 
for(i in 1:length(unique(names_to_redo))){
  sp_name <- unique(names_to_redo)[i]
  sp_occs <- ds_parq |> 
    filter(species == sp_name) |> 
    collect()
  if(nrow(sp_occs) > 0){
    fwrite(sp_occs, paste0("/blue/guralnick/millerjared/BoCP/data/processed/pnw-redo-occs/", gsub(" ", "_", sp_name), ".csv"))
  } else{
    print(paste("no data for", sp_name))
  }
  print(paste("finished", sp_name, "taxa", i, "out of", length(unique(names_to_redo))))
}


## Grab Genome Size Taxa
genome_size_names <- fread("/blue/guralnick/millerjared/BoCP/data/processed/genome-na-names.txt", header = FALSE)
genome_size_names <- genome_size_names %>% 
  mutate(name = paste(V1, V2))
wcvp_taxonomy <- fread("/blue/guralnick/millerjared/PlantSweepeR/data/processed/wcvp-ncbi-alignment-na.csv")
wcvp_taxonomy <- wcvp_taxonomy %>% 
  select(name, alignedParentName)
genome_size_names <- genome_size_names %>% 
  left_join(wcvp_taxonomy, by = "name")

all_wcvp_taxonomy <- fread("/blue/guralnick/millerjared/PlantSweepeR/data/processed/wcvp-ncbi-alignment-2025.csv")
all_wcvp_taxonomy <- all_wcvp_taxonomy %>% 
  select(name, alignedParentName)
genome_size_names_notNA <- genome_size_names %>% 
  filter(is.na(alignedParentName)) %>% 
  left_join(all_wcvp_taxonomy, by = "name")
fwrite(genome_size_names_notNA, "/blue/guralnick/millerjared/BoCP/data/processed/cvals-names-not-NA.csv")

genome_size_names_na <- genome_size_names %>% 
  filter(!is.na(alignedParentName)) %>% 
  rename(species = alignedParentName)


## Remove names that have been covered in previous models 
genome_names_fin <- genome_size_names_na %>% 
  filter(species %in% pnw_sp | species %in% FL_2025_natives$alignedParentName)

genome_names_nonfin <- genome_size_names_na %>% 
  filter(!species %in% genome_names_fin$species)

## Grab these species if they exist, compile datasets for models for those native to NA
for(i in 1:length(unique(genome_names_nonfin$species))){
  sp_name <- unique(genome_names_nonfin$species)[i]
  sp_occs <- ds_parq |> 
    filter(species == sp_name) |> 
    collect()
  if(any(sp_occs$wcvpRangeStatus == "native")){
  
  if(nrow(sp_occs) > 0){
    fwrite(sp_occs, paste0("/blue/guralnick/millerjared/BoCP/data/processed/cval-occs/", gsub(" ", "_", sp_name), ".csv"))
  } else{
    print(paste("no data for", sp_name))
  }
  } else{
    print(paste(sp_name, "Has No Native Occs in NA"))
  }
  print(paste("finished", sp_name, "taxa", i, "out of", length(unique(genome_names_nonfin$species))))
}
names_w_data <- list.files("/blue/guralnick/millerjared/BoCP/data/processed/cval-occs/")
names_w_data <- gsub(".csv", "", names_w_data)
names_w_data <- gsub("_", " ", names_w_data)
likely_not_na <- setdiff(unique(genome_names_nonfin$species), names_w_data)


### Stage data for Deserts (Mojave and Sonoran To start)
desert_taxa <- fread("/blue/soltis/millerjared/desert-diversity/data/processed/desert_taxon_taxa_table.csv")

mojave_sonoran_taxa <- desert_taxa %>% 
  filter(desert %in% c("Mojave Basin and Range", "Sonoran Desert")) %>% 
  filter(wcvpRangeStatus == "native") %>% 
  distinct(species)

# check if any of these species have been modeled
cval_occs <- list.files("/blue/guralnick/millerjared/BoCP/data/processed/cval-occs/")
fl_occs <- list.files("/blue/guralnick/millerjared/BoCP/data/processed/FL-occs/")
fl_redo_occs <- list.files("/blue/guralnick/millerjared/BoCP/data/processed/FL-data-redo/")
pnw_occs <- list.files("/blue/guralnick/millerjared/BoCP/data/processed/pnw-occs/")
pnw_redo_occs <- list.files("/blue/guralnick/millerjared/BoCP/data/processed/pnw-redo-occs/")
all_sp_done_so_far <- c(cval_occs, fl_occs, fl_redo_occs, pnw_occs, pnw_redo_occs)
all_sp_done_so_far <- gsub(".csv", "", all_sp_done_so_far)
all_sp_done_so_far <- gsub("_", " ", all_sp_done_so_far)

mojave_sonoran_all_ready_fin <- mojave_sonoran_taxa %>% 
  filter(species %in% all_sp_done_so_far)

mojave_sonoran_to_model <- mojave_sonoran_taxa %>% 
  filter(!species %in% all_sp_done_so_far)

for(i in 1:length(unique(mojave_sonoran_to_model$species))){
  sp_name <- unique(mojave_sonoran_to_model$species)[i]
  sp_occs <- ds_parq |> 
    filter(species == sp_name) |> 
    collect()
  if(any(sp_occs$wcvpRangeStatus == "native")){
    
    if(nrow(sp_occs) > 0){
      fwrite(sp_occs, paste0("/blue/guralnick/millerjared/BoCP/data/processed/moj-son-occs/", gsub(" ", "_", sp_name), ".csv"))
    } else{
      print(paste("no data for", sp_name))
    }
  } else{
    print(paste(sp_name, "Has No Native Occs in NA"))
  }
  print(paste("finished", sp_name, "taxa", i, "out of", length(unique(mojave_sonoran_to_model$species))))
}

## Next chkpt 
# check if any of these species have been modeled
cval_occs <- list.files("/blue/guralnick/millerjared/BoCP/data/processed/cval-occs/")
fl_occs <- list.files("/blue/guralnick/millerjared/BoCP/data/processed/FL-occs/")
fl_redo_occs <- list.files("/blue/guralnick/millerjared/BoCP/data/processed/FL-data-redo/")
pnw_occs <- list.files("/blue/guralnick/millerjared/BoCP/data/processed/pnw-occs/")
pnw_redo_occs <- list.files("/blue/guralnick/millerjared/BoCP/data/processed/pnw-redo-occs/")
moj_son_occs <- list.files("/blue/guralnick/millerjared/BoCP/data/processed/moj-son-occs/")

all_sp_done_so_far <- c(cval_occs, fl_occs, fl_redo_occs, pnw_occs, pnw_redo_occs, moj_son_occs)
all_sp_done_so_far <- gsub(".csv", "", all_sp_done_so_far)
all_sp_done_so_far <- gsub("_", " ", all_sp_done_so_far)

## Load in Californias taxa 
ca_natives <- fread("/blue/guralnick/millerjared/BoCP/data/processed/CA_names_w_wcvpStatus.csv")
ca_natives <- ca_natives %>% filter(CA_species == TRUE & CA_native == TRUE)
unfin_ca_names <- setdiff(ca_natives$alignedParentName, all_sp_done_so_far)

for(i in 1:length(unique(unfin_ca_names))){
  sp_name <- unique(unfin_ca_names)[i]
  sp_occs <- ds_parq |> 
    filter(species == sp_name) |> 
    collect()
  if(any(sp_occs$wcvpRangeStatus == "native")){
    
    if(nrow(sp_occs) > 0){
      fwrite(sp_occs, paste0("/blue/guralnick/millerjared/BoCP/data/processed/ca-occs/", gsub(" ", "_", sp_name), ".csv"))
    } else{
      print(paste("no data for", sp_name))
    }
  } else{
    print(paste(sp_name, "Has No Native Occs in NA"))
  }
  print(paste("finished", sp_name, "taxa", i, "out of", length(unique(unfin_ca_names))))
}

## Next chkpt 
# check if any of these species have been modeled
cval_occs <- list.files("/blue/guralnick/millerjared/BoCP/data/processed/cval-occs/")
fl_occs <- list.files("/blue/guralnick/millerjared/BoCP/data/processed/FL-occs/")
fl_redo_occs <- list.files("/blue/guralnick/millerjared/BoCP/data/processed/FL-data-redo/")
pnw_occs <- list.files("/blue/guralnick/millerjared/BoCP/data/processed/pnw-occs/")
pnw_redo_occs <- list.files("/blue/guralnick/millerjared/BoCP/data/processed/pnw-redo-occs/")
moj_son_occs <- list.files("/blue/guralnick/millerjared/BoCP/data/processed/moj-son-occs/")
ca_occs <- list.files("/blue/guralnick/millerjared/BoCP/data/processed/ca-occs/")
all_sp_done_so_far <- c(cval_occs, fl_occs, fl_redo_occs, pnw_occs, pnw_redo_occs, moj_son_occs, ca_occs)
all_sp_done_so_far <- gsub(".csv", "", all_sp_done_so_far)
all_sp_done_so_far <- gsub("_", " ", all_sp_done_so_far)

# all NA
NA_taxa <- fread("/blue/guralnick/millerjared/BoCP/data/processed/NA-names-2025-update.csv")
NA_natives <- NA_taxa %>% filter(NA_native == TRUE)
NA_native_sp <- unique(NA_natives$alignedParentName)

NA_natives_unfin <- setdiff(NA_native_sp, all_sp_done_so_far)

for(i in 1:length(unique(NA_natives_unfin))){
  sp_name <- unique(NA_natives_unfin)[i]
  sp_occs <- ds_parq |> 
    filter(species == sp_name) |> 
    collect()
  if(any(sp_occs$wcvpRangeStatus == "native")){
    
    if(nrow(sp_occs) > 0){
      fwrite(sp_occs, paste0("/blue/guralnick/millerjared/BoCP/data/processed/NA-occs/", gsub(" ", "_", sp_name), ".csv"))
    } else{
      print(paste("no data for", sp_name))
    }
  } else{
    print(paste(sp_name, "Has No Native Occs in NA"))
  }
  print(paste("finished", sp_name, "taxa", i, "out of", length(unique(NA_natives_unfin))))
}

## Extra Stuff for other proj
overlap_table <- fread("/blue/soltis/millerjared/SpatioTemporalTradeoffs/data/processed/legume-spatial-overlaps.csv")
wcvp_na_table <- fread("/blue/guralnick/millerjared/PlantSweepeR/data/processed/wcvp-ncbi-alignment-na.csv")

plant_sp <- overlap_table %>% distinct(plant_sp) %>% rename(name = plant_sp) %>% left_join(wcvp_na_table, by = "name")

plant_sp_aligned_name <- unique(plant_sp$alignedParentName)

for(i in 1:length(plant_sp_aligned_name)){
  sp_name <- unique(plant_sp_aligned_name)[i]
  sp_occs <- ds_parq |> 
    filter(species == sp_name) |> 
    collect()
  if(any(sp_occs$wcvpRangeStatus == "native")){
    
    if(nrow(sp_occs) > 0){
      fwrite(sp_occs, paste0("/blue/soltis/millerjared/Boyd-Deep-Models/SDM_pipeline/Legumes-data/sp-occs/", gsub(" ", "_", sp_name), ".csv"))
    } else{
      print(paste("no data for", sp_name))
    }
  } else{
    print(paste(sp_name, "Has No Native Occs in NA"))
  }
  print(paste("finished", sp_name, "taxa", i, "out of", length(unique(plant_sp_aligned_name))))
}
