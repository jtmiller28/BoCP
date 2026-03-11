### Model Computation Time Estimates
### Project: BoCP
### Author: JT Miller
### Date: 04/03/2025

### taking the number of 5km grid cells a species occupies, we're going to decide how to sample taxa

# load libraries
library(data.table)
library(tidyverse)

# load data
sp_occ <- fread("/blue/guralnick/millerjared/BoCP/outputs/5km_cells_occupied.csv")

# clean up header issues, make class correct
sp_occ <- sp_occ %>% filter(species != "species", num_5km_cells != "num_5km_cells") %>% mutate(num_5km_cells = as.numeric(num_5km_cells))

# remove 0s as these wont be modeled regardless
sp_model_occs <- sp_occ %>% 
  filter(num_5km_cells > 10)
# build a simple historgram for num of 5km cells
ggplot(sp_model_occs , mapping = aes(x = num_5km_cells)) + 
  geom_histogram(aes(y = ..density..), binwidth = 200, fill = "steelblue", color = "black") + 
  geom_density(alpha = 0.5, fill = "darkred") #+ 
  #xlim(0, 100000)

# now bin into quantiles of data quantity. 
sp_model_occs$num_cell_percentile <- cut(sp_model_occs$num_5km_cells, 
                                         breaks = quantile(sp_model_occs$num_5km_cells, probs = seq(0,1,0.25), na.rm = TRUE), 
                                         include.lowest = TRUE, 
                                         labels = c("0-25th", "26-50th", "51-75th", "76-100th"))

# list of species that my labmates and I thought we'd like to see
# Liatris ohlingerae
# Olneya tesota 
# Rafinesquia neomexicana
# Calochortus albus
# Crinum americanum
set.seed(999)
selected_0_25 <- sp_model_occs %>% filter(num_cell_percentile == "0-25th") %>% slice_sample(n = 25)
selected_26_50 <- sp_model_occs %>% filter(num_cell_percentile == "26-50th") %>% slice_sample(n = 24)
selected_51_75 <- sp_model_occs %>% filter(num_cell_percentile == "51-75th") %>% slice_sample(n = 25)
selected_76_100 <- sp_model_occs %>% filter(num_cell_percentile == "76-100th") %>% slice_sample(n = 21)
selected_manual <- sp_model_occs %>% filter(species %in% c("Liatris ohlingerae", "Olneya tesota", "Rafinesquia neomexicana", 
                                                           "Calochortus albus", "Crinum americanum"))

selected_sp <- bind_rows(selected_0_25, selected_26_50, selected_51_75, selected_76_100, selected_manual)
fwrite(selected_sp, "./outputs/select-sp-test-table.csv")
# grab these selected species files, do our full spatial flagging clean, then write to a trial dir
# grab names 
sp_vector <- selected_sp$species
sp_vector <- gsub(" ", "-", sp_vector)

# for loop read in names 
for(i in 1:length(sp_vector)){
  # read in the species table 
  print(paste("retrieving data for", sp_vector[i]))
  
  test_df <- fread(paste0("/blue/guralnick/millerjared/BoCP/data/processed/fully-flagged-data/", sp_vector[i], ".csv"))
  test_df_filtered <- test_df %>% 
    filter(taxonomicExactMatch == TRUE) %>% 
    filter(wgs84Datum == TRUE) %>% 
    filter(coordinateIssue == FALSE) %>% 
    filter(trueCoordsWithheld == FALSE) %>% 
    filter(wcvpRangeStatus == "native" | wcvpRangeStatus == "introduced") %>% 
    filter(validRecord == TRUE) %>% 
    filter(equalLatLon == FALSE) %>% 
    filter(zeroCoords == FALSE) %>% 
    filter(capitalCoord == FALSE) %>% 
    filter(centroidCoord == FALSE) %>% 
    filter(inOceanCoord == FALSE) %>% 
    filter(inGBIFHeadquarters == FALSE) %>% 
    filter(inInstitutionBounds == FALSE)
  
  
  test_df_filtered <- test_df_filtered %>%
    filter(!amongAggDuplicate | is.na(AggDuplicateGroupID)) %>% # Keep non-duplicates
    bind_rows( # Add in the lowest-ranked duplicates
      test_df_filtered %>%
        filter(amongAggDuplicate) %>%
        group_by(AggDuplicateGroupID) %>%
        filter(AggDuplicateRank == min(AggDuplicateRank, na.rm = TRUE)) %>%
        slice(1) %>%  # if tied for lowest rank, just take the first
        ungroup()
    )
  
  test_df_filtered <- test_df_filtered %>%
    filter(!specimenDuplicate | is.na(specimenDuplicateGroupID)) %>% # Keep non-duplicates
    bind_rows( # Add in the lowest-ranked duplicates
      test_df_filtered %>%
        filter(specimenDuplicate) %>%
        group_by(specimenDuplicateGroupID) %>%
        filter(specimenDuplicateRank == min(specimenDuplicateRank, na.rm = TRUE)) %>%
        slice(1) %>%  # if tied for lowest rank, just take the first
        ungroup()
    )
  fwrite(test_df_filtered, paste0("/blue/guralnick/millerjared/BoCP/data/processed/test-sdm-prep/", sp_vector[i], ".csv"))
}