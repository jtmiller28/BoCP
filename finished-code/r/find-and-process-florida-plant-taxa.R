### Title: Find and Process Florida Plant Taxa
### Author: JT Miller

library(arrow) # for processing parquet files
library(dplyr) # normal data wrangling 
library(ggplot2) # visuals

# use arrow to open a connection to the parquet files (does not load data into mem)
ds_na_occs <- arrow::open_dataset("/blue/soltis/share/sp-partitioned-occs-flagged_parquet/")

# these are the tdwg level3 botanical regions shapefiles. If you want florida, just search for FLA in the lvl 3 Code fields (requries the shapefiles, you can download these if needed: https://github.com/tdwg/wgsrpd)
#level3_regions <- sf::read_sf("/blue/guralnick/millerjared/BoCP/data/raw/level3-wgsrpd/level3.shp")
#ggplot() + geom_sf(level3_regions, mapping = aes()) + geom_sf(filter(level3_regions, LEVEL3_COD == "FLA"), mapping = aes(), fill = "darkgreen")

# We can use arrow and dplyr to preform wrangling tasks prior to the data being read into memory (collect) 
## When building the sp datasets, I already gave each occurrence a label of what level3 region their point resides in (if georef), 
## with this we can generate a list of species for Florida. Note that LEVEL3_COD = area_code_l3
system.time({ # reports time it takes 
  FL_all_possible_taxa <- ds_na_occs |> 
    filter(area_code_l3 == "FLA") |> # filter to only data within FL
    #mutate(wcvpRangeStatus = as.character(wcvpRangeStatus)) |> 
    distinct(species, wcvpRangeStatus) |>  # find distinct species and their range status (native, introduced, or undocumented)
    group_by(species) |> # group by species
    summarize(FL_status = ifelse(any(wcvpRangeStatus == "native"), "native", "introduced or undocumented")) |> # delimits native, or if its introduced/undocumented
    collect() # reads into memory
})
