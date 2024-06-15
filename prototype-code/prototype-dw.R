install.packages("gatoRs")

# Load packages
library(gatoRs)
library(data.table)
library(dplyr)

# Read in the name alignment
name_alignment <- fread("/home/jtmiller/my_elements/jtmiller/BoCP/data/processed/finalized_name_alignment.csv")

# Structure these names in a list with acceptedName + synonym(s)



# Testing gatoRs changes 
bench::bench_time({
gators_unedited <- gatoRs::gators_download("Olneya tesota")
})
bench::bench_time({
gators_edited <- gators_download_edited("Olneya tesota", call.idigbio = TRUE, call.gbif = FALSE)
})