# 001E Processed GBIF data parquet formatting
# Author: JT Miller 
# Date: 03-06-2024
# Project: BoCP 

# Load packages 
library(tidyverse)
library(data.table)
library(arrow)
library(DBI)
library(duckdb)

# read in the verbatim gbif cut up snapshots
processed1 <- fread("/home/millerjared/blue_guralnick/millerjared/BoCP/data/raw/gbif-tracheophyte-dw/select-occurrence1.txt")
processed2 <- fread("/home/millerjared/blue_guralnick/millerjared/BoCP/data/raw/gbif-tracheophyte-dw/select-occurrence2.txt")
processed3 <- fread("/home/millerjared/blue_guralnick/millerjared/BoCP/data/raw/gbif-tracheophyte-dw/select-occurrence3.txt")
processed4 <- fread("/home/millerjared/blue_guralnick/millerjared/BoCP/data/raw/gbif-tracheophyte-dw/select-occurrence4.txt")

processed <- processed1 %>% 
  rbind(processed2) %>% 
  rbind(processed3) %>% 
  rbind(processed4)

write_parquet(processed, "/home/millerjared/blue_guralnick/millerjared/BoCP/data/raw/parquet-occs/gbif_processed_occ.parquet")