# 001D Symbiota data parqueting
# Author: JT Miller 
# Date: 03-06-2024
# Project: BoCP 

# Load packages 
library(tidyverse)
library(data.table)
library(arrow)
library(DBI)
library(duckdb)

### Load symbiota data

# read in data on vascular plants obtained from Ed Gilbert at Symbiota. Append colnames as those didnt quite make it on the data
colnames_df <- fread("/blue/guralnick/millerjared/BoCP/data/symbiota-dws/symbiota/field_heading.csv")
new_cols <- names(colnames_df)
cch2 <- fread("/blue/guralnick/millerjared/BoCP/data/symbiota-dws/symbiota/cch2_export.csv")
colnames(cch2) <- new_cols
cnh <- fread("/blue/guralnick/millerjared/BoCP/data/symbiota-dws/symbiota/cnh_export.csv")
colnames(cnh) <- new_cols
sienet <- fread("/blue/guralnick/millerjared/BoCP/data/symbiota-dws/symbiota/seinet_export.csv")
colnames(sienet) <- new_cols

# combine data
symbiota_data <- rbind(cch2, cnh, sienet)

# restructure NAs (currently are \Ns)
symbiota_data[symbiota_data == "\\N"] <- NA