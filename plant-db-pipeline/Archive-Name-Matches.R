### Archive Known Name Matches in Major DBs (Archive All Possible Matches)
### Author: JT Miller
### 08/28/2025

# Purpose: Build Fuzzy Name Relation Table to later use in Exact Match Schema to speed up building species flagged tables 

## Load Libraries
library(duckdb)
library(DBI)
library(dplyr)
library(data.table)

## get SLURM info
start_num <- as.numeric(Sys.getenv("START_NUM"))
array_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
cat("Running array task:", array_id, "\n")
task_id <- as.numeric(start_num)
# Set up connection with table
con <- dbConnect(duckdb::duckdb(), 
                 dbdir = "/blue/guralnick/millerjared/BoCP/plant-db-pipeline/data/archive_name_relations.duckdb", 
                 read_only = TRUE) # make sure this is TRUE, write causes errors

## For each accepted_name, match all patterns that are built by the name variants
name_list <- readRDS("/blue/guralnick/millerjared/BoCP/data/processed/name_list.rds")
#name_list <- readRDS("/blue/guralnick/millerjared/BoCP/plant-db-pipeline/unfinished_names.rds")
#name_list <- readRDS("/blue/guralnick/millerjared/BoCP/plant-db-pipeline/zero-byte-namelist.rds") # rerun, for some reason parquet files were 0 bites for these names...
accepted_target <- name_list[[task_id]][1]
accepted_name_filestyle <- gsub(" ", "_", accepted_target)

## find parquet files for aggregate data
parquet_dir <- "/blue/guralnick/millerjared/BoCP/plant-db-pipeline/parquet-files/"
all_parquet_files <- file.path(parquet_dir, "*.parquet")
## SQL: JOIN archive table to parquet using LIKE, annotating with match (can return multiple matches per pattern)
query <- sprintf("SELECT
                    nra.accepted_name,
                    nra.like_pattern,
                    occ.verbatimScientificName AS like_pattern_matched
                 FROM
                    name_relation_archive AS nra
                 LEFT JOIN
                    parquet_scan('%s') AS occ
                 ON
                  UPPER(occ.verbatimScientificName) LIKE nra.like_pattern
                 WHERE
                  nra.accepted_name = $1
                  AND occ.verbatimScientificName IS NOT NULL
                 GROUP BY
                  nra.accepted_name,
                  nra.like_pattern,
                  occ.verbatimScientificName
                 ", all_parquet_files)

result <- dbGetQuery(con, query, params=list(accepted_target))
result <- result %>% dplyr::distinct()

## Write results to a per-name parquet file if not empty
if(nrow(result) > 0){
  matches_parquet_file <- file.path("/blue/guralnick/millerjared/BoCP/plant-db-pipeline/archived-name-matches/",
                                    paste0(accepted_name_filestyle, "_matches.parquet")) 
  arrow::write_parquet(result, matches_parquet_file)
  cat("Wrote", nrow(result), "matches to", matches_parquet_file, "\n")
} else {
  cat("No matches found for", accepted_target, "\n")
}

cat("Finished and Closing DB connection...\n")
dbDisconnect(con, shutdown = TRUE)
