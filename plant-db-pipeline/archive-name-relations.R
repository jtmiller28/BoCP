### Archive Known Name Relations in Major DBs (Archive All Possible Matches)
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
db_path <- Sys.getenv("TASK_DB")
con <- dbConnect(duckdb::duckdb(), dbdir = db_path)
# Make output table name unique per array task
target_table <- paste0("name_relation_archive_matches_", array_id)
cat("Writing results to table:", target_table, "\n")

## Prepare a target table for storing matches if not exists
#DBI::dbExecute(con, sprintf("DROP TABLE IF EXISTS %s;", target_table))
DBI::dbExecute(con, sprintf("
  CREATE TABLE IF NOT EXISTS %s (
    accepted_name VARCHAR,
    like_pattern VARCHAR,
    like_pattern_matched VARCHAR
  );", target_table))

## For each accepted_name, match all patterns that are built by the name variants
name_list <- readRDS("/blue/guralnick/millerjared/BoCP/data/processed/name_list.rds")
accepted_target <- name_list[[task_id]][1]

# gbif_file_path <- "/blue/guralnick/millerjared/BoCP/plant-db-pipeline/parquet-files/gbif_full_occ.parquet"
# idigbio_file_path <- "/blue/guralnick/millerjared/BoCP/plant-db-pipeline/parquet-files/idigio_full_occ.parquet"
# symbiota_file_path <- "/blue/guralnick/millerjared/BoCP/plant-db-pipeline/parquet-files/symbiota_occ.parquet"
parquet_dir <- "/blue/guralnick/millerjared/BoCP/plant-db-pipeline/parquet-files/"
all_parquet_files <- file.path(parquet_dir, "*.parquet")
## SQL: JOIN archive table to parquet using LIKE, annotating with match (can return multiple matches per pattern)
## Use GROUP BY to remove redundant name queries (DISTINCT is also an option, albiet slower)
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

## Prepare a target table for storing matches if not exists
# DBI::dbExecute(con, "DROP TABLE IF EXISTS name_relation_archive_matches;")
# dbExecute(con, "
#           CREATE TABLE IF NOT EXISTS name_relation_archive_matches (
#           accepted_name VARCHAR,
#           like_pattern VARCHAR, 
#           like_pattern_matched VARCHAR
#           )")
## insert results 
# if(nrow(result) > 0){
#   dbWriteTable(con, "name_relation_archive_matches", result, append = TRUE)
# }
if(nrow(result) > 0){
  dbWriteTable(con, target_table, result, append = TRUE)
}

cat("Finished and Closing DB connection...\n")
dbDisconnect(con, shutdown = TRUE)
