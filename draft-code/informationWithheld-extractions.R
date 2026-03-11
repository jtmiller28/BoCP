library(DBI)
library(duckdb)
library(arrow)
idigbio_file_path <- "/blue/guralnick/millerjared/BoCP/data/raw/parquet-occs/idigio_full_occ.parquet"
gbif_file_path <- "/blue/guralnick/millerjared/BoCP/data/raw/parquet-occs/gbif_full_occ.parquet"
symbiota_file_path <- "/blue/guralnick/millerjared/BoCP/data/raw/parquet-occs/symbiota_occ.parquet"
# Connect to DuckDB
connect <- dbConnect(duckdb::duckdb())

# Query to get unique values and their counts
query <- "
  SELECT informationWithheld, COUNT(*) AS count
  FROM parquet_scan(?)
  GROUP BY informationWithheld
  ORDER BY count DESC
"

# Execute the query
df_counts <- dbGetQuery(connect, query, params = list(symbiota_file_path))

# Print results
print(df_counts)

# Disconnect from DuckDB
dbDisconnect(connect, shutdown = TRUE)
