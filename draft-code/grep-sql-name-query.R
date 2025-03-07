library(duckdb)
# Load necessary libraries
library(DBI)
library(duckdb)
library(tidyverse)

# Connect to DuckDB
connect <- dbConnect(duckdb::duckdb())

# Example data for names_to_retrieve and accepted_name_v
# names_to_retrieve <- list(c("species1", "species2", "species3"))
# accepted_name_v <- c("accepted_name1", "accepted_name2")

# Retrieve the names to filter for
names_to_retrieve <- c("Acer Rubrum", "Acer Rubanium")
names_to_retrieve <- toupper(names_to_retrieve)
names_to_retrieve <- gsub(" ", "%", names_to_retrieve)

# Define the file path
file_path <- "./data/raw/parquet-occs/idigbio_occ.parquet"

# Construct the SQL query with LIKE for partial matching
query <- paste0("SELECT verbatimScientificName, decimaLlatitude, decimalLongitude FROM '", file_path, "' WHERE UPPER(verbatimScientificName) LIKE ", 
                paste0("'%", names_to_retrieve, "%'", collapse = " OR UPPER(verbatimScientificName) LIKE "))

# Execute the query
df <- dbGetQuery(connect, query)


# Print the result
print(df)

# Close the connection
dbDisconnect(connect)

query <- "SELECT COUNT(*) FROM '/home/millerjared/blue_guralnick/millerjared/BoCP/data/raw/parquet-occs/idigbio_occ.parquet'"
df_count <- dbGetQuery(connect, query)

# Print the result
print(df_count)

# checking some extraneous info on occurrenceID
query <- "SELECT occurrenceID, COUNT(*) AS occurrence_count 
          FROM '/home/millerjared/blue_guralnick/millerjared/BoCP/data/raw/parquet-occs/idigbio_raw_occ.parquet' 
          GROUP BY occurrenceID 
          ORDER BY occurrence_count DESC"

df_occurrence_counts <- dbGetQuery(connect, query)

# Print the first few rows to inspect
head(df_occurrence_counts)

# Testing a particular occID duplicate 
# Construct the SQL query with LIKE for partial matching

occurrence_id_to_check <- "urn:uuid:1f9fd506-3e7e-41b0-a868-c89e1d0195df"  # Replace with the actual ID

query <- paste0("SELECT * 
                 FROM '/home/millerjared/blue_guralnick/millerjared/BoCP/data/raw/parquet-occs/idigbio_raw_occ.parquet' 
                 WHERE occurrenceID = '", occurrence_id_to_check, "'")

df_raw_occurrence <- dbGetQuery(connect, query)

query <- paste0("SELECT * 
                 FROM '/home/millerjared/blue_guralnick/millerjared/BoCP/data/raw/parquet-occs/idigbio_processed_occ.parquet' 
                 WHERE occurrenceID = '", occurrence_id_to_check, "'")

df_processed_occurrence <- dbGetQuery(connect, query)
