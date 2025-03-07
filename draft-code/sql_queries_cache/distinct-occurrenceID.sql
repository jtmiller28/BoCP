-- !preview conn=DBI::dbConnect(RSQLite::SQLite())
-- Method for taking labeling occurrenceIDs, then ranking ones to filter out. Standard 'query' format
WITH occurrence_counts AS (
    -- Step 1: Identify duplicate types
    SELECT 
        occurrenceID, 
        COUNT(*) AS occ_count, 
        COUNT(DISTINCT verbatimScientificName) AS unique_names
    FROM '/home/millerjared/blue_guralnick/millerjared/BoCP/data/raw/parquet-occs/idigbio_raw_occ.parquet'
    GROUP BY occurrenceID
),

classified_duplicates AS (
    -- Step 2: Assign dupID classification
    SELECT 
        o.*,  
        CASE 
            WHEN oc.occ_count = 1 THEN NULL  -- No duplicate, set dupID to NULL
            WHEN oc.unique_names = 1 THEN 'organismalDup'  -- All verbatimScientificName values match
            ELSE 'misAppliedDup'  -- verbatimScientificName values do not match
        END AS dupID
    FROM '/home/millerjared/blue_guralnick/millerjared/BoCP/data/raw/parquet-occs/idigbio_raw_occ.parquet' o
    JOIN occurrence_counts oc 
    ON o.occurrenceID = oc.occurrenceID
),

ranked_records AS (
    -- Step 3: Rank records based on the number of non-NULL fields
    SELECT *, 
        ( 
            -- Count non-null fields dynamically
            (CASE WHEN verbatimDecimalLatitude IS NOT NULL AND TRIM(verbatimDecimalLatitude) <> '' THEN 1 ELSE 0 END) +
            (CASE WHEN verbatimDecimalLongitude IS NOT NULL AND TRIM(verbatimDecimalLongitude) <> '' THEN 1 ELSE 0 END) +
            (CASE WHEN verbatimYear IS NOT NULL AND TRIM(verbatimYear) <> '' THEN 1 ELSE 0 END) +
            (CASE WHEN verbatimMonth IS NOT NULL AND TRIM(verbatimMonth) <> '' THEN 1 ELSE 0 END) +
            (CASE WHEN verbatimDay IS NOT NULL AND TRIM(verbatimDay) <> '' THEN 1 ELSE 0 END) +
            (CASE WHEN verbatimGeodeticDatum IS NOT NULL AND TRIM(verbatimGeodeticDatum) <> '' THEN 1 ELSE 0 END) +
            (CASE WHEN verbatimCoordinateUncertaintyInMeters IS NOT NULL AND TRIM(verbatimCoordinateUncertaintyInMeters) <> '' THEN 1 ELSE 0 END) +
            (CASE WHEN verbatimInformationWithheld IS NOT NULL AND TRIM(verbatimInformationWithheld) <> '' THEN 1 ELSE 0 END) 
        ) AS non_null_count,
        ROW_NUMBER() OVER (PARTITION BY occurrenceID ORDER BY 
            -- Prefer records with the most non-null fields
            (CASE WHEN verbatimDecimalLatitude IS NOT NULL AND TRIM(verbatimDecimalLatitude) <> '' THEN 1 ELSE 0 END) +
            (CASE WHEN verbatimDecimalLongitude IS NOT NULL AND TRIM(verbatimDecimalLongitude) <> '' THEN 1 ELSE 0 END) +
            (CASE WHEN verbatimYear IS NOT NULL AND TRIM(verbatimYear) <> '' THEN 1 ELSE 0 END) +
            (CASE WHEN verbatimMonth IS NOT NULL AND TRIM(verbatimMonth) <> '' THEN 1 ELSE 0 END) +
            (CASE WHEN verbatimDay IS NOT NULL AND TRIM(verbatimDay) <> '' THEN 1 ELSE 0 END) +
            (CASE WHEN verbatimGeodeticDatum IS NOT NULL AND TRIM(verbatimGeodeticDatum) <> '' THEN 1 ELSE 0 END) +
            (CASE WHEN verbatimCoordinateUncertaintyInMeters IS NOT NULL AND TRIM(verbatimCoordinateUncertaintyInMeters) <> '' THEN 1 ELSE 0 END) +
            (CASE WHEN verbatimInformationWithheld IS NOT NULL AND TRIM(verbatimInformationWithheld) <> '' THEN 1 ELSE 0 END) DESC
        ) AS rank
    FROM classified_duplicates
)
SELECT * FROM ranked_records
WHERE 
    CASE 
        WHEN dupID IS NULL THEN TRUE  -- Keep records without duplicates
        WHEN dupID = 'organismalDup' AND rank = 1 THEN TRUE  -- Keep only the best record for organismal duplicates
        WHEN dupID NOT IN ('misAppliedDup') AND rank = 1 THEN TRUE  -- Ensure ranked filtering applies
        ELSE FALSE 
    END
