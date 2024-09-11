# gatoRs functions: Authors: Natalie Patten & Shelly Gaynor. Original documentation for these can be found here: https://github.com/nataliepatten/gatoRs
## edited by JT Miller for BoCP purposes 2024
## Set-up dependent functions for Gators Download, these remain unchanged ############################################################################################
#' @importFrom dplyr bind_rows distinct rename select case_when
#' @importFrom rgbif occ_data name_backbone
#' @importFrom magrittr "%>%"
#' @importFrom utils write.csv
get_gbif <- function(synonyms.list, gbif.match = "fuzzy", gbif.prov = FALSE, limit = 100000){
  if (gbif.match != "fuzzy" & gbif.match != "code") {
    stop("Invalid value for argument: gbif.match. Value for gbif.match must equal 'fuzzy' or 'code'.")
  }
  if (length(synonyms.list) == 0 | any(is.na(synonyms.list))) {
    stop("Invalid argument: synonyms.list. The argument synonyms.list must be non-empty.")
  }
  
  colNames <- c("scientificName",
                "genus",
                "specificEpithet",
                "infraspecificEpithet",
                "key",
                "occurrenceID",
                "basisOfRecord",
                "eventDate",
                "year",
                "month",
                "day",
                "institutionCode",
                "recordedBy",
                "country",
                "county",
                "stateProvince",
                "locality",
                "occurrenceRemarks",
                "verbatimLocality",
                "decimalLatitude",
                "verbatimLatitude",
                "decimalLongitude",
                "verbatimLongitude",
                "coordinateUncertaintyInMeters",
                "informationWithheld",
                "habitat",
                "geodeticDatum")
  numerical <- c("decimalLatitude",
                 "verbatimLatitude",
                 "decimalLongitude",
                 "verbatimLongitude",
                 "year",
                 "month",
                 "day",
                 "coordinateUncertaintyInMeters")
  
  query_gbif <- data.frame(matrix(ncol = length(colNames), nrow = 0))
  colnames(query_gbif) <- colNames
  # fix columns to be of type character
  for (i in 1:NCOL(query_gbif)) {
    if (colNames[i] %in% numerical) {
      query_gbif[,i] <- as.numeric(query_gbif[,i])
    }
    else {
      query_gbif[,i] <- as.character(query_gbif[,i])
    }
  }
  
  if (gbif.match == "code") {
    for (i in 1:length(synonyms.list)) {
      key <- suppressWarnings(rgbif::name_backbone(name = synonyms.list[i])$speciesKey)
      if (!is.null(key)) {
        temp <- rgbif::occ_data(taxonKey = key, limit = limit, curlopts=list(http_version=2))
        temp <- temp$data
        query_gbif <- dplyr::bind_rows(query_gbif, temp)
      }
    }
  }
  else if (gbif.match == "fuzzy") {
    for (i in 1:length(synonyms.list)) {
      temp <- rgbif::occ_data(scientificName = synonyms.list[i], limit = limit, curlopts=list(http_version=2))
      temp <- temp$data
      # use bind_rows() to account for different number of columns
      query_gbif <- dplyr::bind_rows(query_gbif, temp)
    }
  }
  
  # if no results found
  if (NROW(query_gbif) == 0) {
    message("No results found in GBIF.")
    return(query_gbif)
  }
  
  if (gbif.prov) {
    query_gbif <- dplyr::distinct(query_gbif, key, .keep_all = TRUE)
    query_gbif$key <- as.numeric(query_gbif$key)
    query_gbif <- rgbif::occ_get_verbatim(key = query_gbif$key, fields = colNames, curlopts=list(http_version=2))
  }
  
  temp <- data.frame(matrix(NA, ncol = 0, nrow = 0))
  tempColNames <- colnames(temp)
  
  for (i in 1:length(colNames)) {
    if (!(colNames[i] %in% colnames(query_gbif))) {
      temp <- data.frame(matrix(NA, ncol = NCOL(temp) + 1, nrow = 0))
      colnames(temp) <- c(tempColNames, colNames[i])
      tempColNames <- colnames(temp)
    }
  }
  
  if (NCOL(temp) > 0) {
    for (i in 1:NCOL(temp)) {
      if (colNames[i] %in% numerical) {
        temp[,i] <- as.numeric(temp[,i])
      }
      else {
        temp[,i] <- as.character(temp[,i])
      }
    }
  }
  
  query_gbif <- dplyr::bind_rows(temp, query_gbif)
  query_gbif <- dplyr::rename(query_gbif, ID = "key",
                              latitude = "decimalLatitude",
                              verbatimLatitude = "verbatimLatitude",
                              longitude = "decimalLongitude",
                              verbatimLongitude = "verbatimLongitude")
  query_gbif <- suppressWarnings(correct_class(query_gbif))
  
  # Add occurrenceRemarks, verbatimLocality to locality column
  query_gbif$locality <- paste("locality: ", query_gbif$locality)
  query_gbif$locality <- paste(query_gbif$locality, query_gbif$occurrenceRemarks, sep = ", occurrenceRemarks: ")
  query_gbif$locality <- paste(query_gbif$locality, query_gbif$verbatimLocality, sep = ", verbatimLocality: ")
  
  # if decimal lat/lon columns are empty, replace with verbatim lat/lon columns
  temp <- query_gbif[is.na(query_gbif$latitude), ]
  if (NROW(temp) > 0) {
    query_gbif <- query_gbif[!(is.na(query_gbif$latitude)), ]
    for (i in 1:NROW(temp)) {
      if (!is.na(temp$verbatimLatitude[i]))
        temp$latitude[i] <- temp$verbatimLatitude[i]
    }
    query_gbif <- rbind(query_gbif, temp)
  }
  
  temp <- query_gbif[is.na(query_gbif$longitude), ]
  if (NROW(temp) > 0) {
    query_gbif <- query_gbif[!(is.na(query_gbif$longitude)), ]
    for (i in 1:NROW(temp)) {
      if (!is.na(temp$verbatimLongitude[i]))
        temp$longitude[i] <- temp$verbatimLongitude[i]
    }
    query_gbif <- rbind(query_gbif, temp)
  }
  
  query_gbif$aggregator <- "GBIF"
  query_gbif <- query_gbif %>%
    dplyr::select("scientificName",
                  "genus",
                  "specificEpithet",
                  "infraspecificEpithet",
                  "ID",
                  "occurrenceID",
                  "basisOfRecord",
                  "eventDate",
                  "year",
                  "month",
                  "day",
                  "institutionCode",
                  "recordedBy",
                  "country",
                  "county",
                  "stateProvince",
                  "locality",
                  "latitude",
                  "longitude",
                  "coordinateUncertaintyInMeters",
                  "informationWithheld",
                  "habitat",
                  "aggregator") %>%
    suppressWarnings(correct_class())
  
  return(query_gbif)
}

correct_class <- function(df, scientific.name = "scientificName", genus = "genus",
                          species = "specificEpithet", infraspecific.epithet = "infraspecificEpithet",
                          id = "ID", occ.id = "occurrenceID",
                          basis.of.record = "basisOfRecord", event.date = "eventDate",
                          year = "year", month = "month", day = "day",
                          inst.code = "institutionCode", recorded.by = "recordedBy",
                          country = "country", county = "county", state = "stateProvince",
                          locality = "locality", latitude = "latitude",
                          longitude = "longitude",
                          coord.uncertainty = "coordinateUncertaintyInMeters",
                          info.withheld = "informationWithheld", habitat = "habitat",
                          aggregator = "aggregator"){
  
  df[[scientific.name]] <- dplyr::case_when(df[[scientific.name]] == "" ~ NA, .default = as.character(df[[scientific.name]]))
  df[[genus]] <- dplyr::case_when(df[[genus]] == "" ~ NA, .default = as.character(df[[genus]]))
  df[[species]] <- dplyr::case_when(df[[species]] == "" ~ NA, .default = as.character(df[[species]]))
  df[[basis.of.record]] <- dplyr::case_when(df[[basis.of.record]] == "" ~ NA, .default = as.character(df[[basis.of.record]]))
  df[[event.date]] <- dplyr::case_when(df[[event.date]] == "" ~ NA, .default = as.character(df[[event.date]]))
  df[[year]] <- dplyr::case_when(df[[year]] == "" ~ NA, .default = as.character(df[[year]]))
  df[[month]] <- dplyr::case_when(df[[month]] == "" ~ NA, .default = as.character(df[[month]]))
  df[[day]] <- dplyr::case_when(df[[day]] == "" ~ NA, .default = as.character(df[[day]]))
  df[[inst.code]] <- dplyr::case_when(df[[inst.code]] == "" ~ NA, .default = as.character(df[[inst.code]]))
  df[[country]] <- dplyr::case_when(df[[country]] == "" ~ NA, .default = as.character(df[[country]]))
  df[[county]] <- dplyr::case_when(df[[county]] == "" ~ NA, .default = as.character(df[[county]]))
  df[[state]] <- dplyr::case_when(df[[state]] == "" ~ NA, .default = as.character(df[[state]]))
  df[[locality]] <- dplyr::case_when(df[[locality]] == "" ~ NA, .default = as.character(df[[locality]]))
  df[[latitude]] <- dplyr::case_when(df[[latitude]] == "" ~ NA, .default = as.numeric(df[[latitude]]))
  df[[longitude]]<- dplyr::case_when(df[[longitude]] == "" ~ NA, .default = as.numeric(df[[longitude]]))
  df[[id]] <- dplyr::case_when(df[[id]] == "" ~ NA, .default = as.character(df[[id]]))
  df[[coord.uncertainty]] <- dplyr::case_when(df[[coord.uncertainty]] == "" ~ NA, .default = as.character(df[[coord.uncertainty]]))
  df[[info.withheld]] <- dplyr::case_when(df[[info.withheld]] == "" ~ NA, .default = as.character(df[[info.withheld]]))
  df[[habitat]] <- dplyr::case_when(df[[habitat]] == "" ~ NA, .default = as.character(df[[habitat]]))
  df[[occ.id]] <- dplyr::case_when(df[[occ.id]] == "" ~ NA, .default = as.character(df[[occ.id]]))
  
  if (infraspecific.epithet %in% colnames(df)){
    df[[infraspecific.epithet]] <- dplyr::case_when(df[[infraspecific.epithet]] == "" ~ NA, .default = as.character(df[[infraspecific.epithet]]))
  } else {
    df[[infraspecific.epithet]] <- NA
  }
  
  if (aggregator %in% colnames(df)) {
    df[[aggregator]] <- dplyr::case_when(df[[aggregator]] == "" ~ NA, .default = as.character(df[[aggregator]]))
  }
  
  if (recorded.by %in% colnames(df)) {
    df[[recorded.by]] <- dplyr::case_when(df[[recorded.by]] == "" ~ NA, .default = as.character(df[[recorded.by]]))
  }
  
  return(df)
}
######################################################################################################################################################################
## Gators Download: For our purposes we're only going to query one API idigbio, as gbifdb exists. My edits will just be turning off the gbif api search to reduce time 
gators_download_edited <- function(synonyms.list, write.file = FALSE, filename = NA,
                                   gbif.match = "fuzzy", gbif.prov = FALSE,
                                   idigbio.filter = TRUE, limit = 100000, call.idigbio = TRUE, call.gbif = TRUE) { # adding call idigbio and gbif boolean arguments
  
  # Check for valid arguments
  if (length(synonyms.list) == 0 | any(is.na(synonyms.list))) {
    stop("Invalid argument: synonyms.list. The argument synonyms.list must be non-empty.")
  }
  
  if (gbif.match != "fuzzy" & gbif.match != "code") {
    stop("Invalid value for argument: gbif.match. Value for gbif.match must equal 'fuzzy' or 'code'.")
  }
  
  if (idigbio.filter != TRUE & idigbio.filter != FALSE) {
    stop("Invalid value for argument: idigbio.filter. Value for idigbio.filter must equal 'TRUE' or 'FALSE'.")
  }
  
  if (write.file != TRUE & write.file != FALSE) {
    stop("Invalid value for argument: write.file. Value for write.file must equal 'TRUE' or 'FALSE'.")
  }
  else if (write.file) {
    if (is.na(filename)) {
      stop("Invalid value for argument: filename. The location and name of the output file is not specified.")
    }
    
    if (grepl(".csv", filename) == FALSE) {
      stop("Invalid value for argument: filename. The output file name must end in '.csv'.")
    }
  }
  else if (! is.na(filename)) {
    message("Warning: No output file will be written; the filename argument will be ignored.\nTo write to an output file, set write.file = TRUE.")
  }
  
  output <- data.frame() # Initialize outpus as an empty df (possibly get rid of this)
  
  if(call.idigbio){ # make calling idigbio a conditional on argument value 
    query_idigbio <- fix_names(get_idigbio(synonyms.list, limit = limit)) # initial download, fix capitalization
    if (NROW(query_idigbio) > 0) query_idigbio <- dplyr::distinct(query_idigbio, ID, .keep_all = TRUE) # Remove duplicates - records that share UUIDs or KEYs
    query_idigbio <- fix_names(fix_columns(query_idigbio)) # fill out remaining taxon columns, and fix capitalization again
    
    
    if (idigbio.filter) {
      query_idigbio <- filter_fix_names(query_idigbio, synonyms.list)
    }
    else {
      message("Warning: iDigBio search will return all records where any column has a matching string to the provided scientific names.")
    }
    output <- query_idigbio # maybe get rid of this 
  }
  
  if(call.gbif){
    query_gbif <- fix_names(get_gbif(synonyms.list, gbif.match = gbif.match, gbif.prov = gbif.prov, limit = limit)) # initial download, fix capitalization
    if (NROW(query_gbif) > 0) query_gbif <- dplyr::distinct(query_gbif, ID, .keep_all = TRUE) # Remove duplicates - records that share UUIDs or KEYs
    query_gbif <- fix_names(fix_columns(query_gbif)) # fill out remaining taxon columns, and fix capitalization again
    
    output <- if(nrow(output) == 0){ # set up a conditional to operate based on whether call.idigbio happened
      query_gbif
    } else {
      dplyr::bind_rows(output, query_gbif)
    }
  }
  
  if(nrow(output) == 0){
    stop("No Records Found.")
  }
  
  if (write.file) {
    utils::write.csv(output, filename, row.names = FALSE)
  }
  return(output)
}

### Correct Class: https://github.com/nataliepatten/gatoRs/blob/main/R/correct_class.R
correct_class <- function(df, scientific.name = "scientificName", scientific.name.authorship = "scientificNameAuthorship",
                          id = "ID", occ.id = "occurrenceID",
                          basis.of.record = "basisOfRecord", event.date = "eventDate",
                          year = "year", month = "month", day = "day",
                          inst.code = "institutionCode", recorded.by = "recordedBy",
                          country = "country", county = "county", state = "stateProvince",
                          locality = "locality", latitude = "latitude",
                          longitude = "longitude",
                          coord.uncertainty = "coordinateUncertaintyInMeters",
                          info.withheld = "informationWithheld", habitat = "habitat",
                          aggregator = "aggregator"){
  
  df[[scientific.name]] <- dplyr::case_when(df[[scientific.name]] == "" ~ NA, .default = as.character(df[[scientific.name]]))
  df[[scientific.name.authorship]] <- dplyr::case_when(df[[scientific.name.authorship]] == "" ~ NA, .default = as.character(df[[scientific.name.authorship]]))
  df[[basis.of.record]] <- dplyr::case_when(df[[basis.of.record]] == "" ~ NA, .default = as.character(df[[basis.of.record]]))

  df[[event.date]] <- dplyr::case_when(
    trimws(as.character(df[[event.date]])) == "" ~ NA_character_,
    TRUE ~ as.character(as.Date(df[[event.date]], format = "%Y-%m-%d"))
  )
  df[[year]] <- dplyr::case_when(df[[year]] == "" ~ NA, .default = as.character(df[[year]]))
  df[[month]] <- dplyr::case_when(df[[month]] == "" ~ NA, .default = as.character(df[[month]]))
  df[[day]] <- dplyr::case_when(df[[day]] == "" ~ NA, .default = as.character(df[[day]]))
  df[[inst.code]] <- dplyr::case_when(df[[inst.code]] == "" ~ NA, .default = as.character(df[[inst.code]]))
  df[[country]] <- dplyr::case_when(df[[country]] == "" ~ NA, .default = as.character(df[[country]]))
  df[[county]] <- dplyr::case_when(df[[county]] == "" ~ NA, .default = as.character(df[[county]]))
  df[[state]] <- dplyr::case_when(df[[state]] == "" ~ NA, .default = as.character(df[[state]]))
  df[[locality]] <- dplyr::case_when(df[[locality]] == "" ~ NA, .default = as.character(df[[locality]]))
  df[[latitude]] <- dplyr::case_when(df[[latitude]] == "" ~ NA, .default = as.numeric(df[[latitude]]))
  df[[longitude]]<- dplyr::case_when(df[[longitude]] == "" ~ NA, .default = as.numeric(df[[longitude]]))
  df[[id]] <- dplyr::case_when(df[[id]] == "" ~ NA, .default = as.character(df[[id]]))
  df[[coord.uncertainty]] <- dplyr::case_when(df[[coord.uncertainty]] == "" ~ NA, .default = as.character(df[[coord.uncertainty]]))
  df[[info.withheld]] <- dplyr::case_when(df[[info.withheld]] == "" ~ NA, .default = as.character(df[[info.withheld]]))
  df[[habitat]] <- dplyr::case_when(df[[habitat]] == "" ~ NA, .default = as.character(df[[habitat]]))
  df[[occ.id]] <- dplyr::case_when(df[[occ.id]] == "" ~ NA, .default = as.character(df[[occ.id]]))
  if (aggregator %in% colnames(df)) {
    df[[aggregator]] <- dplyr::case_when(df[[aggregator]] == "" ~ NA, .default = as.character(df[[aggregator]]))
  }
  
  if (recorded.by %in% colnames(df)) {
    df[[recorded.by]] <- dplyr::case_when(df[[recorded.by]] == "" ~ NA, .default = as.character(df[[recorded.by]]))
  }
  # Convert eventDate to Date type
  # df[[event.date]] <- as.Date(df[[event.date]], format = "%Y-%m-%d")
  return(df)
}

### Modifying remove duplicates to include collector. 
#' @return Return data frame with duplicates removed.
#' Information about the columns in the returned data frame can be found in the documentation for `gators_download()`.
#'
#' @importFrom parsedate parse_iso_8601 format_iso_8601
#' @importFrom dplyr distinct mutate select row_number filter
#' @importFrom stats na.omit
#'
#' @export

remove_duplicates_mod <- function(df, event.date = "eventDate",
                              aggregator = "aggregator", id = "ID", occ.id = "occurrenceID",
                              year = "year", month = "month", day = "day",
                              roundedLatitude = "roundedLatitude", roundedLongitude = "roundedLongitude", recordedBy = "recordedBy",
                              remove.NA.occ.id = FALSE, remove.NA.date = FALSE,
                              remove.unparseable = FALSE){
  
  if (NROW(df) == 0) return(df)
  
  if (remove.NA.occ.id == TRUE) {
    df <- df[!is.na(df[[occ.id]]), ]
  }
  if (remove.NA.date == TRUE) {
    df <- df[!is.na(df[[event.date]]), ]
  }
  
  # Remove within aggregation duplicates based on ID (UUID or KEY)
  ## This is done within our download function as well
  ag <- unique(df[[aggregator]])
  tempdf <- c()
  for(i in 1:length(ag)){
    tempdf[[i]] <- df[df$aggregator == ag[i], ] # edited this as the functionality is currently weird
    if(nrow(tempdf[[i]]) == 0) next # SKIP (added by JT possibly remove)
    if( (length(unique(stats::na.omit(tempdf[[i]][[id]])))) != nrow(stats::na.omit(tempdf[[i]][id])) ) {
      tempdf[[i]] <- dplyr::distinct(tempdf[[i]], .data[[id]], .keep_all = TRUE)
    }
  }
  df <- do.call(rbind, tempdf)
  # Remove specimen duplicates
  
  has_date_cols <- function(df, i) {
    # individual year, month, day columns don't even exist
    if (length(df[[year]]) == 0 & length(df[[month]]) == 0 & length(df[[day]]) == 0) {
      return(FALSE)
    } else if (is.na(df[[year]][i]) | is.na(df[[month]][i])) {
      return(FALSE)
    } else {
      return(TRUE)
    }
  }
  
  get_temp_date <- function(date, remove.unparseable) {
    
    tryCatch (
      return(as.character(parsedate::parse_iso_8601(parsedate::format_iso_8601(date)))),
      error=function(e) {
        message("Event date cannot be automatically parsed for date: ", date)
        if (remove.unparseable) {
          return("remove")
        } else {
          temp_message <- paste0("Date cannot be automatically parsed for the eventDate: ", date)
          temp_year <- readline(prompt = paste(temp_message, "Please enter the year in YYYY format, or NA if not provided: ", sep = "\n"))
          temp_month <- readline(prompt = "Please enter the month in MM format, or NA if not provided: ")
          temp_day <- readline(prompt = "Please enter the day in DD format or NA if not provided: ")
          temp_date <- paste(temp_year, temp_month, temp_day, sep="-")
          return(temp_date)
        }
      }
    )
  }
  
  new_df <- df
  to_remove <- c()
  for(i in 1:nrow(df)) {
    # If year, month, day are not available, but eventDate is, attempt to parse
    if(!is.na(df[[event.date]][i]) & !has_date_cols(df, i)) {
      temp_date <- get_temp_date(df[[event.date]][i], remove.unparseable)
      
      if (temp_date == "remove" || is.na(temp_date)) { # im adding the is.na(temp_date) check as weird invalid date forms like "0-09-05" will cause NA rather than be caught by the tryCatch error during parsing. 
        to_remove <- append(to_remove, i)
      } else {
        if (substr(temp_date, 1, 4) == "NA") {
          new_df[[year]][i] <- NA
        } else{
          new_df[[year]][i] <- substr(temp_date, 1, 4)
        }
        if (substr(temp_date, 6, 7) == "NA") {
          new_df[[month]][i] <- NA
        } else{
          new_df[[month]][i] <- substr(temp_date, 6, 7)
        }
        if (substr(temp_date, 9, 10) == "NA") {
          new_df[[day]][i] <- NA
        } else{
          new_df[[day]][i] <- substr(temp_date, 9, 10)
        }
      }
    } else { 
      if (length(df[[year]]) == 0) {new_df[[year]][i] <- NA} else {new_df[[year]][i] <- df[[year]][i]}
      if (length(df[[year]]) == "") {new_df[[year]][i] <- NA} else {new_df[[year]][i] <- df[[year]][i]}
      if (length(df[[month]]) == 0) {new_df[[month]][i] <- NA} else {new_df[[month]][i] <- df[[month]][i]}
      if (length(df[[month]]) == "") {new_df[[month]][i] <- NA} else {new_df[[month]][i] <- df[[month]][i]}
      if (length(df[[day]]) == 0) {new_df[[day]][i] <- NA} else {new_df[[day]][i] <- df[[day]][i]}
      if (length(df[[day]]) == "") {new_df[[day]][i] <- NA} else {new_df[[day]][i] <- df[[day]][i]}
    }
  }
  
  
  if (remove.unparseable == TRUE) {new_df <- dplyr::filter(new_df, !(row_number() %in% to_remove))}
  
  new_df[[year]] <- as.numeric(new_df[[year]])
  new_df[[month]] <- as.numeric(new_df[[month]])
  new_df[[day]] <- as.numeric(new_df[[day]])
  
  # Take the df and organize in a way where we can take prority in distinct for coordUncertainty
  new_df <- new_df %>% 
    group_by(recordedBy, year, month, day, roundedLatitude, roundedLongitude) %>% 
    arrange(coordinateUncertaintyInMeters)
  
  # If occurrence ID column exists, remove rows with identical occurrence ID and date
  if (length(df[[occ.id]]) == 0) {
    # Create a temporary column with unique values so NAs in dates are kept
    # https://stackoverflow.com/questions/66537554/keeping-all-nas-in-dplyr-distinct-function
    new_df <- dplyr::mutate(new_df, temp = dplyr::row_number() * (is.na(year) & is.na(month) & is.na(day)))
    # Adding recordedBy as the collector field here.
    new_df <- dplyr::distinct(new_df, .data[[recordedBy]], .data[[year]], .data[[month]], .data[[day]], .data[[roundedLatitude]], .data[[roundedLongitude]], temp, .keep_all = TRUE)
    new_df <- dplyr::select(new_df, -temp)
  } else if (length(df[[occ.id]]) > 0) {
    new_df <- dplyr::mutate(new_df, temp = row_number() * (is.na(year) & is.na(month) & is.na(day)))
    # Adding recordedBy as the collector field here. 
    new_df <- dplyr::distinct(new_df, .data[[recordedBy]], .data[[year]], .data[[month]],
                              .data[[day]], .data[[occ.id]],  .data[[roundedLatitude]], .data[[roundedLongitude]], temp, .keep_all = TRUE)
    new_df <- dplyr::select(new_df, -temp)
  }
  
  return(new_df)
}
