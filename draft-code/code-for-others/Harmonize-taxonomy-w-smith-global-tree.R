### Harmonize Taxonomy with Smith & Brown Tree
### Author: JT Miller
### Date: 05/22/2025

# load packages
library(data.table) # my preferred way to read/write files and some large tabular operations
library(rgnparser) # parse genus + specificEpithet from author from scientific name 
library(tidyverse) # suite of fxns for general data wrangling and tidying

#### Set up pathing to use gnparser (may not be a required step depending on machine, but may be needed if your system doesnt automatically point to the install)
# my_path <- Sys.getenv("PATH") # grab our path
# Sys.setenv(PATH = paste0(my_path, "/home/millerjared/gnparser"))


# step 1) Load tables
# Load you taxonomic relational table (for our purposes this will be the table I build that harmonizes all plant names in WCVP and attaches relevant NCBI info)
taxa_table <- fread("./data/processed/wcvp-ncbi-alignment.csv")

# take a look at the data
head(taxa_table) # I structured this as a relational table, you should be able to take any plant name and link it up with this table to harmonize with the most current name.

# load your study species table of names (I dont have one here, so Im going to throw the original that was given to me at the start of BoCP)
bocp_name_table <- fread("./data/raw/USACANMEX_plantNamesFieldsAdded.csv")

# take a look at these data 
head(bocp_name_table) # note that there are three author fields due to comma seperated values, this could be true for any dataset so heres some simple logic to collapse these to one field 
bocp_name_table <- bocp_name_table %>% 
  rename(ExteriorNameAuthorship3 = V5) %>%
  mutate(fullSN = case_when(
    ExteriorNameAuthorship == "" & ExteriorNameAuthorship2 == "" & ExteriorNameAuthorship3 == "" ~ Name, # if just Name
    ExteriorNameAuthorship != "" & ExteriorNameAuthorship2 == "" & ExteriorNameAuthorship3 == "" ~  # if name + first ExternalAuthor
      paste(Name, ExteriorNameAuthorship, sep = ", "),
    ExteriorNameAuthorship != "" & ExteriorNameAuthorship2 != "" & ExteriorNameAuthorship3 == "" ~ # if name + first + second external author
      paste(Name, ExteriorNameAuthorship, ExteriorNameAuthorship2, sep = ", "),
    ExteriorNameAuthorship != "" & ExteriorNameAuthorship2 != "" & ExteriorNameAuthorship3 != "" ~  # if name + first + second + third external authro
      paste(Name, ExteriorNameAuthorship, ExteriorNameAuthorship2, ExteriorNameAuthorship3, sep = ", ")
  )) %>% 
  select(fullSN) # we only care about having this field now...

# now we'll parse the name into the genus + specificEpithet , then the authors
# parse the scientificName field into the name and authorship fields
bocp_parsed_names <- gn_parse_tidy(bocp_name_table$fullSN)  # break up names
# remove duplicates where present 
bocp_parsed_names <- bocp_parsed_names %>% 
  distinct(, .keep_all=TRUE) # sometimes dupes show up on parse step, remove them

# note if any names failed to parse
failed_to_parse <- bocp_parsed_names %>% 
  filter(cardinality == 0) # in this case the '?' in the author name caused issues, in this case the issue arose in the author so it truncated the name there and left the author as NA. This is fine for our purposes (ive removed author matching from our methods by filtering out multiple mapping names)

# now merge names with the harmonized taxonomy 
merged_df <- merge(taxa_table, bocp_parsed_names, by.x = 'name', by.y = 'canonicalfull') # note that 'name' field is the field you should join your verbatim names on, its all the possible names known to wcvp/ncbi taxanomic databases. Extra note!!! Subsp and varities make things...more challenging. I knew this table did not include these, but they could exist I'll illustrate one below.
harmonized_table <- merged_df %>% select(-year, -quality) # I usually toss these 
# now you can check what names failed to align 
failed_names <- setdiff(unique(bocp_parsed_names$canonicalfull), unique(merged_df$name))

# Why did these fail? Most that arent hybrids (cannot be aligned) are names that have multiple names possible with authorship
# Stephen and I elected to remove these names from the dataset considering only 5% of the total names fall under this category
# Theoretically you could use these names as long as you match authors as well, but this adds complications (authors are not as standardized making this a large chore to ensure congruence)
# Additionally, in terms of spatial analyses you will likely have many challenges since most occurrence data does NOT contain authorship info, meaning you would need to filter much of the data by exact matching authors...
# I have code if you absolutely need to do this for other projects, (intially this is how I build the pipelines for BoCP) however I cannot recommend that method unless it will remove focal taxa. 

# Also a note on subsp varieties and forms, the taxa table has these as optional names that you can harmonize onto
hybrid_names <- data.frame(fullSN = c("Agastache pallidiflora var. neomexicana","Agastache pallidiflora var neomexicana", "Agastache pallidiflora variety neomexicana", "Agastache pallidiflora subsp. typica", "Agastache pallidiflora subsp typica", "Agastache pallidiflora subspecies typica", "Agastache pallidiflora ssp. typica", "Agastache pallidiflora spp typica", "Agathisanthemum angolense"))
hybrids_parsed <- gn_parse_tidy(hybrid_names$fullSN)  # break up names
head(hybrids_parsed) # note that this is why we use canonicalfull as our merge field, it incorporates subsp. and var. 
# Also notice that alternative forms of notation for subspecies and varities are taken and though there is some issues occassionally with standardization (subspecies is given its own value instead of being corrected to subsp.)
hybrids_parsed <- hybrids_parsed %>% 
  mutate(canonicalfull = gsub("subspecies", "subsp.", canonicalfull)) # just use a simple gsub to fix this case prior to merge
# now merge
merged_df <- merge(taxa_table, hybrids_parsed, by.x = 'name', by.y = 'canonicalfull')
harmonized_table <- merged_df %>% select(-year, -quality)
# Otherwise thats it!!! Its rather simple when using nicely strucutured data and you dont deal with the edge cases. Edge cases can be: multiple mappers, fuzzy matching (names that are spelled incorrectly according to the db)