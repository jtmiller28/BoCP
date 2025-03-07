library(data.table)
library(dplyr)

name_alignment <- fread("./data/processed/finalized_name_alignment_wcvp.csv")
tax_clean <- fread("./data/processed/cleaning-summaries/taxa_clean_info_tb.csv")
distinct_clean <- fread("./data/processed/cleaning-summaries/distinct_clean_info_tb.csv")
soft_sp_clean <- fread("./data/processed/cleaning-summaries/soft_spatial_clean_info_tb.csv")
full_sp_clean <- fread("./data/processed/cleaning-summaries/full-spatial-clean-info-tb.csv")
tax_clean <- distinct(tax_clean, acceptedNameParent, .keep_all = TRUE)
distinct_clean <- distinct(distinct_clean, name, .keep_all = TRUE)
summary_df <- name_alignment %>% 
  left_join(tax_clean, by = "acceptedNameParent")

distinct_clean <- rename(distinct_clean, acceptedNameParent = name)
soft_sp_clean <- soft_sp_clean %>% 
  rename(acceptedNameParent = name) %>% 
  distinct(acceptedNameParent, .keep_all = TRUE)

full_sp_clean <- full_sp_clean %>% 
  rename(acceptedNameParent = name) %>% 
  distinct(acceptedNameParent, .keep_all = TRUE)

summary_df <- summary_df %>% 
  left_join(distinct_clean, by = "acceptedNameParent")

summary_df <- summary_df %>% 
  left_join(soft_sp_clean, by = "acceptedNameParent")

summary_df <- summary_df %>% 
  left_join(full_sp_clean, by = "acceptedNameParent")

regioned_data <- fread("./data/processed/regioned_summaries/region-summary-tb.csv")

regioned_data <- regioned_data %>% 
  rename(acceptedNameParent = species) %>% 
  distinct(acceptedNameParent, .keep_all = TRUE)

summary_df <- summary_df %>% 
  left_join(regioned_data, by = "acceptedNameParent")

summary_df_distinct <- summary_df %>% 
  distinct(acceptedNameParent, .keep_all = TRUE)

fwrite(summary_df_distinct, "./outputs/example-cleaning-summary.csv")
