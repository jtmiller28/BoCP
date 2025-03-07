### Appending Synonyms to the Native Florida Plants list 
# 11-21-2024



fl_native_names <- fread("./data/processed/fl-plants-project/name-summary/fl-plant-summary.csv")

name_alignment <- fread("./data/processed/finalized_name_alignment_wcvp.csv")


fl_native_names_w_synonyms <- merge(fl_native_names, name_alignment, by.x = "parentTaxon", by.y = "acceptedNameParent", all.x = TRUE)
