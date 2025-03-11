# 003C Finalize North American Taxon List, Prune Stephen's Tree 
# Author: JT Miller 
# Date: 03-11-2024
# Project: BoCP 

# load packages
library(data.table)
library(tidyverse)
library(sf)
library(duckdb)
library(arrow)
library(DBI)
library(rgnparser)
library(ape)

# Load North American names 
wcvp_na_name_alignment <- fread("/blue/guralnick/millerjared/BoCP/data/processed/wcvp-ncbi-alignment-na.csv")
names_to_check_ncbi <- fread("/blue/guralnick/millerjared/BoCP/data/processed/wcvp-ncbi-alignment-ncbi-needs-na-check.csv")
names_checked_ncbi <- fread("/blue/guralnick/millerjared/BoCP/outputs/ncbi_na_check/na_name_occ_check_bigmem.csv")

# filter for our checked ncbi names to only include those that fall within North America
ncbi_names_na <- names_checked_ncbi %>% filter(inNA == TRUE)
ncbi_na_alignment <- names_to_check_ncbi %>% filter(alignedParentName %in% ncbi_names_na$acceptedParentName) %>% select(-duplicateSourceForAlignedParentName)

# rebind the name alignments for wcvp and ncbi 
na_name_alignment <- rbind(wcvp_na_name_alignment, ncbi_na_alignment)
fwrite(na_name_alignment, "/blue/guralnick/millerjared/BoCP/data/processed/na-taxonomy-harmonized.csv")

# build an rds object of the list of names for retrieval 
# generate list of names, with the accepted name parent as the first name, with synonyms as the alternatives
accepted_name_v <- unique(na_name_alignment$alignedParentName)
name_list <- list() # initialize an empty list
# reorganize into a nested list
bench::bench_time({
  for(i in 1:length(accepted_name_v)){
    p <- na_name_alignment[alignedParentName == accepted_name_v[[i]]] # create a for-loop that goes through by the accepted name and grabs synonyms, storing them both as a vector in a list. 
    s <- p[!is.na(name)]
    x <- print(s$name)
    a <- print(unique(p$alignedParentName))
    name_list[[i]] <- unique(c(a,x))
  }
})

write_rds(name_list, "/blue/guralnick/millerjared/BoCP/data/processed/name_list.rds")

## Prune stephens tree according to our current knowledge.

# load stephens trees
current_tree_molc_only <- ape::read.tree("/blue/guralnick/millerjared/BoCP/data/processed/tree-outputs/smith-trees-feb-2025/polypod_acro_angio_dated.tre.fam_ord.wcvpmatched")
current_tree_taxa_inf <- ape::read.tree("/blue/guralnick/millerjared/BoCP/data/processed/tree-outputs/smith-trees-feb-2025/polypod_acro_angio_dated.wcvp.tre.fam_ord.wcvpmatched")

# create a vector of 
na_species <- unique(na_name_alignment$alignedParentName)
current_tree_molc_only$tip_labels_trimmed <- gsub("^[^_]+_[^_]+_", "", current_tree_molc_only$tip.label)
current_tree_taxa_inf$tip_labels_trimmed <- gsub("^[^_]+_[^_]+_", "", current_tree_taxa_inf$tip.label)
current_tree_molc_only$tip_labels_trimmed <- gsub("_", " ", current_tree_molc_only$tip_labels_trimmed)
current_tree_taxa_inf$tip_labels_trimmed <- gsub("_", " ", current_tree_taxa_inf$tip_labels_trimmed)
# Find indices to keep
match_molc_indices <- match(na_species, current_tree_molc_only$tip_labels_trimmed)
match_inf_indices <- match(na_species, current_tree_taxa_inf$tip_labels_trimmed)
# Remove NA values from match results
match_molc_indices <- match_molc_indices[!is.na(match_molc_indices)]
match_inf_indices <- match_inf_indices[!is.na(match_inf_indices)]
# Drop tips that are not in na_species
pruned_molc_tree <- drop.tip(current_tree_molc_only, current_tree_molc_only$tip.label[-match_molc_indices])
pruned_inf_tree <- drop.tip(current_tree_taxa_inf, current_tree_taxa_inf$tip.label[-match_inf_indices])
# write these files as outputs to check
write.tree(pruned_molc_tree, "/blue/guralnick/millerjared/BoCP/outputs/pruned-molc-tree.tre")
write.tree(pruned_inf_tree, "/blue/guralnick/millerjared/BoCP/outputs/pruned-tax-inf-tree.tre")


