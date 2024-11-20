### A script for identifying outliers
# Author: JT Miller 
# Date: 09-26-2024
# Project: BoCP 

# Load packages
library(data.table)
library(sf)
library(dplyr)
library(ggplot2)
library(CoordinateCleaner)
library(spThin)
library(gridExtra)

# Set up array logic
start_num <- as.numeric(Sys.getenv("START_NUM"))
task_id <- as.numeric(start_num)
part <- paste0("part", task_id)

# set dir
setwd("/home/millerjared/blue_guralnick/millerjared/BoCP/")


# use name alignment to read in data
na_plants_alignment <- fread("./data/processed/finalized_name_alignment_wcvp.csv")

# load shapefile
level3_regions <- read_sf("./data/raw/level3-wgsrpd/level3.shp")

# Identify some candidates for outlier testing purposes 
regioned_cleaned_summary <- fread("./data/processed/cleaning-summaries/native-range-summary-tb.csv")

# filter if all range counts are zero (and therefore not good for modeling)
regioned_cleaned_summary <- regioned_cleaned_summary %>% 
  filter(!(nativeRangeCount == 0 & !introducedRangeCount == 0))

# plot the distribution of data...
ggplot(regioned_cleaned_summary, aes(x = log(nativeRangeCount + 1))) + 
  geom_histogram(fill = "lightblue", color = "black")

# quantile the data for splits 
regioned_cleaned_summary <- regioned_cleaned_summary %>% 
  mutate(logNativeRangeCount = log(nativeRangeCount + 1))

lower_quantile <- 0.2 # puts us at 7 records, which is probably the minimum we'd want to employ any model...
upper_quantile <- 0.9 # in the 2000s 

quantiles <- quantile(regioned_cleaned_summary$logNativeRangeCount, probs = c(lower_quantile, upper_quantile), na.rm = TRUE)

lower_quantile_label <- paste0(lower_quantile*100, "th Percentile")
upper_quantile_label <- paste0(upper_quantile*100, "th Percentile")

ggplot(regioned_cleaned_summary, aes(x = log(nativeRangeCount + 1))) + 
  geom_histogram(binwidth = function(x) { 2 * IQR(x) / length(x)^(1/3) }, fill = "lightblue", color = "black") + # uses freedman-diaconis Rule to adjust binwidth
  geom_vline(xintercept = quantiles[1], linetype = "dashed", color = "red") + 
  geom_vline(xintercept = quantiles[2], linetype = "dashed", color = "red") + 
  annotate("text", x = quantiles[1], y = Inf, label = lower_quantile_label, vjust = 2, color = "red") +
  annotate("text", x = quantiles[2], y = Inf, label = upper_quantile_label, vjust = 2, color = "red") +
  theme_bw() + 
  ggtitle("Distribution of Occurrence Data Counts with Quantiles") +
  ylab("Frequency") + 
  xlab("Log(Number of Records)")

# subset by quantiles
lower_quantile_subset <- regioned_cleaned_summary %>% 
  filter(logNativeRangeCount <= quantiles[1])
middle_quantile_subset <- regioned_cleaned_summary %>% 
  filter(logNativeRangeCount > quantiles[1] & logNativeRangeCount <= quantiles[2])
upper_quantile <- regioned_cleaned_summary %>% 
  filter(logNativeRangeCount > quantiles[2])

# take random samples from each of these subsets
set.seed(100)
sample_lower <- sample_n(lower_quantile_subset,  25)
sample_middle <- sample_n(middle_quantile_subset, 50)
sample_upper <- sample_n(upper_quantile, 25)

sample_lower$quantile <- "lower"
sample_middle$quantile <- "middle"
sample_upper$quantile <- "upper"

random_sample <- bind_rows(sample_lower, sample_middle, sample_upper)

# build a vector of names to index by, then index according to array tasks
accepted_name_v <- unique(na_plants_alignment$acceptedNameParent)
accepted_name_v <- intersect(accepted_name_v, random_sample$parentTaxon)
accepted_name <- accepted_name_v[task_id] # index by task ID
accepted_name_filestyle <- gsub(" ", "-", accepted_name)

dir.create(paste0("./outputs/testing-outlier-removal/", accepted_name_filestyle))
dir.create(paste0("./outputs/testing-outlier-removal/", accepted_name_filestyle, "/all-data/"))
dir.create(paste0("./outputs/testing-outlier-removal/", accepted_name_filestyle, "/all-data/250km"))
dir.create(paste0("./outputs/testing-outlier-removal/", accepted_name_filestyle, "/native-data"))
dir.create(paste0("./outputs/testing-outlier-removal/", accepted_name_filestyle, "/native-data/250km"))
dir.create(paste0("./outputs/testing-outlier-removal/", accepted_name_filestyle, "/native-plus-introduced-data/250km/"))
dir.create(paste0("./outputs/testing-outlier-removal/", accepted_name_filestyle, "/native-plus-introduced-data"))
dir.create(paste0("./outputs/testing-outlier-removal/", accepted_name_filestyle, "/native-plus-introduced-data/250km"))
### Build region plots for viewing purposes ####################################################################

# load wcvp distribution data
wcvp_distribution_tb <- fread("./data/raw/wcvp_distribution.csv")
wcvp_name_tb <- fread("./data/raw/wcvp_names.csv")
# select for relevant info 
wcvp_names <- wcvp_name_tb %>% 
  select(plant_name_id, taxon_name, taxon_authors, taxon_status)

# merge names with distribution info.
wcvp_names_w_dist <- merge(wcvp_names, wcvp_distribution_tb, by = "plant_name_id")

na_wcvp_names_w_dist <- wcvp_names_w_dist %>% 
  filter(taxon_name %in% na_plants_alignment$acceptedName) # use acceptedName as this will create the full range for the parentTaxon

# select for relevant info
na_wcvp_names_w_dist <- select(na_wcvp_names_w_dist, plant_name_id, taxon_name, continent, region_code_l2, region, area_code_l3,
                               area, introduced, extinct, location_doubtful)
# determine parentTaxon using the name alignment's criteria 
na_wcvp_parentTaxon_w_dist <- na_wcvp_names_w_dist %>% 
  left_join(na_plants_alignment, by = join_by(taxon_name == acceptedName)) %>% 
  rename(parentTaxon = acceptedNameParent) %>% 
  select(names(na_wcvp_names_w_dist), parentTaxon)

# determine the occurrence's native/introdced/novel statuses 
statuses_for_name <- na_wcvp_parentTaxon_w_dist[parentTaxon == accepted_name, ]

## deal with collasping the subspecific names into the parentTaxon concerning distribution
statuses_for_name <- statuses_for_name %>% 
  distinct(across(-c(plant_name_id, taxon_name)), .keep_all = TRUE) %>% # Take distinct values for parentTaxon rather than the taxon
  group_by(parentTaxon, area_code_l3) %>% # group up values of the parentTaxon and the location
  summarize(across(everything(), ~ first(.)), # Priortize natives vs introduced values per grouping  
            introduced = min(introduced)) %>% # take the 0 value instead of 1 if it exists
  ungroup()


## determine the natives vs introduced for color fill choice
true_native_regions <- statuses_for_name %>% 
  filter(introduced == 0)
true_introduced_regions <- statuses_for_name %>% 
  filter(introduced == 1)

region_fills <- level3_regions %>% 
  mutate(regionStatus = case_when(
    LEVEL3_COD %in% true_native_regions$area_code_l3 ~ "native range", 
    LEVEL3_COD %in% true_introduced_regions$area_code_l3 ~ "introduced range",
    TRUE ~ "undocumented"
  ))

################################################################################################################
thinning_dist <- c(1,5,10) # thinning by 1km, 5km, and 10km
for(i in 1:length(thinning_dist)){
# Now take into account the distance between points....
our_tdi <- 250 # whatever we're declaring as the TDI for this run on outlier analysis 
# read in these names # read in these names TRUE
if(file.exists(paste0("./data/processed/regioned-data/native/", accepted_name_filestyle, ".csv"))){
  
  # load the data
  occur_data <- fread(paste0("./data/processed/regioned-data/native/", accepted_name_filestyle, ".csv"))
} else { # If the name does not exist at this step, include this in the final output table. 
  print(paste("No occur data for", accepted_name))
} # end of if statement

# Produce a spatial thinned version of the dataset to avoid overclustering skewing the outlier analysis
spThin::thin(occur_data, 
             lat.col = "roundedLatitude", 
             long.col = "roundedLongitude", 
             spec.col = "parentTaxon",
             reps = 1,
             write.files = TRUE,
             out.base = "thinned_data",
             out.dir = paste0("./outputs/testing-outlier-removal/", accepted_name_filestyle, "/native-data/250km/", "thinning-", thinning_dist[i], "km"),
             thin.par = thinning_dist[i]
)
occur_data <- fread(paste0("./outputs/testing-outlier-removal/", accepted_name_filestyle, "/native-data/250km/", "thinning-", thinning_dist[i], "km/",  "thinned_data_thin1.csv"))
occur_data <- mutate(occur_data, id = 1:n())
id_dist <- CoordinateCleaner::cc_outl(occur_data, 
                                      lon = "roundedLongitude",
                                      lat = "roundedLatitude", 
                                      species = "parentTaxon", 
                                      method = "distance", # options are "distance", "quantile", and "mad"
                                      mltpl = 5, # The multiplier of the interquantile range (for method = "quantile" or median absolute deviation (for method = "mad)
                                      tdi = our_tdi, # the minimum absolute distance of a record to all other records of a species to be IDed as an outlier in km, change to 250, 500
                                      value = "flagged", # determines type of output being cleaned dataset ("clean") or flagged "flagged"
                                      sampling_thresh = 0, # the cut off threshold for the sampling correction. Indicates the quantile of sampling in which outliers should be ignored. 
                                      verbose = TRUE, 
                                      min_occs = 7, # minimum number of geographically unique datapoints needed for a species to be tested. 
                                      thinning = FALSE, # forces a raster approximation for the distance calc. Useful for species with more than 10,000 records to speed up computation or for small datasets that have very uneven sampling in their datasets. 
                                      thinning_res = 0.5 # the resolution for the spatial thinning in decimal degrees, default = 0.5
)

id_mad <- CoordinateCleaner::cc_outl(occur_data, 
                                     lon = "roundedLongitude",
                                     lat = "roundedLatitude", 
                                     species = "parentTaxon", 
                                     method = "mad", # options are "distance", "quantile", and "mad"
                                     mltpl = 5, # The multiplier of the interquantile range (for method = "quantile" or median absolute deviation (for method = "mad)
                                     tdi = our_tdi, # the minimum absolute distance of a record to all other records of a species to be IDed as an outlier in km
                                     value = "flagged", # determines type of output being cleaned dataset ("clean") or flagged "flagged"
                                     sampling_thresh = 0, # the cut off threshold for the sampling correction. Indicates the quantile of sampling in which outliers should be ignored. 
                                     verbose = TRUE, 
                                     min_occs = 7, # minimum number of geographically unique datapoints needed for a species to be tested. 
                                     thinning = FALSE, # forces a raster approximation for the distance calc. Useful for species with more than 10,000 records to speed up computation or for small datasets that have very uneven sampling in their datasets. 
                                     thinning_res = 0.5 # the resolution for the spatial thinning in decimal degrees, default = 0.5
)

id_quantile <- CoordinateCleaner::cc_outl(occur_data, 
                                          lon = "roundedLongitude",
                                          lat = "roundedLatitude", 
                                          species = "parentTaxon", 
                                          method = "mad", # options are "distance", "quantile", and "mad"
                                          mltpl = 5, # The multiplier of the interquantile range (for method = "quantile" or median absolute deviation (for method = "mad)
                                          tdi = our_tdi, # the minimum absolute distance of a record to all other records of a species to be IDed as an outlier in km
                                          value = "flagged", # determines type of output being cleaned dataset ("clean") or flagged "flagged"
                                          sampling_thresh = 0, # the cut off threshold for the sampling correction. Indicates the quantile of sampling in which outliers should be ignored. 
                                          verbose = TRUE, 
                                          min_occs = 7, # minimum number of geographically unique datapoints needed for a species to be tested. 
                                          thinning = FALSE, # forces a raster approximation for the distance calc. Useful for species with more than 10,000 records to speed up computation or for small datasets that have very uneven sampling in their datasets. 
                                          thinning_res = 0.5 # the resolution for the spatial thinning in decimal degrees, default = 0.5
)
id_dist <- !id_dist # because its annoying to have them ordered opposite...
id_mad <- !id_mad
id_quantile <- !id_quantile
occur_data$distOutlierFlagged <- id_dist # append the outleirs as a boolean field 
occur_data$madOutlierFlagged <- id_mad
occur_data$quantileOutlierFlagged <- id_quantile
random_sample_subset <- filter(random_sample, parentTaxon == accepted_name)
occur_data$quantile <- random_sample_subset$quantile

occur_data_sf <- st_as_sf(occur_data, coords = c("roundedLongitude", "roundedLatitude"), crs = 4326)

# Create convex hull and bounding box
hull <- st_convex_hull(st_union(occur_data_sf))
hull_bbox <- st_bbox(hull)

## Calculate some informative criteria regarding these tests ###########################################
# for dist method
occur_data_sf_buffered<- st_buffer(occur_data_sf, dist = our_tdi*1000) 
# for mad method 
dist <- geosphere::distm(occur_data[, c("roundedLongitude", "roundedLatitude")], 
                         fun = geosphere::distHaversine) / 1000 
dist[dist==0] <- NA # set diagonal to zero 
mins <- apply(dist, 1, mean, na.rm = TRUE) # take minmum dist to any next point
quo <- stats::median(mins, na.rm = TRUE) # obtain the median
tester <- stats::mad(mins, na.rm = TRUE) # obtain the mad stat
mltpl <- 5 # define mltpl
outlier_indices <- which(mins > (quo + tester * mltpl)) # find outlrs
occur_data <- occur_data %>% 
  mutate(is_outlier = id %in% outlier_indices)
cutoff_upper <- quo+tester*mltpl # not gonna bother with the lower as negative distance isnt possible
# for quantile method
# quo <- quantile(mins, c(0.25, 0.75), na.rm = TRUE) # construct quantiles
# iqr <- IQR(mins, na.rm = TRUE) 
# out <- which(mins > (quo[2] + iqr*mltpl)) # flag outliers using the upper bound criterion
# id_quantile <- !seq_along(mins) %in% out # id flags
# occur_data$plotting_quantileOutlierFlagged <- id_quantile # apply
# occur_data$mins <- mins
#######################################################################################################

quantile_id <- unique(occur_data$quantile)

# set-up region fills 
color_mapping <- c("native range" = "darkgreen", "introduced range" = "goldenrod", "undocumented" = "grey")


# Plot spatial data with points marked for removal in red
distance_plot <- ggplot() +
  geom_sf(region_fills, mapping = aes(fill = regionStatus), color = "black") +
  scale_fill_manual(values = color_mapping) +
  geom_sf(data = occur_data_sf, mapping = aes(color = distOutlierFlagged)) +  # Map fill based on whether the point is to be removed
  scale_color_manual(name = "Points Removal", values = c("FALSE" = "black", "TRUE" = "red" )) +  # Use red for points flagged for removal
  coord_sf(xlim = c(hull_bbox[1], hull_bbox[3]), ylim = c(hull_bbox[2], hull_bbox[4])) +  # Zoom to convex hull
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  ggtitle(paste0("Native Distribution for ", accepted_name, " of ", quantile_id, " quantile", "\n using Distance Method with", our_tdi, "km outlier threshold"))

buffer_fills <- c("green", "red")
extra_mapping <- c("FALSE"="green", "TRUE" = "red", color_mapping)
distance_w_buffer_plot <- ggplot() +
  geom_sf(region_fills, mapping = aes(fill = regionStatus), color = "black") +
  geom_sf(data = occur_data_sf_buffered, aes(fill = distOutlierFlagged), alpha = 0.2) +  # Use aes for mapping colors
  geom_sf(data = occur_data_sf, aes(color = distOutlierFlagged)) +  # Map color for points flagged for removal'
  scale_fill_manual(values = color_mapping) +
  #scale_fill_manual(name = "Points Removal", values = c("FALSE" = "green", "TRUE" = "red")) +  # Fill for buffer circles
  scale_fill_manual(name = "region and point status", values = extra_mapping) +
  scale_color_manual(name = "Points Removal", values = c("FALSE" = "black", "TRUE" = "red")) +  # Color for points
  coord_sf(xlim = c(hull_bbox[1], hull_bbox[3]), ylim = c(hull_bbox[2], hull_bbox[4])) +  # Zoom to convex hull
  theme_minimal() +
  ggtitle(paste0("Native Distribution for ", accepted_name, " of ", quantile_id, " quantile", "\n using Distance Method with", our_tdi, "km outlier threshold"))

ggsave(plot = distance_plot, filename = paste0("./outputs/testing-outlier-removal/", accepted_name_filestyle, "/native-data/250km/", accepted_name_filestyle, "-native-data-distribution-distance-method-plot.jpg"), height = 10, width = 10)

mad_plot <- ggplot() +
  geom_sf(region_fills, mapping = aes(fill = regionStatus), color = "black") +
  scale_fill_manual(values = color_mapping) +
  geom_sf(data = occur_data_sf, mapping = aes(color = madOutlierFlagged)) +  # Map fill based on whether the point is to be removed
  scale_color_manual(name = "Points Removal", values = c("FALSE" = "black", "TRUE" = "red" )) +  # Use red for points flagged for removal
  coord_sf(xlim = c(hull_bbox[1], hull_bbox[3]), ylim = c(hull_bbox[2], hull_bbox[4])) +  # Zoom to convex hull
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  ggtitle(paste0("Native Distribution for ", accepted_name, " of ", quantile_id, " quantile", "\n using Mean Absolute Distance Method with", our_tdi, "km outlier threshold"))
# Plot histogram of minimum distances
mad_hist_plot <- ggplot(data = data.frame(min_distance = mins), aes(x = min_distance)) +
  geom_histogram(binwidth = 10, color = "black", fill = "lightblue") +
  labs(title = paste("Histogram of Minimum Distances for", accepted_name, "Mean Absolute Distance Method"),
       x = "Minimum Distance (km)",
       y = "Frequency") +
  geom_vline(xintercept = cutoff_upper, linetype = "dashed", color = "red") +
  annotate("text", x = cutoff_upper, y = 1, label = paste("Upper Cutoff: ", round(cutoff_upper, 2)), 
           vjust = -1, color = "red") +
  theme_minimal()

ggsave(plot = mad_plot, filename = paste0("./outputs/testing-outlier-removal/", accepted_name_filestyle, "/native-data/250km/", accepted_name_filestyle, "-native-data-distribution-mad-method-plot.jpg"), height = 10, width = 10)
# quantile_plot <- ggplot() +
#   geom_sf(region_fills, mapping = aes(fill = regionStatus), color = "black") +
#   scale_fill_manual(values = color_mapping) +
#   geom_sf(data = occur_data_sf, mapping = aes(color = quantileOutlierFlagged)) +  # Map fill based on whether the point is to be removed
#   scale_color_manual(name = "Points Removal", values = c("FALSE" = "black", "TRUE" = "red" )) +  # Use red for points flagged for removal
#   coord_sf(xlim = c(hull_bbox[1], hull_bbox[3]), ylim = c(hull_bbox[2], hull_bbox[4])) +  # Zoom to convex hull
#   theme_minimal() + 
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
#   ggtitle(paste0("Native Distribution for ", accepted_name, " of ", quantile_id, " quantile", "\n using Quantiles Method with", our_tdi, "km outlier threshold"))
# ggsave(plot = quantile_plot, filename = paste0("./outputs/testing-outlier-removal/", accepted_name_filestyle, "/native-data/", accepted_name_filestyle, "-native-data-distribution-quantile-method-plot.jpg"), height = 10, width = 10)
# 
# boxplot_quantile_plot <- ggplot(data = occur_data, aes(y = mins)) +
#   geom_boxplot(outlier.shape = 16, outlier.size = 2, color = "black") +
#   labs(title = "Boxplot of Minimum Distances",
#        y = "Minimum Distance (km)") +
#   theme_minimal() +
#   geom_hline(yintercept = quo[2] + iqr * mltpl, linetype = "dashed", color = "red") +
#   annotate("text", x = 0.5, y = quo[2] + iqr * mltpl + 1, 
#            label = paste("Upper Cutoff: ", round(quo[2] + iqr * mltpl, 2), "km"), 
#            vjust = -0.5, hjust = 1.5, color = "red")
# 

native_outlier_removal_summary <- data.frame(
  nativeData = nrow(occur_data), 
  distOutlierRemoval = nrow(subset(occur_data, occur_data$distOutlierFlagged == TRUE)),
  madOutlierRemoval = nrow(subset(occur_data, occur_data$madOutlierFlagged == TRUE))#, 
  #quantileOutlierRemoval = nrow(subset(occur_data, occur_data$quantileOutlierFlagged == TRUE))
)

# create a pdf output
pdf(file = paste0("./outputs/outlier-pdfs/native-data/250km/", accepted_name_filestyle),  width = 30, height = 30)
grid.arrange(distance_plot, distance_w_buffer_plot,  mad_plot, mad_hist_plot, ncol = 1)
dev.off()




fwrite(native_outlier_removal_summary, paste0("./outputs/testing-outlier-removal/", accepted_name_filestyle, "/native-data/250km/", accepted_name_filestyle, "-native-data-outlier-removal-summary.csv"))

#### Repeat for ALL records ######################################################################################
if(file.exists(paste0("./data/processed/regioned-data/all/", accepted_name_filestyle, ".csv"))){
  
  # load the data
  occur_data <- fread(paste0("./data/processed/regioned-data/all/", accepted_name_filestyle, ".csv"))
} else { # If the name does not exist at this step, include this in the final output table. 
  print(paste("No occur data for", accepted_name))
} # end of if statement

# Produce a spatial thinned version of the dataset to avoid overclustering skewing the outlier analysis
spThin::thin(occur_data, 
             lat.col = "roundedLatitude", 
             long.col = "roundedLongitude", 
             spec.col = "parentTaxon",
             reps = 1,
             write.files = TRUE,
             out.base = "thinned_data",
             out.dir = paste0("./outputs/testing-outlier-removal/", accepted_name_filestyle, "/all-data/250km/", "thinning-", thinning_dist[i], "km"),
             thin.par = thinning_dist[i]
)
occur_data <- fread(paste0("./outputs/testing-outlier-removal/", accepted_name_filestyle, "/all-data/250km/", , "thinning-", thinning_dist[i], "km/",  "thinned_data_thin1.csv"))
occur_data <- mutate(occur_data, id = 1:n())

id_dist <- CoordinateCleaner::cc_outl(occur_data, 
                                      lon = "roundedLongitude",
                                      lat = "roundedLatitude", 
                                      species = "parentTaxon", 
                                      method = "distance", # options are "distance", "quantile", and "mad"
                                      mltpl = 5, # The multiplier of the interquantile range (for method = "quantile" or median absolute deviation (for method = "mad)
                                      tdi = our_tdi, # the minimum absolute distance of a record to all other records of a species to be IDed as an outlier in km
                                      value = "flagged", # determines type of output being cleaned dataset ("clean") or flagged "flagged"
                                      sampling_thresh = 0, # the cut off threshold for the sampling correction. Indicates the quantile of sampling in which outliers should be ignored. 
                                      verbose = TRUE, 
                                      min_occs = 7, # minimum number of geographically unique datapoints needed for a species to be tested. 
                                      thinning = FALSE, # forces a raster approximation for the distance calc. Useful for species with more than 10,000 records to speed up computation or for small datasets that have very uneven sampling in their datasets. 
                                      thinning_res = 0.5 # the resolution for the spatial thinning in decimal degrees, default = 0.5
)

id_mad <- CoordinateCleaner::cc_outl(occur_data, 
                                     lon = "roundedLongitude",
                                     lat = "roundedLatitude", 
                                     species = "parentTaxon", 
                                     method = "mad", # options are "distance", "quantile", and "mad"
                                     mltpl = 5, # The multiplier of the interquantile range (for method = "quantile" or median absolute deviation (for method = "mad)
                                     tdi = our_tdi, # the minimum absolute distance of a record to all other records of a species to be IDed as an outlier in km
                                     value = "flagged", # determines type of output being cleaned dataset ("clean") or flagged "flagged"
                                     sampling_thresh = 0, # the cut off threshold for the sampling correction. Indicates the quantile of sampling in which outliers should be ignored. 
                                     verbose = TRUE, 
                                     min_occs = 7, # minimum number of geographically unique datapoints needed for a species to be tested. 
                                     thinning = FALSE, # forces a raster approximation for the distance calc. Useful for species with more than 10,000 records to speed up computation or for small datasets that have very uneven sampling in their datasets. 
                                     thinning_res = 0.5 # the resolution for the spatial thinning in decimal degrees, default = 0.5
)

id_quantile <- CoordinateCleaner::cc_outl(occur_data, 
                                          lon = "roundedLongitude",
                                          lat = "roundedLatitude", 
                                          species = "parentTaxon", 
                                          method = "mad", # options are "distance", "quantile", and "mad"
                                          mltpl = 5, # The multiplier of the interquantile range (for method = "quantile" or median absolute deviation (for method = "mad)
                                          tdi = our_tdi, # the minimum absolute distance of a record to all other records of a species to be IDed as an outlier in km
                                          value = "flagged", # determines type of output being cleaned dataset ("clean") or flagged "flagged"
                                          sampling_thresh = 0, # the cut off threshold for the sampling correction. Indicates the quantile of sampling in which outliers should be ignored. 
                                          verbose = TRUE, 
                                          min_occs = 7, # minimum number of geographically unique datapoints needed for a species to be tested. 
                                          thinning = FALSE, # forces a raster approximation for the distance calc. Useful for species with more than 10,000 records to speed up computation or for small datasets that have very uneven sampling in their datasets. 
                                          thinning_res = 0.5 # the resolution for the spatial thinning in decimal degrees, default = 0.5
)
id_dist <- !id_dist # because its annoying to have them ordered opposite...
id_mad <- !id_mad
id_quantile <- !id_quantile
occur_data$distOutlierFlagged <- id_dist # append the outleirs as a boolean field 
occur_data$madOutlierFlagged <- id_mad
occur_data$quantileOutlierFlagged <- id_quantile

random_sample_subset <- filter(random_sample, parentTaxon == accepted_name)
occur_data$quantile <- random_sample_subset$quantile

occur_data_sf <- st_as_sf(occur_data, coords = c("roundedLongitude", "roundedLatitude"), crs = 4326)

# Create convex hull and bounding box
hull <- st_convex_hull(st_union(occur_data_sf))
hull_bbox <- st_bbox(hull)

## Calculate some informative criteria regarding these tests ###########################################
# for dist method
occur_data_sf_buffered<- st_buffer(occur_data_sf, dist = our_tdi*1000) 
# for mad method 
dist <- geosphere::distm(occur_data[, c("roundedLongitude", "roundedLatitude")], 
                         fun = geosphere::distHaversine) / 1000 
dist[dist==0] <- NA # set diagonal to zero 
mins <- apply(dist, 1, mean, na.rm = TRUE) # take minmum dist to any next point
quo <- stats::median(mins, na.rm = TRUE) # obtain the median
tester <- stats::mad(mins, na.rm = TRUE) # obtain the mad stat
mltpl <- 5 # define mltpl
outlier_indices <- which(mins > (quo + tester * mltpl)) # find outlrs
occur_data <- occur_data %>% 
  mutate(is_outlier = id %in% outlier_indices)
cutoff_upper <- quo+tester*mltpl # not gonna bother with the lower as negative distance isnt possible
# for quantile method
# quo <- quantile(mins, c(0.25, 0.75), na.rm = TRUE) # construct quantiles
# iqr <- IQR(mins, na.rm = TRUE) 
# out <- which(mins > (quo[2] + iqr*mltpl)) # flag outliers using the upper bound criterion
# id_quantile <- !seq_along(mins) %in% out # id flags
# occur_data$plotting_quantileOutlierFlagged <- id_quantile # apply
# occur_data$mins <- mins
#######################################################################################################

# Must be turned off for map subsetting to operate correctly 
if(hull_bbox[1] == 0 & hull_bbox[2] == -90 & hull_bbox[3] == 0 & hull_bbox[4] == -90){
  sf_use_s2(FALSE)
  hull_bbox <- st_bbox(st_convex_hull(level3_regions))
}

sf_use_s2(TRUE)


quantile_id <- unique(occur_data$quantile)

# Plot spatial data with points marked for removal in red
distance_plot <- ggplot() +
  geom_sf(region_fills, mapping = aes(fill = regionStatus), color = "black") +
  scale_fill_manual(values = color_mapping) +
  geom_sf(data = occur_data_sf, mapping = aes(color = distOutlierFlagged)) +  # Map fill based on whether the point is to be removed
  scale_color_manual(name = "Points Removal", values = c("FALSE" = "black", "TRUE" = "red" )) +  # Use red for points flagged for removal
  coord_sf(xlim = c(hull_bbox[1], hull_bbox[3]), ylim = c(hull_bbox[2], hull_bbox[4])) +  # Zoom to convex hull
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  ggtitle(paste0("All Distribution for ", accepted_name, " of ", quantile_id, " quantile", "\n using Distance Method with", our_tdi, "km outlier threshold"))
ggsave(plot = distance_plot, filename = paste0("./outputs/testing-outlier-removal/", accepted_name_filestyle, "/all-data/250km/", accepted_name_filestyle, "-all-data-distribution-distance-method-plot.jpg"), height = 10, width = 10)

buffer_fills <- c("green", "red")
extra_mapping <- c("FALSE"="green", "TRUE" = "red", color_mapping)
distance_w_buffer_plot <- ggplot() +
  geom_sf(region_fills, mapping = aes(fill = regionStatus), color = "black") +
  geom_sf(data = occur_data_sf_buffered, aes(fill = distOutlierFlagged), alpha = 0.2) +  # Use aes for mapping colors
  geom_sf(data = occur_data_sf, aes(color = distOutlierFlagged)) +  # Map color for points flagged for removal'
  scale_fill_manual(values = color_mapping) +
  #scale_fill_manual(name = "Points Removal", values = c("FALSE" = "green", "TRUE" = "red")) +  # Fill for buffer circles
  scale_fill_manual(name = "region and point status", values = extra_mapping) +
  scale_color_manual(name = "Points Removal", values = c("FALSE" = "black", "TRUE" = "red")) +  # Color for points
  coord_sf(xlim = c(hull_bbox[1], hull_bbox[3]), ylim = c(hull_bbox[2], hull_bbox[4])) +  # Zoom to convex hull
  theme_minimal() +
  ggtitle(paste0("All Distribution for ", accepted_name, " of ", quantile_id, " quantile", "\n using Distance Method with", our_tdi, "km outlier threshold"))





mad_plot <- ggplot() +
  geom_sf(region_fills, mapping = aes(fill = regionStatus), color = "black") +
  scale_fill_manual(values = color_mapping) +
  geom_sf(data = occur_data_sf, mapping = aes(color = madOutlierFlagged)) +  # Map fill based on whether the point is to be removed
  scale_color_manual(name = "Points Removal", values = c("FALSE" = "black", "TRUE" = "red" )) +  # Use red for points flagged for removal
  coord_sf(xlim = c(hull_bbox[1], hull_bbox[3]), ylim = c(hull_bbox[2], hull_bbox[4])) +  # Zoom to convex hull
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  ggtitle(paste0("All Distribution for ", accepted_name, " of ", quantile_id, " quantile", "\n using Mean Absolute Distance Method with", our_tdi, "km outlier threshold"))
ggsave(plot = mad_plot, filename = paste0("./outputs/testing-outlier-removal/", accepted_name_filestyle, "/all-data/250km/", accepted_name_filestyle, "-all-data-distribution-mad-method-plot.jpg"), height = 10, width = 10)

# Plot histogram of minimum distances
mad_hist_plot <- ggplot(data = data.frame(min_distance = mins), aes(x = min_distance)) +
  geom_histogram(binwidth = 10, color = "black", fill = "lightblue") +
  labs(title = paste("Histogram of Minimum Distances for", accepted_name, "Mean Absolute Distance Method"),
       x = "Minimum Distance (km)",
       y = "Frequency") +
  geom_vline(xintercept = cutoff_upper, linetype = "dashed", color = "red") +
  annotate("text", x = cutoff_upper, y = 1, label = paste("Upper Cutoff: ", round(cutoff_upper, 2)), 
           vjust = -1, color = "red") +
  theme_minimal()


# quantile_plot <- ggplot() +
#   geom_sf(region_fills, mapping = aes(fill = regionStatus), color = "black") +
#   scale_fill_manual(values = color_mapping) +
#   geom_sf(data = occur_data_sf, mapping = aes(color = quantileOutlierFlagged)) +  # Map fill based on whether the point is to be removed
#   scale_color_manual(name = "Points Removal", values = c("FALSE" = "black", "TRUE" = "red" )) +  # Use red for points flagged for removal
#   coord_sf(xlim = c(hull_bbox[1], hull_bbox[3]), ylim = c(hull_bbox[2], hull_bbox[4])) +  # Zoom to convex hull
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
#   ggtitle(paste0("All Distribution for ", accepted_name, " of ", quantile_id, " quantile", "\n using Quantiles Method with", our_tdi, "km outlier threshold"))
# ggsave(plot = quantile_plot, filename = paste0("./outputs/testing-outlier-removal/", accepted_name_filestyle, "/all-data/", accepted_name_filestyle, "-all-data-distribution-quantile-method-plot.jpg"), height = 10, width = 10)
all_outlier_removal_summary <- data.frame(
  allData = nrow(occur_data), 
  distOutlierRemoval = nrow(subset(occur_data, occur_data$distOutlierFlagged == TRUE)),
  madOutlierRemoval = nrow(subset(occur_data, occur_data$madOutlierFlagged == TRUE))#, 
  #quantileOutlierRemoval = nrow(subset(occur_data, occur_data$quantileOutlierFlagged == TRUE))
)

# create a pdf output
pdf(file = paste0("./outputs/outlier-pdfs/all-data/250km/", accepted_name_filestyle),  width = 30, height = 30)
grid.arrange(distance_plot, distance_w_buffer_plot,  mad_plot, mad_hist_plot, ncol = 1)
dev.off()

fwrite(all_outlier_removal_summary, paste0("./outputs/testing-outlier-removal/", accepted_name_filestyle, "/all-data/250km/", accepted_name_filestyle, "-all-data-outlier-removal-summary.csv"))

#####################################################################################################################

# Using native + introduced.
# read in these names # read in these names TRUE
if(file.exists(paste0("./data/processed/regioned-data/native/", accepted_name_filestyle, ".csv"))){
  
  # load the data
  native_occur_data <- fread(paste0("./data/processed/regioned-data/native/", accepted_name_filestyle, ".csv"))
} else { # If the name does not exist at this step, include this in the final output table. 
  print(paste("No occur data for", accepted_name))
} # end of if statement

if(file.exists(paste0("./data/processed/regioned-data/introduced/", accepted_name_filestyle, ".csv"))){
  
  # load the data
  introduced_occur_data <- fread(paste0("./data/processed/regioned-data/introduced/", accepted_name_filestyle, ".csv"))
} else { # If the name does not exist at this step, include this in the final output table. 
  print(paste("No occur data for", accepted_name))
} # end of if statement

# Create a data list 
data_list <- list(native_occur_data, introduced_occur_data)
# Filter out anything that contains nothing.
non_null_data_list <- Filter(Negate(is.null), data_list)
# Combine all non null dataframes
if(length(non_null_data_list) > 0){
  occur_data <- do.call(rbind, non_null_data_list)
} else { 
  occur_data <- NULL
}

# Produce a spatial thinned version of the dataset to avoid overclustering skewing the outlier analysis
spThin::thin(occur_data, 
             lat.col = "roundedLatitude", 
             long.col = "roundedLongitude", 
             spec.col = "parentTaxon",
             reps = 1,
             write.files = TRUE,
             out.base = "thinned_data",
             out.dir = paste0("./outputs/testing-outlier-removal/", accepted_name_filestyle, "/native-plus-introduced-data/250km/", "thinning-", thinning_dist[i], "km"),
             thin.par = thinning_dist[i]
)
occur_data <- fread(paste0("./outputs/testing-outlier-removal/", accepted_name_filestyle, "/native-plus-introduced-data/250km/", , "thinning-", thinning_dist[i], "km/",  "thinned_data_thin1.csv"))

occur_data <- mutate(occur_data, id = 1:n())
id_dist <- CoordinateCleaner::cc_outl(occur_data, 
                                      lon = "roundedLongitude",
                                      lat = "roundedLatitude", 
                                      species = "parentTaxon", 
                                      method = "distance", # options are "distance", "quantile", and "mad"
                                      mltpl = 5, # The multiplier of the interquantile range (for method = "quantile" or median absolute deviation (for method = "mad)
                                      tdi = our_tdi, # the minimum absolute distance of a record to all other records of a species to be IDed as an outlier in km
                                      value = "flagged", # determines type of output being cleaned dataset ("clean") or flagged "flagged"
                                      sampling_thresh = 0, # the cut off threshold for the sampling correction. Indicates the quantile of sampling in which outliers should be ignored. 
                                      verbose = TRUE, 
                                      min_occs = 7, # minimum number of geographically unique datapoints needed for a species to be tested. 
                                      thinning = FALSE, # forces a raster approximation for the distance calc. Useful for species with more than 10,000 records to speed up computation or for small datasets that have very uneven sampling in their datasets. 
                                      thinning_res = 0.5 # the resolution for the spatial thinning in decimal degrees, default = 0.5
)

id_mad <- CoordinateCleaner::cc_outl(occur_data, 
                                     lon = "roundedLongitude",
                                     lat = "roundedLatitude", 
                                     species = "parentTaxon", 
                                     method = "mad", # options are "distance", "quantile", and "mad"
                                     mltpl = 5, # The multiplier of the interquantile range (for method = "quantile" or median absolute deviation (for method = "mad)
                                     tdi = our_tdi, # the minimum absolute distance of a record to all other records of a species to be IDed as an outlier in km
                                     value = "flagged", # determines type of output being cleaned dataset ("clean") or flagged "flagged"
                                     sampling_thresh = 0, # the cut off threshold for the sampling correction. Indicates the quantile of sampling in which outliers should be ignored. 
                                     verbose = TRUE, 
                                     min_occs = 7, # minimum number of geographically unique datapoints needed for a species to be tested. 
                                     thinning = FALSE, # forces a raster approximation for the distance calc. Useful for species with more than 10,000 records to speed up computation or for small datasets that have very uneven sampling in their datasets. 
                                     thinning_res = 0.5 # the resolution for the spatial thinning in decimal degrees, default = 0.5
)

id_quantile <- CoordinateCleaner::cc_outl(occur_data, 
                                          lon = "roundedLongitude",
                                          lat = "roundedLatitude", 
                                          species = "parentTaxon", 
                                          method = "mad", # options are "distance", "quantile", and "mad"
                                          mltpl = 5, # The multiplier of the interquantile range (for method = "quantile" or median absolute deviation (for method = "mad)
                                          tdi = our_tdi, # the minimum absolute distance of a record to all other records of a species to be IDed as an outlier in km
                                          value = "flagged", # determines type of output being cleaned dataset ("clean") or flagged "flagged"
                                          sampling_thresh = 0, # the cut off threshold for the sampling correction. Indicates the quantile of sampling in which outliers should be ignored. 
                                          verbose = TRUE, 
                                          min_occs = 7, # minimum number of geographically unique datapoints needed for a species to be tested. 
                                          thinning = FALSE, # forces a raster approximation for the distance calc. Useful for species with more than 10,000 records to speed up computation or for small datasets that have very uneven sampling in their datasets. 
                                          thinning_res = 0.5 # the resolution for the spatial thinning in decimal degrees, default = 0.5
)
id_dist <- !id_dist # because its annoying to have them ordered opposite...
id_mad <- !id_mad
id_quantile <- !id_quantile
occur_data$distOutlierFlagged <- id_dist # append the outleirs as a boolean field 
occur_data$madOutlierFlagged <- id_mad
occur_data$quantileOutlierFlagged <- id_quantile

random_sample_subset <- filter(random_sample, parentTaxon == accepted_name)
occur_data$quantile <- random_sample_subset$quantile

occur_data_sf <- st_as_sf(occur_data, coords = c("roundedLongitude", "roundedLatitude"), crs = 4326)

## Calculate some informative criteria regarding these tests ###########################################
# for dist method
occur_data_sf_buffered<- st_buffer(occur_data_sf, dist = our_tdi*1000) 
# for mad method 
dist <- geosphere::distm(occur_data[, c("roundedLongitude", "roundedLatitude")], 
                         fun = geosphere::distHaversine) / 1000 
dist[dist==0] <- NA # set diagonal to zero 
mins <- apply(dist, 1, mean, na.rm = TRUE) # take minmum dist to any next point
quo <- stats::median(mins, na.rm = TRUE) # obtain the median
tester <- stats::mad(mins, na.rm = TRUE) # obtain the mad stat
mltpl <- 5 # define mltpl
outlier_indices <- which(mins > (quo + tester * mltpl)) # find outlrs
occur_data <- occur_data %>% 
  mutate(is_outlier = id %in% outlier_indices)
cutoff_upper <- quo+tester*mltpl # not gonna bother with the lower as negative distance isnt possible
# for quantile method
# quo <- quantile(mins, c(0.25, 0.75), na.rm = TRUE) # construct quantiles
# iqr <- IQR(mins, na.rm = TRUE) 
# out <- which(mins > (quo[2] + iqr*mltpl)) # flag outliers using the upper bound criterion
# id_quantile <- !seq_along(mins) %in% out # id flags
# occur_data$plotting_quantileOutlierFlagged <- id_quantile # apply
# occur_data$mins <- mins
#######################################################################################################


# Create convex hull and bounding box
hull <- st_convex_hull(st_union(occur_data_sf))
hull_bbox <- st_bbox(hull)

quantile_id <- unique(occur_data$quantile)

# set-up region fills 
color_mapping <- c("native range" = "darkgreen", "introduced range" = "goldenrod", "undocumented" = "grey")


# Plot spatial data with points marked for removal in red
distance_plot <- ggplot() +
  geom_sf(region_fills, mapping = aes(fill = regionStatus), color = "black") +
  scale_fill_manual(values = color_mapping) +
  geom_sf(data = occur_data_sf, mapping = aes(color = distOutlierFlagged)) +  # Map fill based on whether the point is to be removed
  scale_color_manual(name = "Points Removal", values = c("FALSE" = "black", "TRUE" = "red" )) +  # Use red for points flagged for removal
  coord_sf(xlim = c(hull_bbox[1], hull_bbox[3]), ylim = c(hull_bbox[2], hull_bbox[4])) +  # Zoom to convex hull
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  ggtitle(paste0("Native and Introduced Distribution for ", accepted_name, " of ", quantile_id, " quantile", "\n using Distance Method with", our_tdi, "km outlier threshold"))
ggsave(plot = distance_plot, filename = paste0("./outputs/testing-outlier-removal/", accepted_name_filestyle, "/native-plus-introduced-data/250km/", accepted_name_filestyle, "-native-plus-introduced-data-distribution-distance-method-plot.jpg"), height = 10, width = 10)

buffer_fills <- c("green", "red")
extra_mapping <- c("FALSE"="green", "TRUE" = "red", color_mapping)
distance_w_buffer_plot <- ggplot() +
  geom_sf(region_fills, mapping = aes(fill = regionStatus), color = "black") +
  geom_sf(data = occur_data_sf_buffered, aes(fill = distOutlierFlagged), alpha = 0.2) +  # Use aes for mapping colors
  geom_sf(data = occur_data_sf, aes(color = distOutlierFlagged)) +  # Map color for points flagged for removal'
  scale_fill_manual(values = color_mapping) +
  #scale_fill_manual(name = "Points Removal", values = c("FALSE" = "green", "TRUE" = "red")) +  # Fill for buffer circles
  scale_fill_manual(name = "region and point status", values = extra_mapping) +
  scale_color_manual(name = "Points Removal", values = c("FALSE" = "black", "TRUE" = "red")) +  # Color for points
  coord_sf(xlim = c(hull_bbox[1], hull_bbox[3]), ylim = c(hull_bbox[2], hull_bbox[4])) +  # Zoom to convex hull
  theme_minimal() +
  ggtitle(paste0("Native & Introduced Distribution for ", accepted_name, " of ", quantile_id, " quantile", "\n using Distance Method with", our_tdi, "km outlier threshold"))



mad_plot <- ggplot() +
  geom_sf(region_fills, mapping = aes(fill = regionStatus), color = "black") +
  scale_fill_manual(values = color_mapping) +
  geom_sf(data = occur_data_sf, mapping = aes(color = madOutlierFlagged)) +  # Map fill based on whether the point is to be removed
  scale_color_manual(name = "Points Removal", values = c("FALSE" = "black", "TRUE" = "red" )) +  # Use red for points flagged for removal
  coord_sf(xlim = c(hull_bbox[1], hull_bbox[3]), ylim = c(hull_bbox[2], hull_bbox[4])) +  # Zoom to convex hull
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  ggtitle(paste0("Native and Introduced Distribution for ", accepted_name, " of ", quantile_id, " quantile", "\n using Mean Absolute Distance Method with", our_tdi, "km outlier threshold"))
ggsave(plot = mad_plot, filename = paste0("./outputs/testing-outlier-removal/", accepted_name_filestyle, "/native-plus-introduced-data/250km/", accepted_name_filestyle, "-native-plus-introduced-data-distribution-mad-method-plot.jpg"), height = 10, width = 10)

# Plot histogram of minimum distances
mad_hist_plot <- ggplot(data = data.frame(min_distance = mins), aes(x = min_distance)) +
  geom_histogram(binwidth = 10, color = "black", fill = "lightblue") +
  labs(title = paste("Histogram of Minimum Distances for", accepted_name, "Mean Absolute Distance Method"),
       x = "Minimum Distance (km)",
       y = "Frequency") +
  geom_vline(xintercept = cutoff_upper, linetype = "dashed", color = "red") +
  annotate("text", x = cutoff_upper, y = 1, label = paste("Upper Cutoff: ", round(cutoff_upper, 2)), 
           vjust = -1, color = "red") +
  theme_minimal()


# quantile_plot <- ggplot() +
#   geom_sf(region_fills, mapping = aes(fill = regionStatus), color = "black") +
#   scale_fill_manual(values = color_mapping) +
#   geom_sf(data = occur_data_sf, mapping = aes(color = quantileOutlierFlagged)) +  # Map fill based on whether the point is to be removed
#   scale_color_manual(name = "Points Removal", values = c("FALSE" = "black", "TRUE" = "red" )) +  # Use red for points flagged for removal
#   coord_sf(xlim = c(hull_bbox[1], hull_bbox[3]), ylim = c(hull_bbox[2], hull_bbox[4])) +  # Zoom to convex hull
#   theme_minimal() + 
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
#   ggtitle(paste0("Native and Introduced Distribution for ", accepted_name, " of ", quantile_id, " quantile", "\n using Quantiles Method with", our_tdi, "km outlier threshold"))
# ggsave(plot = quantile_plot, filename = paste0("./outputs/testing-outlier-removal/", accepted_name_filestyle, "/native-plus-introduced-data/", accepted_name_filestyle, "-native-plus-introduced-data-distribution-quantile-method-plot.jpg"), height = 10, width = 10)
native_and_introduced_outlier_removal_summary <- data.frame(
  nativeAndIntroducedData = nrow(occur_data), 
  distOutlierRemoval = nrow(subset(occur_data, occur_data$distOutlierFlagged == TRUE)),
  madOutlierRemoval = nrow(subset(occur_data, occur_data$madOutlierFlagged == TRUE))#, 
  #quantileOutlierRemoval = nrow(subset(occur_data, occur_data$quantileOutlierFlagged == TRUE))
)

# create a pdf output
pdf(file = paste0("./outputs/outlier-pdfs/native-and-introduced-data/250km/", accepted_name_filestyle),  width = 30, height = 30)
grid.arrange(distance_plot, distance_w_buffer_plot,  mad_plot, mad_hist_plot, ncol = 1)
dev.off()
fwrite(native_and_introduced_outlier_removal_summary, paste0("./outputs/testing-outlier-removal/", accepted_name_filestyle, "/native-plus-introduced-data/250km/", accepted_name_filestyle, "-native-plus-introduced-outlier-removal-summary.csv"))
}