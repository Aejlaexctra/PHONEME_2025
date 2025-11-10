# --- SET RANDOM SEED ---
set.seed(42)

# --- (1) DOWNLOADING DATASETS AND MODULES --- ####

## --- Import necessary packages --- ####
library(nloptr)
library(nlme)
library(parallel)
library(ggplot2)
library(dplyr)
library(tidyr)
library(maps)

## --- Loading Phylogenetic, Spatial Distance and Contact Matrices --- ####
phylomatrix <- read.csv("data/distance_matrices/phylogenetic_covariance_matrix.csv",
                        check.names = FALSE)
rownames(phylomatrix) <- phylomatrix[,1]
phylomatrix <- phylomatrix[,-1]
phylomatrix <- as.matrix(phylomatrix)
spmatrix <- read.csv("data/distance_matrices/spatial_distance_matrix.csv",
                     check.names = FALSE)
rownames(spmatrix) <- spmatrix[,1]
spmatrix <- spmatrix[,-1]
spmatrix <- as.matrix(spmatrix)

## --- Loading NEE22 (3) --- ####
nee22 <- read.csv("data/Global predictors of language endangerment and the future of linguistic diversity Data 2.csv",
                  stringsAsFactors = FALSE)
colnames(nee22)[colnames(nee22) == "island"] <- "Island.Endemic"
# Convert Documention variable from string to integer
nee22$documentation[nee22$documentation == "little or none"] <- "0"
nee22$documentation[nee22$documentation == "basic"] <- "1"
nee22$documentation[nee22$documentation == "detailed"] <- "2"
nee22$documentation <- as.integer(nee22$documentation)

## --- Loading NEE24 (1) --- ####
nee24 <- read.csv("data/Islands are engines of language diversity data.csv",
                  stringsAsFactors = FALSE)
colnames(nee24)[colnames(nee24) == "L1.Population"] <- "L1_pop"
nee24$Island.Endemic <- as.numeric(nee24$Island.Endemic)

## --- Loading Phoible (4) --- ####
var_phon_inv_spec_data <- read.csv("Data/cldf-datasets-inventory-study/PHOIBLE-data.csv")
var_phon_inv_glot_data <- read.csv("Data/cldf-datasets-inventory-study/phoible/cldf/languages.csv")
phoible <- merge(
  var_phon_inv_spec_data[, c("Glottocode", "Sounds", "Latitude", "Longitude")],
  var_phon_inv_glot_data[, c("ISO639P3code", "Glottocode")], by.x = "Glottocode", by.y = "Glottocode", all.x = TRUE)
colnames(phoible)[colnames(phoible) == "Sounds"] <- "Phoneme.Inventory.Size"

# --- (2) FUNCTION DEFINITIONS --- ####

## --- GLS (2) Setup --- ####
best_p <- function (p,formula,data,spmatrix,phylomatrix) {
  spmatrix <- spmatrix/max(spmatrix)
  spmatrix <- exp(-(spmatrix/p[2])^2)
  mat <- as.matrix((p[3]*(1-p[1])*spmatrix+(1-p[1])*(1-p[3])*phylomatrix+p[1]*diag(dim(phylomatrix)[1])))
  res <- try(gls(model=formula,data=data,correlation=corSymm(mat[lower.tri(mat)],fixed=T),method="ML"),silent=T)
  if (inherits(res,"try-error")) {
    out <- -10000
  } else {
    out <- res$logLik
    # if (res$logLik>0) {out <- -10000}
  }
  -out
}

ml_fit <- function (p,formula,data,spmatrix,phylomatrix) {
  spmatrix <- spmatrix/max(spmatrix)
  spmatrix <- exp(-(spmatrix/p[2])^2)
  mat <- as.matrix((p[3]*(1-p[1])*spmatrix+(1-p[1])*(1-p[3])*phylomatrix+p[1]*diag(dim(phylomatrix)[1])))
  res <- try(gls(model=formula,data=data,correlation=corSymm(mat[lower.tri(mat)],fixed=T),method="ML"),silent=T)
  res
}

# --- (3) CLEANING AND FORMATING DATASETS --- ####

## --- A1a --- ####

# --- Merging Datasets A1a ---
A1a <- merge(
  nee22[, 
        c("ISO", "L1_pop", "region", "Island.Endemic")],
  phoible[, c("ISO639P3code", "Glottocode", "Latitude", "Longitude", "Phoneme.Inventory.Size")], by.x = "ISO", by.y = "ISO639P3code", all.x = TRUE)

# --- Data Cleanup A1a ---
# Remove all NA entries
print("Size of initial dataset:")
print(dim(A1a)[1])
old_A1a <- A1a
A1a <- na.omit(A1a) 
print("Size of dataset with NA removed:")
print(dim(A1a)[1])
# Remove iso duplicate entries
print("Duplicate iso entries:")
iso_duplicates <- which(table(A1a$ISO) != 1)
print(which(table(A1a$ISO) != 1))
A1a <- A1a[!(A1a$ISO %in% names(which(table(A1a$ISO) != 1))),]
print("Size of dataset after duplicates removed:")
print(dim(A1a)[1])

# --- Adjusting data with matrices A1a ---
common_ids <- Reduce(intersect, list(
  A1a$ISO,
  rownames(phylomatrix),
  rownames(spmatrix)
))
A1a <- A1a[A1a$ISO %in% common_ids, ]
A1a_phylomatrix <- phylomatrix[common_ids, common_ids]
A1a_spmatrix <- spmatrix[common_ids, common_ids]

## --- A1b --- ####

# --- Merging Datasets A1b ---
A1b <- merge(
  nee22[, 
        c("ISO", "L1_pop", "region")],
  nee24[, c("ISO693.3", "Phoneme.Inventory.Size", "Island.Endemic")], by.x = "ISO", by.y = "ISO693.3", all.x = TRUE)

# --- Data Cleanup A1b ---
# Remove all NA entries
print("Size of initial dataset:")
print(dim(A1b)[1])
old_A1b <- A1b
A1b <- na.omit(A1b) 
print("Size of dataset with NA removed:")
print(dim(A1b)[1])
# Remove iso duplicate entries
print("Duplicate iso entries:")
iso_duplicates <- which(table(A1b$ISO) != 1)
print(which(table(A1b$ISO) != 1))
A1b <- A1b[!(A1b$ISO %in% names(which(table(A1b$ISO) != 1))),]
print("Size of dataset after duplicates removed:")
print(dim(A1b)[1])

# --- Adjusting data with matrices A1b ---
common_ids <- Reduce(intersect, list(
  A1b$ISO,
  rownames(phylomatrix),
  rownames(spmatrix)
))
A1b <- A1b[A1b$ISO %in% common_ids, ]
A1b_phylomatrix <- phylomatrix[common_ids, common_ids]
A1b_spmatrix <- spmatrix[common_ids, common_ids]

## --- A2 --- ####

## Merging Datasets
A2 <- merge(
  nee22[, 
        c("ISO", "region", "L1_pop", "bordering_language_richness", "altitude_range", "area")],
  nee24[, c("ISO693.3", "Phoneme.Inventory.Size", "Island.Endemic")], 
  by.x = "ISO", by.y = "ISO693.3", all.x = TRUE)

# --- Data Cleanup ---
# Remove all NA entries
paste("Size of initial dataset:", dim(A2)[1])
old_A2 <- A2
A2 <- na.omit(A2) 
paste("Size of dataset with NA removed:", dim(A2)[1])
# Remove iso duplicate entries
print("Duplicate iso entries:")
iso_duplicates <- duplicated(A2$ISO)
print(A2$ISO[iso_duplicates])
A2 <- A2[!unlist(iso_duplicates),]
paste("Size of dataset after duplicates removed:", dim(A2)[1])

# --- Adjusting data with matrices ---
common_ids <- Reduce(intersect, list(
  A2$ISO,
  rownames(phylomatrix),
  rownames(spmatrix)
))
A2 <- A2[A2$ISO %in% common_ids, ]
A2_phylomatrix <- phylomatrix[common_ids, common_ids]
A2_spmatrix <- spmatrix[common_ids, common_ids]

## --- A3a --- ####

# --- Prepare dataset ---
# Merge datasets
A3a <- merge(nee22[, c("ISO", "region", "bordering_language_richness")], 
             nee24[, c("ISO693.3", "L1_pop","Phoneme.Inventory.Size", "Island.Endemic", 
                       "Distance.to.Mainland", "Distance.to.Continent", "Range.Size..km2.")], 
             by.x = "ISO", by.y = "ISO693.3", all.x = TRUE)
# Remove all NA entries
A3a <- na.omit(A3a) 
paste("Size of dataset with NA remove:", dim(A3a)[1])
paste("#Island.Endemic: ", sum(A3a$Island.Endemic))
paste("#Island.Endemic to Total Ratio: ", sum(A3a$Island.Endemic) / dim(A3a)[1])
paste("Dataset with all languages",dim(A3a)[[1]])

# --- Adjusting A3a data with matrices ---
common_ids <- Reduce(intersect, list(
  A3a$ISO,
  rownames(phylomatrix),
  rownames(spmatrix)
))
A3a <- A3a[A3a$ISO %in% common_ids, ]
A3a_phylomatrix <- phylomatrix[common_ids, common_ids]
A3a_spmatrix <- spmatrix[common_ids, common_ids]

## --- A3b --- ####

# --- Prepare dataset ---
# Merge datasets
A3b <- merge(nee22[, c("ISO", "region")], 
             nee24[, c("ISO693.3", "L1_pop","Phoneme.Inventory.Size", "Island.Endemic", 
                       "Distance.to.Mainland", "Distance.to.Continent", "Range.Size..km2.")], 
             by.x = "ISO", by.y = "ISO693.3", all.x = TRUE)
# Remove all NA entries
A3b <- na.omit(A3b) 
paste("Size of dataset with NA remove:", dim(A3b)[1])
paste("#Island.Endemic: ", sum(A3b$Island.Endemic))
paste("#Island.Endemic to Total Ratio: ", sum(A3b$Island.Endemic) / dim(A3b)[1])
# Keep only island endemic languages
paste("Dataset with all languages",dim(A3b)[[1]])
A3b <- A3b[A3b$Island.Endemic == 1,]
paste("Dataset with only island endemic languages",dim(A3b)[[1]])

# --- Adjusting A3b data with matrices ---
common_ids <- Reduce(intersect, list(
  A3b$ISO,
  rownames(phylomatrix),
  rownames(spmatrix)
))
A3b <- A3b[A3b$ISO %in% common_ids, ]
A3b_phylomatrix <- phylomatrix[common_ids, common_ids]
A3b_spmatrix <- spmatrix[common_ids, common_ids]

## --- A4 --- ####

# --- Prepare dataset ---
# Merge datasets
A4 <- merge(nee22[, c("ISO", "region", "documentation")], 
            nee24[, c("ISO693.3", "Phoneme.Inventory.Size")], 
            by.x = "ISO", by.y = "ISO693.3", all.x = TRUE)
# Remove all NA entries
A4 <- na.omit(A4) 
paste("Size of dataset with NA remove:", dim(A4)[1])

# --- Adjusting A4 data with matrices ---
common_ids <- Reduce(intersect, list(
  A4$ISO,
  rownames(phylomatrix),
  rownames(spmatrix)
))
A4 <- A4[A4$ISO %in% common_ids, ]
A4_phylomatrix <- phylomatrix[common_ids, common_ids]
A4_spmatrix <- spmatrix[common_ids, common_ids]

# --- (4) EXPORT DATASETS FOR ANALYSIS --- ####

## --- A1a --- ####

# --- Preview and save data A1a ---
paste("Size of dataset adjusted:", dim(A1a)[1])
paste("#Island.Endemic: ", sum(A1a$Island.Endemic))
paste("#Island.Endemic to Total Ratio: ", sum(A1a$Island.Endemic) / dim(A1a)[1])
write.csv(A1a, file = "output/A1a/A1a_adjusted.csv", row.names=FALSE)

# --- Preview and save matrices A1a ---
print(dim(A1a_phylomatrix))
print(dim(A1a_spmatrix))
print(A1a_phylomatrix[1:10,1:10])
print(A1a_spmatrix[1:10,1:10])
save(A1a_phylomatrix, file = "output/A1a/A1a_phylomatrix.RData") 
save(A1a_spmatrix, file = "output/A1a/A1a_spmatrix.RData")

## --- A1b --- ####

# --- Preview and save data A1b ---
paste("Size of dataset adjusted:", dim(A1b)[1])
paste("#Island.Endemic: ", sum(A1b$Island.Endemic))
paste("#Island.Endemic to Total Ratio: ", sum(A1b$Island.Endemic) / dim(A1b)[1])
write.csv(A1b, file = "output/A1b/A1b_adjusted.csv", row.names=FALSE)

# --- Preview and save matrices A1b ---
print(dim(A1b_phylomatrix))
print(dim(A1b_spmatrix))
print(A1b_phylomatrix[1:10,1:10])
print(A1b_spmatrix[1:10,1:10])
save(A1b_phylomatrix, file = "output/A1b/A1b_phylomatrix.RData") 
save(A1b_spmatrix, file = "output/A1b/A1b_spmatrix.RData")

## --- A2 --- ####

# --- Preview and save data ---
paste("Size of dataset adjusted:", dim(A2)[1])
paste("#Island.Endemic: ", sum(A2$Island.Endemic))
paste("#Island.Endemic to Total Ratio: ", sum(A2$Island.Endemic) / dim(A2)[1])
write.csv(A2, file = "output/A2/A2_adjusted.csv", row.names=FALSE)

### Preview and save matrices
print(dim(A2_phylomatrix))
print(dim(A2_spmatrix))
print(A2_phylomatrix[1:10,1:10])
print(A2_spmatrix[1:10,1:10])
save(A2_phylomatrix, file = "output/A2/A2_phylomatrix.RData") 
save(A2_spmatrix, file = "output/A2/A2_spmatrix.RData")

## --- A3a --- ####

# Preview and save data
paste("Size of dataset adjusted:", dim(A3a)[1])
paste("#Island.Endemic: ", sum(A3a$Island.Endemic))
paste("#Island.Endemic to Total Ratio: ", sum(A3a$Island.Endemic) / dim(A3a)[1])
# save data 
write.csv(A3a, file = "output/A3a/A3a_adjusted.csv", row.names=FALSE)

# Preview and save matrices
print(dim(A3a_phylomatrix))
print(dim(A3a_spmatrix))
print(A3a_phylomatrix[1:10,1:10])
print(A3a_spmatrix[1:10,1:10])
save(A3a_phylomatrix, file = "output/A3a/A3a_phylomatrix.RData") 
save(A3a_spmatrix, file = "output/A3a/A3a_spmatrix.RData")

## --- A3b --- ####

# Preview and save data
paste("Size of dataset adjusted:", dim(A3b)[1])
paste("#Island.Endemic: ", sum(A3b$Island.Endemic))
paste("#Island.Endemic to Total Ratio: ", sum(A3b$Island.Endemic) / dim(A3b)[1])
# save data 
write.csv(A3b, file = "output/A3b/A3b_adjusted.csv", row.names=FALSE)

# Preview and save matrices
print(dim(A3b_phylomatrix))
print(dim(A3b_spmatrix))
print(A3b_phylomatrix[1:10,1:10])
print(A3b_spmatrix[1:10,1:10])
save(A3b_phylomatrix, file = "output/A3b/A3b_phylomatrix.RData") 
save(A3b_spmatrix, file = "output/A3b/A3b_spmatrix.RData")

## --- A4 --- ####

# Preview and save data
paste("Size of dataset adjusted:", dim(A4)[1])
# save data 
write.csv(A4, file = "output/A4/A4_adjusted.csv", row.names=FALSE)

# Preview and save matrices
print(dim(A4_phylomatrix))
print(dim(A4_spmatrix))
print(A4_phylomatrix[1:10,1:10])
print(A4_spmatrix[1:10,1:10])
save(A4_phylomatrix, file = "output/A4/A4_phylomatrix.RData") 
save(A4_spmatrix, file = "output/A4/A4_spmatrix.RData")

## CLEAR ALL CURRENT VARIABLES ####
rm(list = ls())

# --- (5) DATASET VISUALISATION --- ####

## Setup ####

region_order <- c("Oceania", "Australia and New Zealand", 
                  "South-Eastern Asia", "Southern Asia", "Asia", 
                  "Europe", "North Africa and Arabia", "Africa", 
                  "Western Africa", "Northern America", "Central America", "South America")

# Retrieve nee22 (if not already loaded)
if(!exists("nee22")){
  nee22 <- read.csv("data/Global predictors of language endangerment and the future of linguistic diversity Data 2.csv",
                    stringsAsFactors = FALSE) 
}

## --- A1a --- ####

# --- Get raw A1a adjusted data and matrices ---
raw_A1a <- read.csv("output/A1a/A1a_adjusted.csv")
load("output/A1a/A1a_phylomatrix.RData")
load("output/A1a/A1a_spmatrix.RData")

# --- Preview raw data A1a ---
head(raw_A1a)
summary(raw_A1a)
# Phoneme Inventory Size distribution
PIS_hist <- ggplot(raw_A1a, aes(x = Phoneme.Inventory.Size)) +
  geom_histogram() +
  labs(
    y = "Number of Languages Observed",
    x = "PISa",
  ) + 
  scale_x_continuous(n.breaks = 15) +
  scale_y_continuous(n.breaks = 15)

ggsave(
  filename = "output/A1a/PIS_hist.png",
  plot = PIS_hist,
  scale = 1
)
print(PIS_hist)

# --- Log transform variables A1a ---
A1a <- raw_A1a
A1a$Phoneme.Inventory.Size <- log(A1a$Phoneme.Inventory.Size) 
A1a$L1_pop <- log(A1a$L1_pop + 0.5)

# --- Mapping phoneme inventory sizes onto world map A1a ---
# Get world map data
world <- map_data("world")
PIS_world_map <-ggplot() +
  # Base world map
  geom_polygon(
    data = world,
    aes(x = long, y = lat, group = group),
    fill = "gray90", color = "gray50"
  ) +
  # Points with size and color based on count
  geom_point(
    data = A1a,
    aes(x = Longitude, y = Latitude, color = log(Phoneme.Inventory.Size)),
    alpha = 0.7
  ) +
  scale_color_viridis_c(option="plasma") +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(x = "", y = "", color = "Log PISa")
ggsave("output/A1a/PIS_world_map.png", plot = PIS_world_map, width = 15, height = 10, dpi = 300)
PIS_world_map

# --- Preview correlations A1a ---
# Plot data
phoneme_L1_scatter = ggplot(A1a, aes(x = L1_pop, y = Phoneme.Inventory.Size, color = factor(region, levels = region_order))) + 
  geom_point() + 
  # stat_ellipse(geom = "polygon", alpha = 0.2, aes(fill = group)) + 
  labs(x = "L1_pop", y = "PISa", color = "Region") + # title = "Phoneme Count ~ L1 Population Size Scatterplot",  
  theme_minimal()
# theme(legend.position = "none") + # Hide the legend for now, often too many groups
# Save graph
ggsave("output/A1a/phoneme_to_L1_scatter.png", plot = phoneme_L1_scatter, width = 15, height = 10, dpi = 100,scale = 0.5)
# Show graph
print(phoneme_L1_scatter)

# --- Testing for region sampling bias A1a ---
# number of repetitions
n_rep = 1e4
# Get region names
region_names <- region_order
# Find number of languages in each region for A1a
A1a_region_tally <- as.data.frame(table(A1a$region))
A1a_region_tally <- A1a_region_tally[match(region_names, A1a_region_tally$Var1), ]
region_sample_tallies <- data.frame(matrix(ncol = length(region_names), nrow = n_rep))
colnames(region_sample_tallies) <- region_names # Set column names
# Repeated n_rep times, randomly sample without replacement 
for(i in 1:n_rep){
  # Get region sample
  region_sample <- sample(nee22$region, size = length(A1a$region), replace = FALSE)
  # Current region tally.
  region_sample_tally = table(region_sample)
  # Add current tally to all tally data
  region_sample_tallies[i, names(region_sample_tally)] <- as.numeric(region_sample_tally)
}
# Reshape the dataframe to long format
region_sample_tallies_long <- pivot_longer(region_sample_tallies, everything(), names_to = c("Group"), values_to = 'Value')
# Set 'Group' as a factor with levels in the desired order, this is to ensure that ggplot doesn't alphabetically sort the regions.
region_sample_tallies_long$Group <- factor(region_sample_tallies_long$Group, levels = region_names)
region_box <- ggplot(region_sample_tallies_long, aes(x = Group, y = Value)) +
  geom_boxplot(varwidth = TRUE) +
  geom_point(data = A1a_region_tally, 
             aes(x = region_names, y = Freq), 
             color = "red",
             size = 2) + 
  ylab("Number of Languages Sampled") +
  xlab("Region") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.title.x.bottom = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y.left = element_text(size = 10))
# Save graph
ggsave("output/A1a/region_sampling.png", plot = region_box, width = 15, height = 15, dpi = 300, scale = 0.5)
print(region_box)

# --- Phoneme.Inventory.Size across regions A1a ---
box_phoneme_regions <- ggplot(raw_A1a, aes(
  x = factor(region, levels = region_order),
  y = Phoneme.Inventory.Size
)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = median(raw_A1a$Phoneme.Inventory.Size), color = "red") +
  labs(
    x = "Region",
    y = "PISa"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(box_phoneme_regions)
ggsave("output/A1a/phoneme_region_box.png", plot = box_phoneme_regions, width = 15, height = 15, dpi = 300, scale = 0.5)

## --- A1b --- ####

# --- Preview raw data A1b ---

# --- Get raw A1b adjusted data and matrices ---
raw_A1b <- read.csv("output/A1b/A1b_adjusted.csv")
load("output/A1b/A1b_phylomatrix.RData")
load("output/A1b/A1b_spmatrix.RData")

head(raw_A1b)
summary(raw_A1b)
# Phoneme Inventory Size distribution
PIS_hist <- ggplot(raw_A1b, aes(x = Phoneme.Inventory.Size)) +
  geom_histogram() +
  labs(
    y = "Number of Languages Observed",
    x = "PISb",
  ) + 
  scale_x_continuous(n.breaks = 15) +
  scale_y_continuous(n.breaks = 15)

ggsave(
  filename = "output/A1b/PIS_hist.png",
  plot = PIS_hist,
  scale = 1
)
print(PIS_hist)

# --- Log transform variables A1b ---
A1b<- raw_A1b
A1b$Phoneme.Inventory.Size <- log(A1b$Phoneme.Inventory.Size) 
A1b$L1_pop <- log(A1b$L1_pop + 0.5)

# --- Preview correlations A1b ---
# Plot data
phoneme_L1_scatter = ggplot(A1b, aes(x = L1_pop, y = Phoneme.Inventory.Size, color = factor(region, levels = region_order))) + 
  geom_point() + 
  # stat_ellipse(geom = "polygon", alpha = 0.2, aes(fill = group)) + 
  labs(x = "L1_pop", y = "PISb", color = "Region") + # title = "Phoneme Count ~ L1 Population Size Scatterplot",  
  theme_minimal()
# theme(legend.position = "none") + # Hide the legend for now, often too many groups
# Save graph
ggsave("output/A1b/phoneme_to_L1_scatter.png", plot = phoneme_L1_scatter, width = 15, height = 10, dpi = 100,scale = 0.5)
# Show graph
print(phoneme_L1_scatter)

# --- Testing for region sampling bias A1b ---
# Retrieve nee22 (if not already loaded)
if(!exists("nee22")){
  nee22 <- read.csv("data/Global predictors of language endangerment and the future of linguistic diversity Data 2.csv",
                    stringsAsFactors = FALSE) 
}
# number of repetitions
n_rep = 1e4
# Get region names
region_names <- region_order
# Find number of languages in each region for A1b
A1b_region_tally <- as.data.frame(table(A1b$region))
A1b_region_tally <- A1b_region_tally[match(region_names, A1b_region_tally$Var1), ]
region_sample_tallies <- data.frame(matrix(ncol = length(region_names), nrow = n_rep))
colnames(region_sample_tallies) <- region_names # Set column names
# Repeated n_rep times, randomly sample without replacement 
for(i in 1:n_rep){
  # Get region sample
  region_sample <- sample(nee22$region, size = length(A1b$region), replace = FALSE)
  # Current region tally.
  region_sample_tally = table(region_sample)
  # Add current tally to all tally data
  region_sample_tallies[i, names(region_sample_tally)] <- as.numeric(region_sample_tally)
}
# Reshape the dataframe to long format
region_sample_tallies_long <- pivot_longer(region_sample_tallies, everything(), names_to = c("Group"), values_to = 'Value')
# Set 'Group' as a factor with levels in the desired order, this is to ensure that ggplot doesn't alphabetically sort the regions.
region_sample_tallies_long$Group <- factor(region_sample_tallies_long$Group, levels = region_names)
region_box <- ggplot(region_sample_tallies_long, aes(x = Group, y = Value)) +
  geom_boxplot(varwidth = TRUE) +
  geom_point(data = A1b_region_tally, 
             aes(x = region_names, y = Freq), 
             color = "red",
             size = 2) + 
  ylab("Number of Languages Sampled") +
  xlab("Region") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.title.x.bottom = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y.left = element_text(size = 10))
# Save graph
ggsave("output/A1b/region_sampling.png", plot = region_box, width = 15, height = 15, dpi = 300, scale = 0.5)
print(region_box)

# --- Phoneme.Inventory.Size across regions A1b ---
box_phoneme_regions <- ggplot(raw_A1b, aes(
  x = factor(region, levels = region_order),
  y = Phoneme.Inventory.Size
)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = median(raw_A1b$Phoneme.Inventory.Size), color = "red") +
  labs(
    x = "Region",
    y = "PISb"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(box_phoneme_regions)
ggsave("output/A1b/phoneme_region_box.png", plot = box_phoneme_regions, width = 15, height = 15, dpi = 300, scale = 0.5)

## --- A2 --- ####

# Get raw A2 adjusted data and matrices
raw_A2 <- read.csv("output/A2/A2_adjusted.csv")
load("output/A2/A2_phylomatrix.RData")
load("output/A2/A2_spmatrix.RData")

# --- Log transform variables ---
A2 <- raw_A2
A2$Phoneme.Inventory.Size <- log(A2$Phoneme.Inventory.Size) 
A2$area <- log(A2$area)
A2$altitude_range <- log(A2$altitude_range + 0.5)
A2$L1_pop <- log(A2$L1_pop + 0.5)
A2$bordering_language_richness <- log(A2$bordering_language_richness + 0.5)

# --- Preview transformed data ---
head(A2)
summary(A2)
PIS_Island <- ggplot(data = A2, aes(x = as.factor(Island.Endemic), y = Phoneme.Inventory.Size)) +
  geom_boxplot()
print(PIS_Island)

# --- Preview correlations ---
pairs(A2[,
         c("Phoneme.Inventory.Size", "area", "L1_pop", 
           "altitude_range", "bordering_language_richness")],
      panel = function(x, y, ...) {
        points(x, y, cex = 0.1, ...) # cex = 0.5 makes points half the default size
      })
cor(A2[,
       c("Phoneme.Inventory.Size", "area", "L1_pop", 
         "altitude_range", "bordering_language_richness")])
cor.test(A2$Phoneme.Inventory.Size,A2$Island.Endemic)

## --- A3a --- ####

# Get raw A3a adjusted data and matrices
raw_A3a <- read.csv("output/A3a/A3a_adjusted.csv")
load("output/A3a/A3a_phylomatrix.RData")
load("output/A3a/A3a_spmatrix.RData")

# --- Log transform variables ---
A3a <- raw_A3a
A3a$Phoneme.Inventory.Size <- log(A3a$Phoneme.Inventory.Size) 
A3a$Range.Size..km2. <- log(A3a$Range.Size..km2.)
A3a$L1_pop <- log(A3a$L1_pop + 0.5)
A3a$bordering_language_richness <- log(A3a$bordering_language_richness + 0.5)
# Do not transform zero entries for distance variables
A3a$Distance.to.Mainland[A3a$Distance.to.Mainland != 0] <- log(
  A3a$Distance.to.Mainland[A3a$Distance.to.Mainland != 0])
A3a$Distance.to.Continent[A3a$Distance.to.Continent != 0] <- log(
  A3a$Distance.to.Continent[A3a$Distance.to.Continent != 0])

# --- Preview transformed data ---
head(A3a)
summary(A3a)
# PIS ~ Range_Size, Grouping by Island Endemics
PIS_Range <- ggplot(data = A3a,aes(x = Range.Size..km2., y = Phoneme.Inventory.Size)) +
  geom_point()
print(PIS_Range)
ggsave(
  filename = "output/A3a/PIS_Range.png",
  plot = PIS_Range,
  scale = 1,
  width=7,
  height=5
)

# --- Preview correlations ---
pairs(A3a[,
          c("Range.Size..km2.",
            "Distance.to.Mainland", "Distance.to.Continent",
            "L1_pop", "bordering_language_richness")],
      panel = function(x, y, ...) {
        points(x, y, cex = 0.1, ...) # cex = 0.5 makes points half the default size
      })
cor(A3a[,
        c("Range.Size..km2.",
          "Distance.to.Mainland", "Distance.to.Continent",
          "L1_pop", "bordering_language_richness")])


## --- A3b --- ####

# Get raw A3b adjusted data and matrices
raw_A3b <- read.csv("output/A3b/A3b_adjusted.csv")
load("output/A3b/A3b_phylomatrix.RData")
load("output/A3b/A3b_spmatrix.RData")

# --- Log transform variables ---
A3b <- raw_A3b
A3b$Phoneme.Inventory.Size <- log(A3b$Phoneme.Inventory.Size) 
A3b$Range.Size..km2. <- log(A3b$Range.Size..km2.)
A3b$L1_pop <- log(A3b$L1_pop + 0.5)
# Do not transform zero entries for distance variables
A3b$Distance.to.Mainland[A3b$Distance.to.Mainland != 0] <- log(
  A3b$Distance.to.Mainland[A3b$Distance.to.Mainland != 0])
A3b$Distance.to.Continent[A3b$Distance.to.Continent != 0] <- log(
  A3b$Distance.to.Continent[A3b$Distance.to.Continent != 0])

# --- Preview transformed data ---
head(A3b)
summary(A3b)
# PIS ~ Range_Size, Grouping by Island Endemics
PIS_Range <- ggplot(data = A3b,aes(x = Range.Size..km2., y = Phoneme.Inventory.Size)) +
  geom_point()
print(PIS_Range)
ggsave(
  filename = "output/A3b/PIS_Range.png",
  plot = PIS_Range,
  scale = 1,
  width=7,
  height=5
)

# --- Preview correlations ---
pairs(A3b[,
          c("Range.Size..km2.",
            "Distance.to.Mainland", "Distance.to.Continent",
            "L1_pop")],
      panel = function(x, y, ...) {
        points(x, y, cex = 0.1, ...) # cex = 0.5 makes points half the default size
      })
cor(A3b[,
        c("Range.Size..km2.",
          "Distance.to.Mainland", "Distance.to.Continent",
          "L1_pop")])

## --- A4 --- ####

# Get raw A4 adjusted data and matrices
raw_A4 <- read.csv("output/A4/A4_adjusted.csv")
load("output/A4/A4_phylomatrix.RData")
load("output/A4/A4_spmatrix.RData")

# --- Log transform variables ---
A4 <- raw_A4
A4$Phoneme.Inventory.Size <- log(A4$Phoneme.Inventory.Size) 

# --- Preview transformed data ---
head(A4)
summary(A4)
print("Number of languages for each documentation level:")
print(table(A4$documentation))

# --- Preview correlations ---
# PIS ~ Documentation
PIS_Doc <- ggplot(data = A4, aes(x = as.factor(documentation), y = Phoneme.Inventory.Size)) +
  geom_boxplot()
print(PIS_Doc)

# --- (6) OLS AND GLS ANALYSIS --- ####

## Setup ####
A1_model <- c(Phoneme.Inventory.Size ~ L1_pop)

## --- A1a --- ####

# --- Ordinary Least Squares A1a ---
A1a_l1 <- lm(formula = Phoneme.Inventory.Size ~ L1_pop, data = A1a)
print(summary(A1a_l1))
print(BIC(A1a_l1))
print(logLik(A1a_l1))

# --- Best p values for all models A1a ---
A1a_p.res <- sbplx(c(0.5, 0.5, 0.5),
                   best_p,
                   formula=A1_model[[1]],
                   data=A1a,
                   spmatrix=A1a_spmatrix,phylomatrix=A1a_phylomatrix,
                   lower=c(0,0,0),upper=c(1,1,1),
                   nl.info = TRUE)
# Save p.res
save(A1a_p.res, file = "output/A1a/A1a_p_res.RData")

# --- Maximum likelihood fits for all models A1a ---
A1a_model <- ml_fit(p=A1a_p.res$par,formula=A1_model[[1]],data=A1a,spmatrix=A1a_spmatrix,phylomatrix=A1a_phylomatrix)
# Save model fit
save(A1a_model, file = "output/A1a/A1a_model.RData") 

## --- A1b --- ####

# --- Ordinary Least Squares A1b ---
A1b_l1 <- lm(formula = Phoneme.Inventory.Size ~ L1_pop, data = A1b)
print(summary(A1b_l1))
print(BIC(A1b_l1))
print(logLik(A1b_l1))

# --- Best p values for all models A1b ---
A1b_p.res <- sbplx(c(0.5, 0.5, 0.5),
                   best_p,
                   formula=A1_model[[1]],
                   data=A1b,
                   spmatrix=A1b_spmatrix,phylomatrix=A1b_phylomatrix,
                   lower=c(0,0,0),upper=c(1,1,1),
                   nl.info = TRUE)
# Save p.res
save(A1b_p.res, file = "output/A1b/A1b_p_res.RData")

# --- Maximum likelihood fits for all models A1b ---
A1b_model <- ml_fit(p=A1b_p.res$par,formula=A1_model[[1]],data=A1b,spmatrix=A1b_spmatrix,phylomatrix=A1b_phylomatrix)
# Save model fit
save(A1b_model, file = "output/A1b/A1b_model.RData") 

## --- A2 --- ####

# --- Ordinary Least Squares ---
A2_l1 <- lm(formula = Phoneme.Inventory.Size ~ L1_pop, data = A2)
print(summary(A2_l1))
print(BIC(A2_l1))
print(logLik(A2_l1))

# --- Model predictor combinations ---
A2_predictors <- c("Island.Endemic", "area",
                "altitude_range", "bordering_language_richness",
                "L1_pop")
response <- "Phoneme.Inventory.Size"
n_pred <- length(A2_predictors)
pred_combs <- sapply(1:n_pred, function(x) combn(A2_predictors, x))
A2_models <- c(Phoneme.Inventory.Size ~ 1)
for (i in 1:n_pred) {
  for (j in 1:dim(pred_combs[[i]])[2]) {
    combination <- pred_combs[[i]][,j]
    # Accounting for interactions
    # if(all(c("Island.Endemic","area") %in% combination)){
    #   int_model <- (paste(response, paste(combination, "Island.Endemic:area", collapse="+"), sep="~"))
    #   int_model <- as.formula(int_model)
    #   int_models <- c(models,int_model)
    # }
    # if(all(c("L1_pop","area") %in% combination)){
    #   int_model <- (paste(response, paste(combination, "L1_pop:area", collapse="+"), sep="~"))
    #   int_model <- as.formula(int_model)
    #   int_models <- c(models,int_model)
    # }
    # if(all(c("bordering_language_richness","area") %in% combination)){
    #   int_model <- (paste(response, paste(combination, "bordering_language_richness:area", collapse="+"), sep="~"))
    #   int_model <- as.formula(int_model)
    #   int_models <- c(models,int_model)
    # }
    model <- (paste(response, paste(combination, collapse="+"), sep="~"))
    model <- as.formula(model)
    A2_models <- c(A2_models,model)
  }
}
unlist(A2_models)

# --- Best p values for all models ---
A2_p.res <- sbplx(c(0.5, 0.5, 0.5),
               best_p,
               formula=Phoneme.Inventory.Size~1,
               data=A2,
               spmatrix=A2_spmatrix,phylomatrix=A2_phylomatrix,
               lower=c(0,0,0),upper=c(1,1,1),
               nl.info = TRUE)
# Save p.res and covariance matrix
save(A2_p.res, file = "output/A2/A2_p_res.RData")
spmatrix_temp <- A2_spmatrix/max(A2_spmatrix)
spmatrix_temp <- exp(-(spmatrix_temp/A2_p.res$par[2])^2)
A2_mat <- as.matrix((A2_p.res$par[3]*(1-A2_p.res$par[1])*spmatrix_temp+(1-A2_p.res$par[1])*(1-A2_p.res$par[3])*A2_phylomatrix+A2_p.res$par[1]*diag(dim(A2_phylomatrix)[1])))
save(A2_mat, file = "output/A2/A2_mat.RData")

# --- Maximum likelihood fits for all models ---
A2_fits <- mclapply(
  X = A2_models,
  FUN = function(f) ml_fit(p=A2_p.res$par,formula=f,data=A2,spmatrix=A2_spmatrix,phylomatrix=A2_phylomatrix),
  mc.cores = 8
)
# Save models
saveRDS(A2_fits, "output/A2/A2_fits.RDS") 

# --- Model Comparison ---
# Make empty NA list of coefficients and their respective p values in tuple form
n <- length(A2_fits)
A2_summaries <- list("Island.Endemic" = as.list(rep(NA,n)),
                     "area" = as.list(rep(NA,n)),
                     "altitude_range" = as.list(rep(NA,n)), 
                     "bordering_language_richness" = as.list(rep(NA,n)),
                     "L1_pop" = as.list(rep(NA,n)))
ml_logLik <- c()
ml_BIC <- c()
ml_MSE <- c()

for (i in 1:n){
  # get current model summary
  ml_summary <- summary(A2_fits[[i]])
  # if intercept only model just get BIC, logLik and MLE
  if(dim(ml_summary$tTable)[1] == 1) {
    ml_logLik <- c(ml_logLik, ml_summary$logLik)
    ml_BIC <- c(ml_BIC, ml_summary$BIC)
    ml_MSE <- c(ml_MSE, sum(A2_fits[[i]]$residuals**2)/dim(A2)[1])
    next
  }
  # Get coefficient values
  coef_val <- ml_summary$tTable[,1]
  # Get corresponding p-values
  coef_p <- ml_summary$tTable[,4]
  # Get coefficient names
  coef_names <- names(ml_summary$tTable[,4])
  # Get number of coefficients
  n_coef <- length(coef_names)
  # Add coef info into model_summaries
  for (j in 2:n_coef){
    if(coef_p[j]<=0.05){
      A2_summaries[[coef_names[j]]][[i]] <- paste0(round(coef_val[j], 5),"***")
    } else {
      A2_summaries[[coef_names[j]]][[i]] <- paste0(round(coef_val[j], 5))
    }
  }
  # Add BIC, logLik and MLE
  ml_logLik <- c(ml_logLik, ml_summary$logLik)
  ml_BIC <- c(ml_BIC, ml_summary$BIC)
  ml_MSE <- c(ml_MSE, sum(A2_fits[[i]]$residuals**2)/dim(A2)[1])
}

# Create new data frame from lists
A2_ml_data <- data.frame(BIC = ml_BIC,
                      logLik = ml_logLik,
                      MSE = ml_MSE,
                      unlist(A2_summaries$Island.Endemic),
                      unlist(A2_summaries$area),
                      unlist(A2_summaries$altitude_range),
                      unlist(A2_summaries$bordering_language_richness),
                      unlist(A2_summaries$L1_pop),
                      stringsAsFactors = FALSE)
colnames(A2_ml_data) <- c("BIC", "logLik", "MSE", A2_predictors)
# Sort based on increasing BIC
A2_ml_data <- A2_ml_data[order(A2_ml_data$BIC), ]
# Add delta BIC column (min_BIC - current_BIC)
min_BIC <-  min(A2_ml_data$BIC)
A2_ml_data$deltaBIC <- A2_ml_data$BIC - min_BIC
# Save output
write.csv(A2_ml_data, file = "output/A2/A2_ml_data.csv", row.names=FALSE)

## --- A3a --- ####

# --- Model predictor combinations ---
A3a_predictors <- c("Range.Size..km2.",
                "Distance.to.Mainland", "Distance.to.Continent",
                "L1_pop", "bordering_language_richness")
response <- "Phoneme.Inventory.Size"
n_pred <- length(A3a_predictors)
pred_combs <- sapply(1:n_pred, function(x) combn(A3a_predictors, x))
A3a_models <- c(Phoneme.Inventory.Size~1)
for (i in 1:n_pred) {
  #print(pred_combs[[i]])
  for (j in 1:dim(pred_combs[[i]])[2]) {
    model <- (paste(response, paste(pred_combs[[i]][,j], collapse="+"), sep="~"))
    model <- as.formula(model)
    A3a_models <- c(A3a_models,model)
  }
}

# --- Best p values for all models ---
A3a_p.res <- sbplx(c(0.5, 0.5, 0.5),
               best_p,
               formula=Phoneme.Inventory.Size~1,
               data=A3a,
               spmatrix=A3a_spmatrix,phylomatrix=A3a_phylomatrix,
               lower=c(0,0,0),upper=c(1,1,1),
               nl.info = TRUE)
# Save p.res and covariance matrix
save(A3a_p.res, file = "output/A3a/A3a_p_res.RData")
spmatrix_temp <- A3a_spmatrix/max(A3a_spmatrix)
spmatrix_temp <- exp(-(spmatrix_temp/A3a_p.res$par[2])^2)
A3a_mat <- as.matrix((A3a_p.res$par[3]*(1-A3a_p.res$par[1])*spmatrix_temp+(1-A3a_p.res$par[1])*(1-A3a_p.res$par[3])*A3a_phylomatrix+A3a_p.res$par[1]*diag(dim(A3a_phylomatrix)[1])))
save(A3a_mat, file = "output/A3a/A3a_mat.RData")

# --- Maximum likelihood fits for all models ---
A3a_fits <- mclapply(
  X = A3a_models,
  FUN = function(f) ml_fit(p=A3a_p.res$par,formula=f,data=A3a,spmatrix=A3a_spmatrix,phylomatrix=A3a_phylomatrix),
  mc.cores = 8
)
# Save models
saveRDS(A3a_fits, "output/A3a/A3a_fits.RDS") 

# --- Model Comparison ---
# Make empty NA list of coefficients and their respective p values in tuple form
n <- length(A3a_fits)
A3a_summaries <- list("bordering_language_richness" = as.list(rep(NA,n)),
                     "Range.Size..km2." = as.list(rep(NA,n)),
                     "Distance.to.Mainland" = as.list(rep(NA,n)), 
                     "Distance.to.Continent" = as.list(rep(NA,n)),
                     "L1_pop" = as.list(rep(NA,n)))
ml_logLik <- c()
ml_BIC <- c()
ml_MSE <- c()

for (i in 1:n){
  # get current model summary
  ml_summary <- summary(A3a_fits[[i]])
  # if intercept only model just get BIC, logLik and MLE
  if(dim(ml_summary$tTable)[1] == 1) {
    ml_logLik <- c(ml_logLik, ml_summary$logLik)
    ml_BIC <- c(ml_BIC, ml_summary$BIC)
    ml_MSE <- c(ml_MSE, sum(A3a_fits[[i]]$residuals**2)/dim(A3a)[1])
    next
  }
  # Get coefficient values
  coef_val <- ml_summary$tTable[,1]
  # Get corresponding p-values
  coef_p <- ml_summary$tTable[,4]
  # Get coefficient names
  coef_names <- names(ml_summary$tTable[,4])
  # Get number of coefficients
  n_coef <- length(coef_names)
  # Add coef info into model_summaries
  for (j in 2:n_coef){
    if(coef_p[j]<=0.05){
      A3a_summaries[[coef_names[j]]][[i]] <- paste0(round(coef_val[j], 5),"***")
    } else {
      A3a_summaries[[coef_names[j]]][[i]] <- paste0(round(coef_val[j], 5))
    }
  }
  # Add BIC, logLik and MLE
  ml_logLik <- c(ml_logLik, ml_summary$logLik)
  ml_BIC <- c(ml_BIC, ml_summary$BIC)
  ml_MSE <- c(ml_MSE, sum(A3a_fits[[i]]$residuals**2)/dim(A3a)[1])
}

# Create new data frame from lists
A3a_ml_data <- data.frame(BIC = ml_BIC,
                      logLik = ml_logLik,
                      MSE = ml_MSE,
                      unlist(A3a_summaries$bordering_language_richness),
                      unlist(A3a_summaries$Range.Size..km2.),
                      unlist(A3a_summaries$Distance.to.Mainland),
                      unlist(A3a_summaries$Distance.to.Continent),
                      unlist(A3a_summaries$L1_pop),
                      stringsAsFactors = FALSE)
colnames(A3a_ml_data) <- c("BIC", "logLik", "MSE", A3a_predictors)
# Sort based on increasing BIC
A3a_ml_data <- A3a_ml_data[order(A3a_ml_data$BIC), ]
# Add delta BIC column (min_BIC - current_BIC)
min_BIC <-  min(A3a_ml_data$BIC)
A3a_ml_data$deltaBIC <- A3a_ml_data$BIC - min_BIC
# Save output
write.csv(A3a_ml_data, file = "output/A3a/A3a_ml_data.csv", row.names=FALSE)

## --- A3b --- ####

# --- Model predictor combinations ---
A3b_predictors <- c("Range.Size..km2.",
                "Distance.to.Mainland", "Distance.to.Continent",
                "L1_pop")
response <- "Phoneme.Inventory.Size"
n_pred <- length(A3b_predictors)
pred_combs <- sapply(1:n_pred, function(x) combn(A3b_predictors, x))
A3b_models <- c(Phoneme.Inventory.Size~1)
for (i in 1:n_pred) {
  #print(pred_combs[[i]])
  for (j in 1:dim(pred_combs[[i]])[2]) {
    model <- (paste(response, paste(pred_combs[[i]][,j], collapse="+"), sep="~"))
    model <- as.formula(model)
    A3b_models <- c(A3b_models,model)
  }
}

# --- Best p values for all models ---
A3b_p.res <- sbplx(c(0.5, 0.5, 0.5),
               best_p,
               formula=Phoneme.Inventory.Size~1,
               data=A3b,
               spmatrix=A3b_spmatrix,phylomatrix=A3b_phylomatrix,
               lower=c(0,0,0),upper=c(1,1,1),
               nl.info = TRUE)
# Save p.res and covariance matrix
save(A3b_p.res, file = "output/A3b/A3b_p_res.RData")
spmatrix_temp <- A3b_spmatrix/max(A3b_spmatrix)
spmatrix_temp <- exp(-(spmatrix_temp/A3b_p.res$par[2])^2)
A3b_mat <- as.matrix((A3b_p.res$par[3]*(1-A3b_p.res$par[1])*spmatrix_temp+(1-A3b_p.res$par[1])*(1-A3b_p.res$par[3])*A3b_phylomatrix+A3b_p.res$par[1]*diag(dim(A3b_phylomatrix)[1])))
save(A3b_mat, file = "output/A3b/A3b_mat.RData")

# --- Maximum likelihood fits for all models ---
A3b_fits <- mclapply(
  X = A3b_models,
  FUN = function(f) ml_fit(p=A3b_p.res$par,formula=f,data=A3b,spmatrix=A3b_spmatrix,phylomatrix=A3b_phylomatrix),
  mc.cores = 8
)
# Save models
saveRDS(A3b_fits, "output/A3b/A3b_fits.RDS") 

# --- Model Comparison ---
# Make empty NA list of coefficients and their respective p values in tuple form
n <- length(A3b_fits)
A3b_summaries <- list("Range.Size..km2." = as.list(rep(NA,n)),
                     "Distance.to.Mainland" = as.list(rep(NA,n)), 
                     "Distance.to.Continent" = as.list(rep(NA,n)),
                     "L1_pop" = as.list(rep(NA,n)))
ml_logLik <- c()
ml_BIC <- c()
ml_MSE <- c()

for (i in 1:n){
  # get current model summary
  ml_summary <- summary(A3b_fits[[i]])
  # if intercept only model just get BIC, logLik and MLE
  if(dim(ml_summary$tTable)[1] == 1) {
    ml_logLik <- c(ml_logLik, ml_summary$logLik)
    ml_BIC <- c(ml_BIC, ml_summary$BIC)
    ml_MSE <- c(ml_MSE, sum(A3b_fits[[i]]$residuals**2)/dim(A3b)[1])
    next
  }
  # Get coefficient values
  coef_val <- ml_summary$tTable[,1]
  # Get corresponding p-values
  coef_p <- ml_summary$tTable[,4]
  # Get coefficient names
  coef_names <- names(ml_summary$tTable[,4])
  # Get number of coefficients
  n_coef <- length(coef_names)
  # Add coef info into model_summaries
  for (j in 2:n_coef){
    if(coef_p[j]<=0.05){
      A3b_summaries[[coef_names[j]]][[i]] <- paste0(round(coef_val[j], 5),"***")
    } else {
      A3b_summaries[[coef_names[j]]][[i]] <- paste0(round(coef_val[j], 5))
    }
  }
  # Add BIC, logLik and MLE
  ml_logLik <- c(ml_logLik, ml_summary$logLik)
  ml_BIC <- c(ml_BIC, ml_summary$BIC)
  ml_MSE <- c(ml_MSE, sum(A3b_fits[[i]]$residuals**2)/dim(A3b)[1])
}

# Create new data frame from lists
A3b_ml_data <- data.frame(BIC = ml_BIC,
                      logLik = ml_logLik,
                      MSE = ml_MSE,
                      unlist(A3b_summaries$Range.Size..km2.),
                      unlist(A3b_summaries$Distance.to.Mainland),
                      unlist(A3b_summaries$Distance.to.Continent),
                      unlist(A3b_summaries$L1_pop),
                      stringsAsFactors = FALSE)
colnames(A3b_ml_data) <- c("BIC", "logLik", "MSE", A3b_predictors)
# Sort based on increasing BIC
A3b_ml_data <- A3b_ml_data[order(A3b_ml_data$BIC), ]
# Add delta BIC column (min_BIC - current_BIC)
min_BIC <-  min(A3b_ml_data$BIC)
A3b_ml_data$deltaBIC <- A3b_ml_data$BIC - min_BIC
# Save output
write.csv(A3b_ml_data, file = "output/A3b/A3b_ml_data.csv", row.names=FALSE)

## --- A4 --- ####

# --- Model predictor combinations ---
A4_model <- c(Phoneme.Inventory.Size ~ documentation)

# --- Best p values for all models ---
A4_p.res <- sbplx(c(0.5, 0.5, 0.5),
               best_p,
               formula=A4_model[[1]],
               data=A4,
               spmatrix=A4_spmatrix,phylomatrix=A4_phylomatrix,
               lower=c(0,0,0),upper=c(1,1,1),
               nl.info = TRUE)
# Save p.res and covariance matrix
save(A4_p.res, file = "output/A4/A4_p_res.RData")
spmatrix_temp <- A4_spmatrix/max(A4_spmatrix)
spmatrix_temp <- exp(-(spmatrix_temp/A4_p.res$par[2])^2)
mat <- as.matrix((A4_p.res$par[3]*(1-A4_p.res$par[1])*spmatrix_temp+(1-A4_p.res$par[1])*(1-A4_p.res$par[3])*A4_phylomatrix+A4_p.res$par[1]*diag(dim(A4_phylomatrix)[1])))

# --- Maximum likelihood fits for all models A1b ---
A4_model <- ml_fit(p=A4_p.res$par,formula=A4_model[[1]],data=A4,spmatrix=A4_spmatrix,phylomatrix=A4_phylomatrix)
# Save model fit
save(A4_model, file = "output/A4/A4_model.RData") 

# --- (7) REFERENCES --- ####
#
# (1) Bromham, L., Yaxley, K.J. & Cardillo, M. Islands are engines of language diversity. Nat Ecol Evol 8, 1991–2002 (2024). https://doi.org/10.1038/s41559-024-02488-4
# 
# (2) Hua, X., Greenhill, S.J., Cardillo, M. et al. The ecological drivers of variation in global language diversity. Nat Commun 10, 2047 (2019). https://doi.org/10.1038/s41467-019-09842-2
# 
# (3) Bromham, L., Dinnage, R., Skirgård, H. et al. Global predictors of language endangerment and the future of linguistic diversity. Nat Ecol Evol 6, 163–173 (2022). https://doi.org/10.1038/s41559-021-01604-y
# 
# (4) Cormac Anderson, Tiago Tresoldi, Simon J Greenhill, Robert Forkel, Russell Gray, Johann-Mattis List, Variation in phoneme inventories: quantifying the problem and improving comparability, Journal of Language Evolution, Volume 8, Issue 2, July 2023, Pages 149–168, https://doi.org/10.1093/jole/lzad011
# 