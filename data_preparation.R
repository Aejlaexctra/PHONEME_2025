# --- Reference Guide ---

# (1) Global predictors of language endangerment and the future of linguistic diversity
# (2) Variation in phoneme inventories: quantifying the problem and improving comparability
# (3) Islands are engines of language diversity

# --- Helper Functions

# Reshaping phylogenetic and spatial matrices to match data
reshape_matrix <- function(matrix, col){
  col <- intersect(col, rownames(matrix))  # Ensure all elements are valid
  matrix <- matrix[col, col, drop = FALSE] # Reorder both rows and cols
  return(matrix)
}

# --- Retrieve all data ---

# Fetch (1) data
global_pred_data <- read.csv("Data/Global predictors of language endangerment and the future of linguistic diversity Data 2.csv")
# Fetch (1) contact matrices
# Fetch (1) distance/contact matrices
load("Data/distance_matrices/Wnb.Rdata")
load("Data/distance_matrices/Wphy.Rdata")
load("Data/distance_matrices/Wsp.Rdata")
conmatrix <- Wnb
phylomatrix <- Wphy
spmatrix <- Wsp
# Fetch (2) data
var_phon_inv_spec_data <- read.csv("Data/cldf-datasets-inventory-study/PHOIBLE-data.csv")
var_phon_inv_glot_data <- read.csv("Data/cldf-datasets-inventory-study/phoible/cldf/languages.csv")
# Fetch (3) data
island_engine_data <- read.csv("Data/Islands are engines of language diversity data.csv")
# # Fetch (3) phylogenetic and spatial distance matrices
# phylomatrix <- read.csv("Data/distance_matrices/phylogenetic_covariance_matrix.csv")
# rownames(phylomatrix) <- phylomatrix[,1]
# phylomatrix <- phylomatrix[,-1]
# phylomatrix <- as.matrix(phylomatrix)
# spmatrix <- read.csv("Data/distance_matrices/spatial_distance_matrix.csv")
# rownames(spmatrix) <- spmatrix[,1]
# spmatrix <- spmatrix[,-1]
# spmatrix <- as.matrix(spmatrix)

# --- Merge predictor variables from datasets ---

# (1) and (2) For glottocodes, isocodes, L1 pop, roughness, altitudinal range
# and language documentation
data <- merge(
  global_pred_data[, c("ISO", "L1_pop", "region", "bordering_language_richness", "roughness", "altitude_range", "area", "documentation")],
  var_phon_inv_glot_data[, c("ISO639P3code", "Glottocode")], by.x = "ISO", by.y = "ISO639P3code", all.x = TRUE)
# (2) For latitude, longitude
# data <- merge(data, var_phon_inv_spec_data[, c("Glottocode", "Latitude", "Longitude")],
#               by.x = "Glottocode", by.y = "Glottocode", all.x = TRUE)
# 

### ERROR ABOVE HERE, var_phon_inv_spec_data

# (2) For Phoneme Counts
data <- merge(data, var_phon_inv_spec_data[, c("Glottocode", "Sounds")], 
                      by.x = "Glottocode", by.y = "Glottocode", all.x = TRUE)
# (3) For Island Endemic, Distance to Mainland, Distance to Continent
data <- merge(data, island_engine_data[, c("ISO693.3", "Island.Endemic", "Distance.to.Mainland", "Distance.to.Continent")], 
              by.x = "ISO", by.y = "ISO693.3", all.x = TRUE)

# --- Clean up and format dataset ---

# Convert Island Endemic variable from boolean to binary
data$Island.Endemic <- as.integer(as.logical(data$Island.Endemic))
# Convert Documention variable from string to integer
data$documentation[data$documentation == "little or none"] <- "0"
data$documentation[data$documentation == "basic"] <- "1"
data$documentation[data$documentation == "detailed"] <- "2"
data$documentation <- as.integer(data$documentation)
# Change Arab region name to Middle East
# data$region[data$region == "Arab"] <- "Middle East"

# Change Sounds column name to phoneme_count
names(data)[names(data) == "Sounds"] <- "phoneme_count"
# Change Island.Endemic to island_endemic
names(data)[names(data) == "Island.Endemic"] <- "island_endemic"
# Change Distance.to.Continent to continent_dist
names(data)[names(data) == "Distance.to.Continent"] <- "continent_dist"
# Change Distance.to.Mainland to mainland_dist
names(data)[names(data) == "Distance.to.Mainland"] <- "mainland_dist"

#Remove all NA value containing rows
paste("Size of initial dataset:", length(data$Glottocode))
paste("#Island_Endemic: ", dim(data[data$island_endemic==1,][1])[1])
paste("#Island_Endemic to Total Ratio: ", dim(data[data$island_endemic==1,])[1] / dim(data[data$island_endemic==0,])[1])
old_data <- data
data <- na.omit(data) 
paste("Number of incomplete (NA) language entries removed:", 
      length(old_data$Glottocode) - length(data$Glottocode))

# Checking for duplicates, as some ISO codes can represent multiple unique glottocodes
dup_group = data$ISO
duplicates <- data[duplicated(dup_group) | duplicated(dup_group, fromLast = TRUE), ]
print(duplicates$Glottocode)
paste("Number of duplicate languages removed:", length(duplicates$Glottocode))
# Remove duplicate rows based on the 'ISO' column
data <- data[!duplicated(data$ISO), ] 
# Sort entries based on alphabetical order in ISO column
data <- data[order(data$ISO), ] 
paste("#Island_Endemic: ", dim(data[data$island_endemic==1,][1])[1])
paste("#Island_Endemic to Total Ratio: ", dim(data[data$island_endemic==1,])[1] / dim(data[data$island_endemic==0,])[1])

# --- Adjust Size of Matrices ---
# spmatrix <- reshape_matrix(spmatrix, data$ISO)
# phylomatrix <- reshape_matrix(phylomatrix, data$ISO)
# conmatrix <- reshape_matrix(conmatrix, data$ISO)
# index <- data$ISO %in% rownames(phylomatrix)
# data <- data[index,] 
common_ids <- Reduce(intersect, list(
  data$ISO,
  rownames(phylomatrix),
  rownames(conmatrix),
  rownames(spmatrix)
))
data <- data[data$ISO %in% common_ids, ]
phylomatrix <- phylomatrix[common_ids, common_ids]
conmatrix   <- conmatrix[common_ids, common_ids]
spmatrix    <- spmatrix[common_ids, common_ids]

paste("Size of final dataset:", length(data$ISO))

# --- Normalise Rows Again ---
conmatrix <- t(apply(conmatrix, 1, function(x) {
  s <- sum(x)
  if (s > 0) x / s else x  # keep row of zeros if sum is 0
}))
spmatrix <- t(apply(spmatrix, 1, function(x) {
  s <- sum(x)
  if (s > 0) x / s else x  # keep row of zeros if sum is 0
}))
phylomatrix <- t(apply(phylomatrix, 1, function(x) {
  s <- sum(x)
  if (s > 0) x / s else x  # keep row of zeros if sum is 0
}))

# --- Export Data ---

write.csv(data, "data/raw_data.csv", row.names=FALSE)
save(spmatrix, file="data/spmatrix.Rdata")
save(phylomatrix, file="data/phylomatrix.Rdata")
save(conmatrix, file="data/conmatrix.Rdata")
