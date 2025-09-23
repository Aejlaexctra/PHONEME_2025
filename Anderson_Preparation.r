# --- Setting random seed ---
set.seed(42)

# --- Loading Phylogenetic, Spatial Distance and Contact Matrices ---
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
load("data/distance_matrices/old_Wnb.Rdata")

# Restoring contact matrix to a binary matrix
Wnb[Wnb != 0] <- 1 

### Preview Matrices
print(phylomatrix[1:5,1:5])
print(spmatrix[1:5,1:5])

# --- Loading NEE24 (1) ---
nee24 <- read.csv("data/Islands are engines of language diversity data.csv",
                  stringsAsFactors = FALSE)

# --- NEE24 Cleanup ---
# Convert variable types
nee24$Island.Endemic <- as.numeric(nee24$Island.Endemic)

# --- Preview data ---
head(nee24)
summary(nee24)

# --- Loading NEE22 (3) ---
nee22 <- read.csv("data/Global predictors of language endangerment and the future of linguistic diversity Data 2.csv",
                  stringsAsFactors = FALSE)

# --- NEE22 Cleanup ---
# Convert Documention variable from string to integer
nee22$documentation[nee22$documentation == "little or none"] <- "0"
nee22$documentation[nee22$documentation == "basic"] <- "1"
nee22$documentation[nee22$documentation == "detailed"] <- "2"
nee22$documentation <- as.integer(nee22$documentation)

# --- Preview data ---
head(nee22)
summary(nee22)

# --- Loading phoible (4) ---
var_phon_inv_spec_data <- read.csv("Data/cldf-datasets-inventory-study/PHOIBLE-data.csv")
var_phon_inv_glot_data <- read.csv("Data/cldf-datasets-inventory-study/phoible/cldf/languages.csv")
phoible <- merge(
  var_phon_inv_spec_data[, c("Glottocode", "Sounds", "Latitude", "Longitude")],
  var_phon_inv_glot_data[, c("ISO639P3code", "Glottocode")], by.x = "Glottocode", by.y = "Glottocode", all.x = TRUE)

# --- Phoible Cleanup ---
# Remove duplicate glottocode entries
print("Duplicate Glottocode entries:")
glotto_duplicates <- duplicated(phoible$Glottocode)
print(phoible$Glottocode[glotto_duplicates])
phoible <- phoible[!unlist(glotto_duplicates),]
paste("Size of dataset after duplicates removed:", dim(phoible)[1])

# rename Sounds column to Phoneme.Inventory.Size
colnames(phoible)[colnames(phoible) == "Sounds"] <- "Phoneme.Inventory.Size"

# --- Preview data ---
head(phoible)
summary(phoible)

## Merging Datasets
data <- merge(
  nee22[, 
        c("ISO", "L1_pop", "region", "bordering_language_richness", 
          "roughness", "altitude_range", "documentation")],
  phoible[, c("ISO639P3code", "Glottocode", "Latitude", "Longitude", "Phoneme.Inventory.Size")], by.x = "ISO", by.y = "ISO639P3code", all.x = TRUE)
data <- merge(data, nee24[, c("ISO693.3", "Island.Endemic", "Distance.to.Mainland", "Distance.to.Continent", "Range.Size..km2.")], 
              by.x = "ISO", by.y = "ISO693.3", all.x = TRUE)

# --- Data Cleanup ---
# Remove all NA entries
paste("Size of initial dataset:", dim(data)[1])
old_data <- data
data <- na.omit(data) 
paste("Size of dataset with NA removed:", dim(data)[1])
# Remove iso duplicate entries
print("Duplicate iso entries:")
iso_duplicates <- duplicated(data$ISO)
print(data$ISO[iso_duplicates])
data <- data[!unlist(iso_duplicates),]
paste("Size of dataset after duplicates removed:", dim(data[1]))

# --- Adjusting data with matrices ---
common_ids <- Reduce(intersect, list(
  data$ISO,
  rownames(phylomatrix),
  rownames(spmatrix)
))
data <- data[data$ISO %in% common_ids, ]
phylomatrix <- phylomatrix[common_ids, common_ids]
spmatrix <- spmatrix[common_ids, common_ids]

# --- Preview and save data ---
paste("Size of dataset adjusted:", dim(data)[1])
paste("#Island_Endemic: ", sum(data$Island.Endemic))
paste("#Island_Endemic to Total Ratio: ", sum(data$Island.Endemic) / dim(data)[1])
write.csv(data, file = "data/anderson_adjusted.csv", row.names=FALSE)

### Preview and save matrices
print(dim(phylomatrix))
print(dim(spmatrix))
print(phylomatrix[1:10,1:10])
print(spmatrix[1:10,1:10])
save(phylomatrix, file = "data/anderson/phylomatrix.RData") 
save(spmatrix, file = "data/anderson/spmatrix.RData")
save(Wnb, file = "data/anderson/Wnb.Rdata")