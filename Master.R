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
# Phylogenetic distance matrix
phylomatrix <- read.csv("data/distance_matrices/phylogenetic_covariance_matrix.csv",
                        check.names = FALSE)
rownames(phylomatrix) <- phylomatrix[,1]
phylomatrix <- phylomatrix[,-1]
phylomatrix <- as.matrix(phylomatrix)
# Spatial distance matrix
spmatrix <- read.csv("data/distance_matrices/spatial_distance_matrix.csv",
                     check.names = FALSE)
rownames(spmatrix) <- spmatrix[,1]
spmatrix <- spmatrix[,-1]
spmatrix <- as.matrix(spmatrix)

## --- Loading NEE22 (3) --- ####
nee22 <- read.csv("data/Global predictors of language endangerment and the future of linguistic diversity Data 2.csv",
                  stringsAsFactors = FALSE)
colnames(nee22)[colnames(nee22) == "island"] <- "Island.Endemic" # Rename Island predictor to Island.Endemic
nee22$documentation[nee22$documentation == "little or none"] <- "0" # Convert Documentation predictor from string to integer
nee22$documentation[nee22$documentation == "basic"] <- "1"
nee22$documentation[nee22$documentation == "detailed"] <- "2"
nee22$documentation <- as.integer(nee22$documentation)
nee22$region[nee22$region == "Arab"] <- "North Africa and Arabia" # Change Region "Arab" to "North Africa and Arabia"

## --- Loading NEE24 (1) --- ####
nee24 <- read.csv("data/Islands are engines of language diversity data.csv",
                  stringsAsFactors = FALSE)
# Rename L1.Population predictor to L1_pop
colnames(nee24)[colnames(nee24) == "L1.Population"] <- "L1_pop"
# Set Island.Endemic predictor as numeric
nee24$Island.Endemic <- as.numeric(nee24$Island.Endemic)

## --- Loading Phoible (4) --- ####
var_phon_inv_spec_data <- read.csv("Data/cldf-datasets-inventory-study/PHOIBLE-data.csv")
var_phon_inv_glot_data <- read.csv("Data/cldf-datasets-inventory-study/phoible/cldf/languages.csv")
phoible <- merge(
  var_phon_inv_spec_data[, c("Glottocode", "Sounds", "Latitude", "Longitude")],
  var_phon_inv_glot_data[, c("ISO639P3code", "Glottocode")], by.x = "Glottocode", by.y = "Glottocode", all.x = TRUE)
colnames(phoible)[colnames(phoible) == "Sounds"] <- "Phoneme.Inventory.Size"
glotto_duplicates <- duplicated(phoible$Glottocode) # Remove duplicate glottocode entries
phoible <- phoible[!unlist(glotto_duplicates),]

# --- (2) FUNCTION DEFINITIONS --- ####

## --- GLS (2) Setup --- ####
best_p <- function (p,formula,data,spmatrix,phylomatrix,opt="SP") {
  spmatrix <- spmatrix/max(spmatrix)
  spmatrix <- exp(-(spmatrix/p[2])^2)
  if(opt=="SP"){ # Accounting for Spatial and Phylogenetic Autocorrelation
    mat <- as.matrix((p[3]*(1-p[1])*spmatrix+(1-p[1])*(1-p[3])*phylomatrix+p[1]*diag(dim(phylomatrix)[1]))) 
  } else if(opt == "S"){ # Accounting for Spatial Autocorrelation
    mat <- as.matrix((1-p[1])*spmatrix+p[1]*diag(dim(spmatrix)[1]))
  } else { # Accounting for Phylogenetic Autocorrelation
    mat <- as.matrix((1-p[1])*phylomatrix+p[1]*diag(dim(phylomatrix)[1]))
  }
  res <- try(gls(model=formula,data=data,correlation=corSymm(mat[lower.tri(mat)],fixed=T),method="ML"),silent=T)
  if (inherits(res,"try-error")) {
    out <- -10000
  } else {
    out <- res$logLik
    # if (res$logLik>0) {out <- -10000} # Positive log likelihoods tolerated, unlike function in (2)
  }
  -out
}

ml_fit <- function (p,formula,data,spmatrix,phylomatrix,opt="SP") {
  spmatrix <- spmatrix/max(spmatrix)
  spmatrix <- exp(-(spmatrix/p[2])^2)
  if(opt=="SP"){ # Accounting for Spatial and Phylogenetic Autocorrelation
    mat <- as.matrix((p[3]*(1-p[1])*spmatrix+(1-p[1])*(1-p[3])*phylomatrix+p[1]*diag(dim(phylomatrix)[1]))) 
  } else if(opt == "S"){ # Accounting for Spatial Autocorrelation
    mat <- as.matrix((1-p[1])*spmatrix+p[1]*diag(dim(spmatrix)[1]))
  } else { # Accounting for Phylogenetic Autocorrelation
    mat <- as.matrix((1-p[1])*phylomatrix+p[1]*diag(dim(phylomatrix)[1]))
  }
  res <- try(gls(model=formula,data=data,correlation=corSymm(mat[lower.tri(mat)],fixed=T),method="ML"),silent=T)
  res
}

# --- (3) CLEANING AND FORMATING DATASETS --- ####

## --- A1a --- ####

# --- Merging Datasets A1a ---
# Dataset for Section 1, with variables:
# PISa, L1_pop
A1a <- merge(
  nee22[, 
        c("ISO", "L1_pop", "region", "Island.Endemic")],
  phoible[, c("ISO639P3code", "Glottocode", "Latitude", "Longitude", "Phoneme.Inventory.Size")], by.x = "ISO", by.y = "ISO639P3code", all.x = TRUE)
A1a <- na.omit(A1a) # Remove all NA entries
iso_duplicates <- which(table(A1a$ISO) != 1) # Remove iso duplicate entries
A1a <- A1a[!(A1a$ISO %in% names(which(table(A1a$ISO) != 1))),]

# --- Adjusting data with matrices A1a ---
common_ids <- Reduce(intersect, list(
  A1a$ISO,
  rownames(phylomatrix),
  rownames(spmatrix)
)) # Get common ISO codes across A1a, phylogenetic and spatial distance matrices
A1a <- A1a[A1a$ISO %in% common_ids, ] # Remove languages with ISO codes not common to all
A1a_phylomatrix <- phylomatrix[common_ids, common_ids]
A1a_spmatrix <- spmatrix[common_ids, common_ids]

## --- A1b --- ####

# --- Merging Datasets A1b ---
# Dataset for Section 1 (Supplementary), with variables:
# PISc, L1_pop (Island predictor not in OLS/GLS analysis but for plots)
A1b <- merge(
  nee22[, 
        c("ISO", "region")],
  nee24[, c("ISO693.3", "Phoneme.Inventory.Size", "L1_pop", "Island.Endemic")], by.x = "ISO", by.y = "ISO693.3", all.x = TRUE)
A1b <- na.omit(A1b) # Remove all NA entries
iso_duplicates <- which(table(A1b$ISO) != 1) # Remove iso duplicate entries
A1b <- A1b[!(A1b$ISO %in% names(which(table(A1b$ISO) != 1))),]

# --- Adjusting data with matrices A1b ---
common_ids <- Reduce(intersect, list(
  A1b$ISO,
  rownames(phylomatrix),
  rownames(spmatrix)
)) # Get common ISO codes across A1b, phylogenetic and spatial distance matrices
A1b <- A1b[A1b$ISO %in% common_ids, ] # Remove languages with ISO codes not common to all
A1b_phylomatrix <- phylomatrix[common_ids, common_ids]
A1b_spmatrix <- spmatrix[common_ids, common_ids]

## --- A2 --- ####

## Merging Datasets
# Dataset for Section 2, with variables:
# PISc, Bordering, Altitude, L1_pop, Island and Area 
A2 <- merge(
  nee22[, 
        c("ISO", "region", "bordering_language_richness", "altitude_range")],
  nee24[, c("ISO693.3", "Phoneme.Inventory.Size", "L1_pop", "Island.Endemic", "Range.Size..km2.")], 
  by.x = "ISO", by.y = "ISO693.3", all.x = TRUE)
A2 <- na.omit(A2) # Remove all NA entries
iso_duplicates <- duplicated(A2$ISO) # Remove iso duplicate entries
A2 <- A2[!unlist(iso_duplicates),]

# --- Adjusting data with matrices ---
common_ids <- Reduce(intersect, list(
  A2$ISO,
  rownames(phylomatrix),
  rownames(spmatrix)
)) # Get common ISO codes across A2, phylogenetic and spatial distance matrices
A2 <- A2[A2$ISO %in% common_ids, ] # Remove languages with ISO codes not common to all
A2_phylomatrix <- phylomatrix[common_ids, common_ids]
A2_spmatrix <- spmatrix[common_ids, common_ids]

## --- A3a --- ####

# --- Prepare dataset ---
## Merging datasets
# Dataset for Section 3, with variables:
# PISc, Bordering, L1_pop, Island, Mainland, Continent andArea
A3a <- merge(nee22[, c("ISO", "region", "bordering_language_richness")], 
             nee24[, c("ISO693.3", "L1_pop","Phoneme.Inventory.Size", "Island.Endemic", 
                       "Distance.to.Mainland", "Distance.to.Continent", "Range.Size..km2.")], 
             by.x = "ISO", by.y = "ISO693.3", all.x = TRUE)
A3a <- na.omit(A3a) # Remove all NA entries

# --- Adjusting A3a data with matrices ---
common_ids <- Reduce(intersect, list(
  A3a$ISO,
  rownames(phylomatrix),
  rownames(spmatrix)
)) # Get common ISO codes across A3a, phylogenetic and spatial distance matrices
A3a <- A3a[A3a$ISO %in% common_ids, ] # Remove languages with ISO codes not common to all
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

## --- A5a --- ####
A5_shared_ISO <- intersect(A1a$ISO,A2$ISO) # Get shared languages between A1a (PISa) and A2 (PISc) datasets

# --- Prepare dataset ---
## Trimming datasets of A1a
# Dataset for Re-analysis of languages in both PISa and PISc , with variables:
A5a <- A1a[A1a$ISO %in% A5_shared_ISO,]

# --- Adjusting A5a data with matrices ---
common_ids <- Reduce(intersect, list(
  A5a$ISO,
  rownames(phylomatrix),
  rownames(spmatrix)
)) # Get common ISO codes across A5a, phylogenetic and spatial distance matrices
A5a <- A5a[A5a$ISO %in% common_ids, ] # Remove languages with ISO codes not common to all
A5a_phylomatrix <- phylomatrix[common_ids, common_ids]
A5a_spmatrix <- spmatrix[common_ids, common_ids]

## --- A5b --- ####

# --- Prepare dataset ---
## Trimming datasets of A2
# Dataset for Re-analysis of languages in both PISa and PISc , with variables:
A5b <- A2[A2$ISO %in% A5_shared_ISO,]

# --- Adjusting A5b data with matrices ---
common_ids <- Reduce(intersect, list(
  A5b$ISO,
  rownames(phylomatrix),
  rownames(spmatrix)
)) # Get common ISO codes across A5b, phylogenetic and spatial distance matrices
A5b <- A5b[A5b$ISO %in% common_ids, ] # Remove languages with ISO codes not common to all
A5b_phylomatrix <- phylomatrix[common_ids, common_ids]
A5b_spmatrix <- spmatrix[common_ids, common_ids]

## --- A5c --- ####

# --- Prepare dataset ---
## Trimming datasets of A3a
# Dataset for Re-analysis of languages in both PISa and PISc , with variables:
A5c <- A3a[A3a$ISO %in% A5_shared_ISO,]

# --- Adjusting A5c data with matrices ---
common_ids <- Reduce(intersect, list(
  A5c$ISO,
  rownames(phylomatrix),
  rownames(spmatrix)
)) # Get common ISO codes across A5c, phylogenetic and spatial distance matrices
A5c <- A5c[A5c$ISO %in% common_ids, ] # Remove languages with ISO codes not common to all
A5c_phylomatrix <- phylomatrix[common_ids, common_ids]
A5c_spmatrix <- spmatrix[common_ids, common_ids]

## --- A5d --- ####

# --- Prepare dataset ---
## Trimming datasets of A4
# Dataset for Re-analysis of languages in both PISa and PISc , with variables:
A5d <- A4[A4$ISO %in% A5_shared_ISO,]

# --- Adjusting A5d data with matrices ---
common_ids <- Reduce(intersect, list(
  A5d$ISO,
  rownames(phylomatrix),
  rownames(spmatrix)
)) # Get common ISO codes across A5d, phylogenetic and spatial distance matrices
A5d <- A5d[A5d$ISO %in% common_ids, ] # Remove languages with ISO codes not common to all
A5d_phylomatrix <- phylomatrix[common_ids, common_ids]
A5d_spmatrix <- spmatrix[common_ids, common_ids]

## --- A6 --- ####

# --- Prepare dataset ---
# Merge datasets
A6 <- merge(nee22[, c("ISO", "region", "documentation")], 
            nee24[, c("ISO693.3", "Phoneme.Inventory.Size","L1_pop")], 
            by.x = "ISO", by.y = "ISO693.3", all.x = TRUE)
# Remove all NA entries
A6 <- na.omit(A6) 
paste("Size of dataset with NA remove:", dim(A6)[1])

# --- Adjusting A4 data with matrices ---
common_ids <- Reduce(intersect, list(
  A4$ISO,
  rownames(phylomatrix),
  rownames(spmatrix)
))
A6 <- A6[A6$ISO %in% common_ids, ]
A6_phylomatrix <- phylomatrix[common_ids, common_ids]
A6_spmatrix <- spmatrix[common_ids, common_ids]

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

## --- A5a --- ####

# Preview and save data
paste("Size of dataset adjusted:", dim(A5a)[1])
# save data 
write.csv(A5a, file = "output/A5a/A5a_adjusted.csv", row.names=FALSE)

# Preview and save matrices
print(dim(A5a_phylomatrix))
print(dim(A5a_spmatrix))
print(A5a_phylomatrix[1:10,1:10])
print(A5a_spmatrix[1:10,1:10])
save(A5a_phylomatrix, file = "output/A5a/A5a_phylomatrix.RData") 
save(A5a_spmatrix, file = "output/A5a/A5a_spmatrix.RData")

## --- A5b --- ####

# Preview and save data
paste("Size of dataset adjusted:", dim(A5b)[1])
# save data 
write.csv(A5b, file = "output/A5b/A5b_adjusted.csv", row.names=FALSE)

# Preview and save matrices
print(dim(A5b_phylomatrix))
print(dim(A5b_spmatrix))
print(A5b_phylomatrix[1:10,1:10])
print(A5b_spmatrix[1:10,1:10])
save(A5b_phylomatrix, file = "output/A5b/A5b_phylomatrix.RData") 
save(A5b_spmatrix, file = "output/A5b/A5b_spmatrix.RData")

## --- A5c --- ####

# Preview and save data
paste("Size of dataset adjusted:", dim(A5c)[1])
# save data 
write.csv(A5c, file = "output/A5c/A5c_adjusted.csv", row.names=FALSE)

# Preview and save matrices
print(dim(A5c_phylomatrix))
print(dim(A5c_spmatrix))
print(A5c_phylomatrix[1:10,1:10])
print(A5c_spmatrix[1:10,1:10])
save(A5c_phylomatrix, file = "output/A5c/A5c_phylomatrix.RData") 
save(A5c_spmatrix, file = "output/A5c/A5c_spmatrix.RData")

## --- A5d --- ####

# Preview and save data
paste("Size of dataset adjusted:", dim(A5d)[1])
# save data 
write.csv(A5d, file = "output/A5d/A5d_adjusted.csv", row.names=FALSE)

# Preview and save matrices
print(dim(A5d_phylomatrix))
print(dim(A5d_spmatrix))
print(A5d_phylomatrix[1:10,1:10])
print(A5d_spmatrix[1:10,1:10])
save(A5d_phylomatrix, file = "output/A5d/A5d_phylomatrix.RData") 
save(A5d_spmatrix, file = "output/A5d/A5d_spmatrix.RData")

## --- A6 --- ####

# Preview and save data
paste("Size of dataset adjusted:", dim(A6)[1])
# save data 
write.csv(A6, file = "output/A6/A6_adjusted.csv", row.names=FALSE)

# Preview and save matrices
print(dim(A6_phylomatrix))
print(dim(A6_spmatrix))
print(A6_phylomatrix[1:10,1:10])
print(A6_spmatrix[1:10,1:10])
save(A6_phylomatrix, file = "output/A6/A6_phylomatrix.RData") 
save(A6_spmatrix, file = "output/A6/A6_spmatrix.RData")

# --- (5) TRANSFORM VARIABLES --- ####

## --- A1a --- ####
A1a$Phoneme.Inventory.Size <- log(A1a$Phoneme.Inventory.Size) 
A1a$L1_pop <- log(A1a$L1_pop + 0.5)

## --- A1b ---  ####
A1b$Phoneme.Inventory.Size <- log(A1b$Phoneme.Inventory.Size) 
A1b$L1_pop <- log(A1b$L1_pop + 0.5)

## --- A2 ---  ####

# --- Log transform variables ---
A2$Phoneme.Inventory.Size <- log(A2$Phoneme.Inventory.Size) 
A2$Range.Size..km2. <- log(A2$Range.Size..km2.)
A2$altitude_range <- log(A2$altitude_range + 0.5)
A2$L1_pop <- log(A2$L1_pop + 0.5)
A2$bordering_language_richness <- log(A2$bordering_language_richness + 0.5)

# --- Standardise predictors ---
A2$Range.Size..km2. <- scale(A2$Range.Size..km2.)
A2$altitude_range <- scale(A2$altitude_range)
A2$L1_pop <- scale(A2$L1_pop)
A2$bordering_language_richness <- scale(A2$bordering_language_richness)

## --- A3a ---  ####

# --- Log transform variables ---
A3a$Phoneme.Inventory.Size <- log(A3a$Phoneme.Inventory.Size) 
A3a$Range.Size..km2. <- log(A3a$Range.Size..km2.)
A3a$L1_pop <- log(A3a$L1_pop + 0.5)
A3a$bordering_language_richness <- log(A3a$bordering_language_richness + 0.5)
# Do not transform zero entries for distance variables
A3a$Distance.to.Mainland[A3a$Distance.to.Mainland != 0] <- log(
  A3a$Distance.to.Mainland[A3a$Distance.to.Mainland != 0])
A3a$Distance.to.Continent[A3a$Distance.to.Continent != 0] <- log(
  A3a$Distance.to.Continent[A3a$Distance.to.Continent != 0])

# --- Standardise predictors ---
A3a$Range.Size..km2. <- scale(A3a$Range.Size..km2.)
A3a$L1_pop <- scale(A3a$L1_pop)
A3a$bordering_language_richness <- scale(A3a$bordering_language_richness)
A3a$Distance.to.Mainland <- scale(A3a$Distance.to.Mainland)
A3a$Distance.to.Continent <- scale(A3a$Distance.to.Continent)

## --- A3b ---  ####

# --- Log transform variables ---
A3b$Phoneme.Inventory.Size <- log(A3b$Phoneme.Inventory.Size) 
A3b$Range.Size..km2. <- log(A3b$Range.Size..km2.)
A3b$L1_pop <- log(A3b$L1_pop + 0.5)
# Do not transform zero entries for distance variables
A3b$Distance.to.Mainland[A3b$Distance.to.Mainland != 0] <- log(
  A3b$Distance.to.Mainland[A3b$Distance.to.Mainland != 0])
A3b$Distance.to.Continent[A3b$Distance.to.Continent != 0] <- log(
  A3b$Distance.to.Continent[A3b$Distance.to.Continent != 0])

# --- Standardise predictors ---
A3b$Range.Size..km2. <- scale(A3b$Range.Size..km2.)
A3b$L1_pop <- scale(A3b$L1_pop)
A3b$Distance.to.Mainland <- scale(A3b$Distance.to.Mainland)
A3b$Distance.to.Continent <- scale(A3b$Distance.to.Continent)

## --- A4 --- ####
A4$Phoneme.Inventory.Size <- log(A4$Phoneme.Inventory.Size) 

## --- A5a --- ####
A5a$Phoneme.Inventory.Size <- log(A5a$Phoneme.Inventory.Size) 
A5a$L1_pop <- log(A5a$L1_pop + 0.5)

## --- A5b ---  ####

# --- Log transform variables ---
A5b$Phoneme.Inventory.Size <- log(A5b$Phoneme.Inventory.Size) 
A5b$Range.Size..km2. <- log(A5b$Range.Size..km2.)
A5b$altitude_range <- log(A5b$altitude_range + 0.5)
A5b$L1_pop <- log(A5b$L1_pop + 0.5)
A5b$bordering_language_richness <- log(A5b$bordering_language_richness + 0.5)

# --- Standardise predictors ---
A5b$Range.Size..km2. <- scale(A5b$Range.Size..km2.)
A5b$altitude_range <- scale(A5b$altitude_range)
A5b$L1_pop <- scale(A5b$L1_pop)
A5b$bordering_language_richness <- scale(A5b$bordering_language_richness)

## --- A5c ---  ####

# --- Log transform variables ---
A5c$Phoneme.Inventory.Size <- log(A5c$Phoneme.Inventory.Size) 
A5c$Range.Size..km2. <- log(A5c$Range.Size..km2.)
A5c$L1_pop <- log(A5c$L1_pop + 0.5)
A5c$bordering_language_richness <- log(A5c$bordering_language_richness + 0.5)
# Do not transform zero entries for distance variables
A5c$Distance.to.Mainland[A5c$Distance.to.Mainland != 0] <- log(
  A5c$Distance.to.Mainland[A5c$Distance.to.Mainland != 0])
A5c$Distance.to.Continent[A5c$Distance.to.Continent != 0] <- log(
  A5c$Distance.to.Continent[A5c$Distance.to.Continent != 0])

# --- Standardise predictors ---
A5c$Range.Size..km2. <- scale(A5c$Range.Size..km2.)
A5c$L1_pop <- scale(A5c$L1_pop)
A5c$bordering_language_richness <- scale(A5c$bordering_language_richness)
A5c$Distance.to.Mainland <- scale(A5c$Distance.to.Mainland)
A5c$Distance.to.Continent <- scale(A5c$Distance.to.Continent)

## --- A5d --- ####
A5d$Phoneme.Inventory.Size <- log(A5d$Phoneme.Inventory.Size) 

## --- A6 --- ####
A6$L1_pop <- log(A6$L1_pop + 0.5) 

# --- (6) DATASET VISUALISATION --- ####

## Setup ####

# Ordering regions to have those close together kept in plots
region_order <- c("Oceania", "Australia and New Zealand", 
                  "South-Eastern Asia", "Southern Asia", "Asia", 
                  "Europe", "North Africa and Arabia", "Africa", 
                  "Western Africa", "Northern America", "Central America", "South America")

# Consistent colours scheme for scatter plots labelling regions
# colour palette from https://colorbrewer2.org
region_colours <- c("Oceania" = "#a6cee3", "Australia and New Zealand" = "#1f78b4",
                    "South-Eastern Asia" = "#b2df8a", "Southern Asia" = "#33a02c",
                    "Asia" = "#fb9a99", "Europe" = "#e31a1c",
                    "North Africa and Arabia" = "#fdbf6f", "Africa" = "#ff7f00",
                    "Western Africa" = "#cab2d6", "Northern America" = "#6a3d9a",
                    "Central America" = "#ffff99", "South America" = "#b15928")

# Retrieve nee22 (if not already loaded)
if(!exists("nee22")){
  nee22 <- read.csv("data/Global predictors of language endangerment and the future of linguistic diversity Data 2.csv",
                    stringsAsFactors = FALSE)
  # Change Region "Arab" to "North Africa and Arabia"
  nee22$region[nee22$region == "Arab"] <- "North Africa and Arabia"
}

## --- A1a --- ####

# --- Preview raw data A1a ---
raw_A1a <- read.csv("output/A1a/A1a_adjusted.csv")
head(raw_A1a)
summary(raw_A1a)
# Phoneme Inventory Size distribution
PIS_hist <- ggplot(raw_A1a, aes(x = Phoneme.Inventory.Size)) +
  geom_histogram() +
  labs(
    y = "Number of Languages",
    x = "PIS_A",
  ) + 
  scale_x_continuous(n.breaks = 15, expand = c(0, 0)) +
  scale_y_continuous(n.breaks = 15, expand = c(0, 0)) + 
  coord_cartesian(xlim = c(0, 150), ylim = c(0, 400)) +
  theme_classic()

ggsave(
  filename = "output/A1a/PIS_hist.png",
  plot = PIS_hist,
  scale = 1
)
print(PIS_hist)

# --- Mapping phoneme inventory sizes onto world map A1a ---
# Get world map data
world <- map_data("world")
world <- world[world$region != "Antarctica",]
PIS_world_map <- ggplot() +
  # Base world map
  geom_polygon(
    data = world,
    aes(x = long, y = lat, group = group),
    fill = "gray80", color = "gray80"
  ) +
  # Points with size and color based on count
  geom_point(
    data = A1a,
    # aes(x = Longitude, y = Latitude, color = Phoneme.Inventory.Size, size = L1_pop),
    aes(x = Longitude, y = Latitude, color = Phoneme.Inventory.Size),
    alpha = 0.7
  ) +
  scale_color_viridis_c(option="plasma") +
  coord_fixed(1.3) +
  theme_minimal(base_size = 20) +
  labs(x = NULL, y = NULL, color = "Log PISa", size = "Log L1_pop")
PIS_world_map
ggsave("output/A1a/PIS_world_map.png", plot = PIS_world_map, width = 15, height = 10, dpi = 300)

# --- Preview correlations A1a ---
# Plot data
phoneme_L1_scatter = ggplot(A1a, aes(x = L1_pop, y = Phoneme.Inventory.Size, color = factor(region, levels = region_order))) + 
  geom_point() + 
  labs(x = "L1_pop", y = "PISa", color = "Region") + 
  scale_color_manual(values = region_colours) +
  theme_minimal()
# Save graph
ggsave("output/A1a/phoneme_to_L1_scatter.png", plot = phoneme_L1_scatter, width = 15, height = 10, dpi = 100,scale = 0.5)
# Show graph
print(phoneme_L1_scatter)

# --- Testing for region sampling bias of A1a (BUT WILL BE USED IN SECTION 4 of PAPER) ---
# number of repetitions
n_rep = 1e4
# Find number of languages in each region for A1a
A1a_region_tally <- as.data.frame(table(A1a$region))
A1a_region_tally <- A1a_region_tally[match(region_order, A1a_region_tally$Var1), ]
region_sample_tallies <- data.frame(matrix(ncol = length(region_order), nrow = n_rep))
colnames(region_sample_tallies) <- region_order # Set column names
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
A1a_region_sample_tallies_long <- pivot_longer(region_sample_tallies, everything(), names_to = c("Group"), values_to = 'Value')
# Set 'Group' as a factor with levels in the desired order, this is to ensure that ggplot doesn't alphabetically sort the regions.
A1a_region_sample_tallies_long$Group <- factor(A1a_region_sample_tallies_long$Group, levels = region_order)
region_box <- ggplot(A1a_region_sample_tallies_long, aes(x = Group, y = Value)) +
  geom_boxplot(varwidth = TRUE) +
  geom_point(data = A1a_region_tally, 
             aes(x = region_order, y = Freq), 
             shape = 4,
             size = 4) + 
  geom_text(
    data = A1a_region_tally,
    aes(x = region_order, y = Freq, label = Freq),
    vjust = 0.1,
    hjust = 1.5,
    size = 3.5
  ) +
  ylab("Number of Languages Sampled") +
  xlab("Region") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.title.x.bottom = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y.left = element_text(size = 10))
# Save graph
ggsave("output/A1a/region_sampling.png", plot = region_box, width = 15, height = 15, dpi = 300, scale = 0.5)
print(region_box)

# Show diffs between number of observed languages sampled and the median, Q1, Q4, min and max of sampled regions from NEE22 (Fig 4c)
A1a_region_sample_diff <- data.frame(region_order, 
                                     sapply(region_order, function(x) 
                                       {median(A1a_region_sample_tallies_long$Value[A1a_region_sample_tallies_long == x])}),
                                     sapply(region_order, function(x) 
                                     {quantile(A1a_region_sample_tallies_long$Value[A1a_region_sample_tallies_long == x])[[2]]}),
                                     sapply(region_order, function(x) 
                                     {quantile(A1a_region_sample_tallies_long$Value[A1a_region_sample_tallies_long == x])[[4]]}),
                                     sapply(region_order, function(x) 
                                     {min(A1a_region_sample_tallies_long$Value[A1a_region_sample_tallies_long == x])}),
                                     sapply(region_order, function(x) 
                                     {max(A1a_region_sample_tallies_long$Value[A1a_region_sample_tallies_long == x])}),
                                     sapply(region_order, function(x) 
                                       {median(A1a$Phoneme.Inventory.Size[A1a$region == x])}),
                                     stringsAsFactors = FALSE)

colnames(A1a_region_sample_diff) <- c("region", "median", "Q1", "Q3", "min", "max", "median_PISa")
A1a_region_sample_diff$median <- A1a_region_tally$Freq - A1a_region_sample_diff$median
A1a_region_sample_diff$Q1 <- A1a_region_tally$Freq - A1a_region_sample_diff$Q1
A1a_region_sample_diff$Q3 <- A1a_region_tally$Freq - A1a_region_sample_diff$Q3
A1a_region_sample_diff$min <- A1a_region_tally$Freq - A1a_region_sample_diff$min
A1a_region_sample_diff$max <- A1a_region_tally$Freq - A1a_region_sample_diff$max

# Make into long format
A1a_region_sample_diff_long <- A1a_region_sample_diff %>%
  pivot_longer(
    cols = c("median", "Q1", "Q3", "min", "max"),
    names_to = "diff_type",
    values_to = "diff_value"
  )
A1a_region_sample_diff_long$diff_type <- factor(
  A1a_region_sample_diff_long$diff_type,
  levels = c("min", "Q1", "median", "Q3", "max"),
  labels = c("Min", "Q1", "Median", "Q3", "Max")
)

A1a_region_sample_plot <- ggplot(data=A1a_region_sample_diff_long, aes(y=diff_value,x=median_PISa,color=diff_type)) +
  geom_point() + 
  labs(x = "Median PISa", y = "Sample Bias", color = "Sample Bias Type") + 
  theme_classic()
print(A1a_region_sample_plot)
ggsave("output/A1a/region_sampling_diff.png",A1a_region_sample_plot)
# Take quick OLS
summary(lm(A1a_region_sample_diff_long$median_PISa ~ A1a_region_sample_diff_long$diff_value))

# --- Phoneme.Inventory.Size across regions A1a (BUT WILL BE USED IN SECTION 4 of PAPER) ---
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
  theme_classic() +
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
    y = "Number of Languages",
    x = "PIS_C",
  ) + 
  scale_x_continuous(n.breaks = 15, expand = c(0, 0)) +
  scale_y_continuous(n.breaks = 15, expand = c(0, 0)) + 
  coord_cartesian(xlim = c(0, 150), ylim = c(0, 400)) +
  theme_classic()

ggsave(
  filename = "output/A1b/PIS_hist.png",
  plot = PIS_hist,
  scale = 1
)
print(PIS_hist)

# --- Preview correlations A1b ---
# Plot data
phoneme_L1_scatter = ggplot(A1b, aes(x = L1_pop, y = Phoneme.Inventory.Size, color = factor(region, levels = region_order))) + 
  geom_point() + 
  labs(x = "L1_pop", y = "PISc", color = "Region") + 
  scale_color_manual(values = region_colours) +
  theme_minimal()
# theme(legend.position = "none") + # Hide the legend for now, often too many groups
# Save graph
ggsave("output/A1b/phoneme_to_L1_scatter.png", plot = phoneme_L1_scatter, width = 15, height = 10, dpi = 100,scale = 0.5)
# Show graph
print(phoneme_L1_scatter)

# --- Testing for region sampling bias A1b (BUT WILL BE USED IN SECTION 4 of PAPER) ---
# Retrieve nee22 (if not already loaded)
if(!exists("nee22")){
  nee22 <- read.csv("data/Global predictors of language endangerment and the future of linguistic diversity Data 2.csv",
                    stringsAsFactors = FALSE) 
}
# number of repetitions
n_rep = 1e4
# Get region names
region_order <- region_order
# Find number of languages in each region for A1b
A1b_region_tally <- as.data.frame(table(A1b$region))
A1b_region_tally <- A1b_region_tally[match(region_order, A1b_region_tally$Var1), ]
region_sample_tallies <- data.frame(matrix(ncol = length(region_order), nrow = n_rep))
colnames(region_sample_tallies) <- region_order # Set column names
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
A1b_region_sample_tallies_long <- pivot_longer(region_sample_tallies, everything(), names_to = c("Group"), values_to = 'Value')
# Set 'Group' as a factor with levels in the desired order, this is to ensure that ggplot doesn't alphabetically sort the regions.
A1b_region_sample_tallies_long$Group <- factor(A1b_region_sample_tallies_long$Group, levels = region_order)
region_box <- ggplot(A1b_region_sample_tallies_long, aes(x = Group, y = Value)) +
  geom_boxplot(varwidth = TRUE) +
  geom_point(data = A1b_region_tally, 
             aes(x = region_order, y = Freq), 
             shape = 4,
             size = 4) + 
  geom_text(
    data = A1b_region_tally,
    aes(x = region_order, y = Freq, label = Freq),
    vjust = 0.1,
    hjust = 1.4,
    size = 3
  ) +
  ylab("Number of Languages Sampled") +
  xlab("Region") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.title.x.bottom = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y.left = element_text(size = 10))
# Save graph
ggsave("output/A1b/region_sampling.png", plot = region_box, width = 15, height = 15, dpi = 300, scale = 0.5)
print(region_box)

# --- Phoneme.Inventory.Size across regions A1b (BUT WILL BE USED IN SECTION 4 of PAPER) ---
box_phoneme_regions <- ggplot(raw_A1b, aes(
  x = factor(region, levels = region_order),
  y = Phoneme.Inventory.Size
)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = median(raw_A1b$Phoneme.Inventory.Size), color = "red") +
  labs(
    x = "Region",
    y = "PISc"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(box_phoneme_regions)
ggsave("output/A1b/phoneme_region_box.png", plot = box_phoneme_regions, width = 15, height = 15, dpi = 300, scale = 0.5)

# Show diff between number of observed languages sampled and the median, Q1, Q4, min and max of sampled regions from NEE22 (Fig 4c)
A1b_region_sample_diff <- data.frame(region_order, 
                                     sapply(region_order, function(x) 
                                     {median(A1b_region_sample_tallies_long$Value[A1b_region_sample_tallies_long == x])}),
                                     sapply(region_order, function(x) 
                                     {quantile(A1b_region_sample_tallies_long$Value[A1b_region_sample_tallies_long == x])[[2]]}),
                                     sapply(region_order, function(x) 
                                     {quantile(A1b_region_sample_tallies_long$Value[A1b_region_sample_tallies_long == x])[[4]]}),
                                     sapply(region_order, function(x) 
                                     {min(A1b_region_sample_tallies_long$Value[A1b_region_sample_tallies_long == x])}),
                                     sapply(region_order, function(x) 
                                     {max(A1b_region_sample_tallies_long$Value[A1b_region_sample_tallies_long == x])}),
                                     sapply(region_order, function(x) 
                                     {median(A1b$Phoneme.Inventory.Size[A1b$region == x])}),
                                     stringsAsFactors = FALSE)

colnames(A1b_region_sample_diff) <- c("region", "median", "Q1", "Q3", "min", "max", "median_PISc")
A1b_region_sample_diff$median <- A1b_region_tally$Freq - A1b_region_sample_diff$median
A1b_region_sample_diff$Q1 <- A1b_region_tally$Freq - A1b_region_sample_diff$Q1
A1b_region_sample_diff$Q3 <- A1b_region_tally$Freq - A1b_region_sample_diff$Q3
A1b_region_sample_diff$min <- A1b_region_tally$Freq - A1b_region_sample_diff$min
A1b_region_sample_diff$max <- A1b_region_tally$Freq - A1b_region_sample_diff$max

# Make into long format
A1b_region_sample_diff_long <- A1b_region_sample_diff %>%
  pivot_longer(
    cols = c("median", "Q1", "Q3", "min", "max"),
    names_to = "diff_type",
    values_to = "diff_value"
  )
A1b_region_sample_diff_long$diff_type <- factor(
  A1b_region_sample_diff_long$diff_type,
  levels = c("min", "Q1", "median", "Q3", "max"),
  labels = c("Min", "Q1", "Median", "Q3", "Max")
)

A1b_region_sample_plot <- ggplot(data=A1b_region_sample_diff_long, aes(y=diff_value,x=median_PISc,color=diff_type)) +
  geom_point() + 
  labs(x = "Median PISc", y = "Sample Bias", color = "Sample Bias Type") + 
  theme_classic()
print(A1b_region_sample_plot)
ggsave("output/A1b/region_sampling_diff.png",A1b_region_sample_plot)
# Take quick OLS
summary(lm(A1b_region_sample_diff_long$median_PISc ~ A1b_region_sample_diff_long$diff_value))

## --- A2 --- ####

# Get raw A2 adjusted data and matrices
raw_A2 <- read.csv("output/A2/A2_adjusted.csv")
load("output/A2/A2_phylomatrix.RData")
load("output/A2/A2_spmatrix.RData")

# --- Preview transformed data ---
head(A2)
summary(A2)
PIS_Island <- ggplot(data = A2, aes(x = as.factor(Island.Endemic), y = Phoneme.Inventory.Size)) +
  geom_boxplot()
print(PIS_Island)

# --- Preview correlations ---
pairs(A2[,
         c("Phoneme.Inventory.Size", "Range.Size..km2.", "L1_pop", 
           "altitude_range", "bordering_language_richness")],
      panel = function(x, y, ...) {
        points(x, y, cex = 0.1, ...) # cex = 0.5 makes points half the default size
      })
cor(A2[,
       c("Phoneme.Inventory.Size", "Range.Size..km2.", "L1_pop", 
         "altitude_range", "bordering_language_richness")])
cor.test(A2$Phoneme.Inventory.Size,A2$Island.Endemic)

## --- A3a --- ####

# Get raw A3a adjusted data and matrices
raw_A3a <- read.csv("output/A3a/A3a_adjusted.csv")
load("output/A3a/A3a_phylomatrix.RData")
load("output/A3a/A3a_spmatrix.RData")

# --- Preview transformed data ---
head(A3a)
summary(A3a)

# --- PISc across regions and island endemism ---
counts_real <- A3a %>%
  group_by(region, Island.Endemic) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(Island.Endemic = factor(Island.Endemic, levels = c(0,1)))  # make factor

# Create all combinations and merge
PIS_Island_Region_counts <- expand_grid(
  region = factor(region_order, levels = region_order),
  Island.Endemic = factor(c(0,1), levels = c(0,1))
) %>%
  left_join(counts_real, by = c("region","Island.Endemic")) %>%
  mutate(
    n = ifelse(is.na(n), 0, n),  # missing combos get 0
    x_comb = factor(paste0(region, "_", Island.Endemic),
                    levels = paste0(rep(region_order, each = 2), "_", 0:1))
  )
A3a_expanded <- A3a %>%
  mutate(
    region = factor(region, levels = region_order),
    Island.Endemic = factor(Island.Endemic, levels = c(0,1))
  ) %>%
  tidyr::complete(region, Island.Endemic, fill = list(Phoneme.Inventory.Size = NA)) %>%
  arrange(region, Island.Endemic) %>%
  mutate(
    x_comb = factor(paste0(region, "_", Island.Endemic),
                    levels = paste0(rep(region_order, each = 2), "_", 0:1))
  )
x_labels <- rep("", length(levels(A3a_expanded$x_comb)))
x_labels[seq(1, length(x_labels), by = 2)] <- region_order
# Step 3: Plot
PIS_Island_Region <- ggplot(A3a_expanded, aes(x = x_comb, y = Phoneme.Inventory.Size, fill = Island.Endemic)) +
  geom_boxplot(na.rm = TRUE, position = position_dodge(width = 0.8)) +
  geom_text(
    data = PIS_Island_Region_counts,
    aes(x = x_comb, y = max(A3a_expanded$Phoneme.Inventory.Size, na.rm = TRUE) * 0.9, 
        label = n),
    position = position_dodge(width = 0.8),
    size = 3
  ) +
  geom_hline(yintercept=median(A3a$Phoneme.Inventory.Size)) + 
  scale_fill_manual(
    values = c("0" = "#1f78b4", "1" = "#fb9a99"),
    labels = c("0" = "Not Island Endemic", "1" = "Island Endemic"),
    name = ""
  ) +
  labs(x = "Region", y = "Phoneme Inventory Size", fill = "Island Endemic") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.ticks.x=element_blank()
  ) +
  scale_x_discrete(labels = x_labels)

print(PIS_Island_Region)
ggsave(
  filename = "output/A3a/PIS_Island_Region.png",
  plot = PIS_Island_Region,
  scale = 1,
  width=8,
  height=5
)

# --- PISc ~ Area ---
# Regions labelled
PIS_Range_Region <- ggplot(A3a, aes(x = Range.Size..km2., 
                                    y = Phoneme.Inventory.Size, 
                                    color = factor(region, levels = region_order))) + 
  geom_point() + 
  labs(x = "Range Size", y = "PISc", color = "Region") + 
  scale_color_manual(values = region_colours) +
  theme_minimal()
print(PIS_Range_Region)
ggsave(
  filename = "output/A3a/PIS_Range_Region.png",
  plot = PIS_Range_Region,
  scale = 1,
  width=7,
  height=5
)
# Island
PIS_Range_Island <- ggplot(A3a, aes(x = Range.Size..km2., 
                                    y = Phoneme.Inventory.Size, 
                                    color = factor(Island.Endemic))) + 
  geom_point() + 
  scale_color_manual(
    values = c("0" = "#1f78b4", "1" = "#fb9a99"),
    labels = c("0" = "Not Island Endemic", "1" = "Island Endemic"),
    name = ""
  ) + 
  labs(x = "Range Size (km^2)", y = "Phoneme Inventory Size", color = "Island Endemic", ) +
  theme_classic()
print(PIS_Range_Island)
ggsave(
  filename = "output/A3a/PIS_Range_Island.png",
  plot = PIS_Range_Island,
  scale = 1,
  width=7,
  height=5
)

# --- Get PISc ~ Area OLS, for each region ---
PIS_Range_Region_Cor <- sapply(region_order, function(x) {
  cor(A3a$Range.Size..km2.[A3a$region==x], A3a$Phoneme.Inventory.Size[A3a$region==x])
})

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

# --- Preview transformed data ---
head(A3b)
summary(A3b)
# PIS ~ Range_Size
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

# --- Preview transformed data ---
head(A4)
summary(A4)
print("Number of languages for each documentation level:")
print(table(A4$documentation))

# --- Preview correlations ---
# PIS ~ Documentation
PIS_Doc <- ggplot(data = A4, aes(x = as.factor(documentation), y = Phoneme.Inventory.Size)) +
  geom_boxplot() +
  labs(x = "Documentation", y = "PISc") + 
  theme_classic()
print(PIS_Doc)
ggsave(filename = "output/A4/PIS_Doc.png", plot = PIS_Doc)

## --- A5a --- ####
# --- Mapping phoneme inventory sizes onto world map A1a ---
# Get world map data
PIS_world_map_2 <- ggplot() +
  # Base world map
  geom_polygon(
    data = world,
    aes(x = long, y = lat, group = group),
    fill = "gray90", color = "gray50"
  ) +
  # Points with size and color based on count
  geom_point(
    data = A5a,
    aes(x = Longitude, y = Latitude, color = Phoneme.Inventory.Size),
    alpha = 0.7
  ) +
  scale_color_viridis_c(option="plasma") +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(x = "", y = "", color = "Log PISa")
ggsave("output/A5a/PIS_world_map_2.png", plot = PIS_world_map_2, width = 15, height = 10, dpi = 300)
PIS_world_map_2

# --- Testing for region sampling bias of A5a (applies for A5b,A5c and A5d too) ---
# number of repetitions
n_rep = 1e4
# Find number of languages in each region for A5a
A5a_region_tally <- as.data.frame(table(A5a$region))
A5a_region_tally <- A5a_region_tally[match(region_order, A5a_region_tally$Var1), ]
region_sample_tallies <- data.frame(matrix(ncol = length(region_order), nrow = n_rep))
colnames(region_sample_tallies) <- region_order # Set column names
# Repeated n_rep times, randomly sample without replacement 
for(i in 1:n_rep){
  # Get region sample
  region_sample <- sample(nee22$region, size = length(A5a$region), replace = FALSE)
  # Current region tally.
  region_sample_tally = table(region_sample)
  # Add current tally to all tally data
  region_sample_tallies[i, names(region_sample_tally)] <- as.numeric(region_sample_tally)
}
# Reshape the dataframe to long format
A5a_region_sample_tallies_long <- pivot_longer(region_sample_tallies, everything(), names_to = c("Group"), values_to = 'Value')
# Set 'Group' as a factor with levels in the desired order, this is to ensure that ggplot doesn't alphabetically sort the regions.
A5a_region_sample_tallies_long$Group <- factor(A5a_region_sample_tallies_long$Group, levels = region_order)
region_box <- ggplot(A5a_region_sample_tallies_long, aes(x = Group, y = Value)) +
  geom_boxplot(varwidth = TRUE) +
  geom_point(data = A5a_region_tally, 
             aes(x = region_order, y = Freq), 
             shape = 4,
             size = 4) + 
  geom_text(
    data = A5a_region_tally,
    aes(x = region_order, y = Freq, label = Freq),
    vjust = 0.1,
    hjust = 1.5,
    size = 3.5
  ) +
  ylab("Number of Languages Sampled") +
  xlab("Region") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.title.x.bottom = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y.left = element_text(size = 10))
# Save graph
ggsave("output/A5a/region_sampling.png", plot = region_box, width = 15, height = 15, dpi = 300, scale = 0.5)
print(region_box)

## --- A6 --- ####

# Get raw A4 adjusted data and matrices
raw_A6 <- read.csv("output/A6/A6_adjusted.csv")
load("output/A6/A6_phylomatrix.RData")
load("output/A6/A6_spmatrix.RData")

# --- Preview transformed data ---
head(A6)
summary(A6)
print("Number of languages for each documentation level:")
print(table(A6$documentation))

# --- Preview correlations ---
# L1_pop ~ Documentation
L1_pop_Doc <- ggplot(data = A6, aes(x = as.factor(documentation), y = L1_pop)) +
  geom_boxplot() +
  labs(x = "Documentation", y = "L1_pop") + 
  theme_classic()
print(L1_pop_Doc)
ggsave(filename = "output/A6/L1_pop_Doc.png", plot = L1_pop_Doc)

# --- (7) OLS AND GLS ANALYSIS --- ####

## --- A1a --- ####

# Models for Phoneme.Inventory.Size ~ L1_pop

# --- Ordinary Least Squares A1a ---
A1a_OLS <- lm(formula = Phoneme.Inventory.Size ~ L1_pop, data = A1a)
A1a_OLS_null <- lm(formula = Phoneme.Inventory.Size ~ 1, data = A1a)
print(summary(A1a_OLS))
print(BIC(A1a_OLS))
print(logLik(A1a_OLS))

# --- Best p values for model A1a (S+P) ---
A1a_GLS_p.res <- sbplx(c(0.5, 0.5, 0.5),
                   best_p,
                   formula=Phoneme.Inventory.Size ~ L1_pop,
                   data=A1a,
                   spmatrix=A1a_spmatrix,phylomatrix=A1a_phylomatrix,
                   lower=c(0,0,0),upper=c(1,1,1),
                   nl.info = TRUE)
# Save p.res
save(A1a_GLS_p.res, file = "output/A1a/A1a_GLS_p_res.RData")

# --- Maximum likelihood fits A1a GLS (S+P) ---
A1a_GLS_model <- ml_fit(p=A1a_GLS_p.res$par,formula=Phoneme.Inventory.Size ~ L1_pop,data=A1a,spmatrix=A1a_spmatrix,phylomatrix=A1a_phylomatrix)
# Save model fit
save(A1a_GLS_model, file = "output/A1a/A1a_GLS_model.RData") 

# --- Best p values for model A1a (S) ---
A1a_GLS_S_p.res <- sbplx(c(0.5, 0.5, 0.5),
                       best_p,
                       formula=Phoneme.Inventory.Size ~ L1_pop,
                       data=A1a,
                       spmatrix=A1a_spmatrix,phylomatrix=A1a_phylomatrix,
                       opt = "S",
                       lower=c(0,0,0),upper=c(1,1,1),
                       nl.info = TRUE)
# Save p.res
save(A1a_GLS_S_p.res, file = "output/A1a/A1a_GLS_S_p_res.RData")

# --- Maximum likelihood fits A1a GLS (S) ---
A1a_GLS_S_model <- ml_fit(p=A1a_GLS_S_p.res$par,formula=Phoneme.Inventory.Size ~ L1_pop,data=A1a,spmatrix=A1a_spmatrix,phylomatrix=A1a_phylomatrix,opt="S")
# Save model fit
save(A1a_GLS_S_model, file = "output/A1a/A1a_GLS_S_model.RData") 

# --- Best p values for model A1a (P) ---
A1a_GLS_P_p.res <- sbplx(c(0.5, 0.5, 0.5),
                         best_p,
                         formula=Phoneme.Inventory.Size ~ L1_pop,
                         data=A1a,
                         spmatrix=A1a_spmatrix,phylomatrix=A1a_phylomatrix,
                         opt = "P",
                         lower=c(0,0,0),upper=c(1,1,1),
                         nl.info = TRUE)
# Save p.res
save(A1a_GLS_P_p.res, file = "output/A1a/A1a_GLS_P_p_res.RData")

# --- Maximum likelihood fits A1a GLS (P) ---
A1a_GLS_P_model <- ml_fit(p=A1a_GLS_P_p.res$par,formula=Phoneme.Inventory.Size ~ L1_pop,data=A1a,spmatrix=A1a_spmatrix,phylomatrix=A1a_phylomatrix,opt="P")
# Save model fit
save(A1a_GLS_P_model, file = "output/A1a/A1a_GLS_P_model.RData")

# --- Best p values for null of model A1a GLS ---
# A1a_GLS_null_p.res <- sbplx(c(0.5, 0.5, 0.5),
#                    best_p,
#                    formula=Phoneme.Inventory.Size ~ 1,
#                    data=A1a,
#                    spmatrix=A1a_spmatrix,phylomatrix=A1a_phylomatrix,
#                    lower=c(0,0,0),upper=c(1,1,1),
#                    nl.info = TRUE)
# # Save p.res
# save(A1a_GLS_null_p.res, file = "output/A1a/A1a_GLS_null_p_res.RData")

# --- Nagelkerke/Cragg-Uhler Pseudo-R^2 of OLS and GLS ---
# (1) Using OLS null against full OLS
A1a_OLS_nll <- - logLik(A1a_OLS) %>% as.numeric # Negative log likelihood of full model
A1a_OLS_null_nll <- -logLik(A1a_OLS_null) %>% as.numeric # Negative loglikelihood of null model
A1a_OLS_R2 <- (1-exp(A1a_OLS_nll - A1a_OLS_null_nll)^(2/dim(A1a)[1]))/(1-exp(-A1a_OLS_null_nll)^(2/dim(A1a)[1]))
# (2) Using OLS null against GLS (S+P)
A1a_GLS_SP_nll <- A1a_GLS_p.res$value # Negative log likelihood of full model
A1a_GLS_SP_R2 <-  (1-exp(A1a_GLS_SP_nll-A1a_OLS_null_nll)^(2/dim(A1a)[1]))/(1-exp(-A1a_OLS_null_nll)^(2/dim(A1a)[1]))
# (3) Using OLS null against GLS (S)
A1a_GLS_S_nll <- A1a_GLS_S_p.res$value # Negative log likelihood of full model
A1a_GLS_S_R2 <-  (1-exp(A1a_GLS_S_nll-A1a_OLS_null_nll)^(2/dim(A1a)[1]))/(1-exp(-A1a_OLS_null_nll)^(2/dim(A1a)[1]))
# (4) Using OLS null against GLS (P)
A1a_GLS_P_nll <- A1a_GLS_P_p.res$value # Negative log likelihood of full model
A1a_GLS_P_R2 <-  (1-exp(A1a_GLS_P_nll-A1a_OLS_null_nll)^(2/dim(A1a)[1]))/(1-exp(-A1a_OLS_null_nll)^(2/dim(A1a)[1]))

## --- A1b --- ####

# --- Ordinary Least Squares A1b ---
A1b_l1 <- lm(formula = Phoneme.Inventory.Size ~ L1_pop, data = A1b)
print(summary(A1b_l1))
print(BIC(A1b_l1))
print(logLik(A1b_l1))

# --- Best p values for all models A1b ---
A1b_p.res <- sbplx(c(0.5, 0.5, 0.5),
                   best_p,
                   formula=Phoneme.Inventory.Size ~ L1_pop,
                   data=A1b,
                   spmatrix=A1b_spmatrix,phylomatrix=A1b_phylomatrix,
                   lower=c(0,0,0),upper=c(1,1,1),
                   nl.info = TRUE)
# Save p.res
save(A1b_p.res, file = "output/A1b/A1b_p_res.RData")

# --- Maximum likelihood fits for all models A1b ---
A1b_model <- ml_fit(p=A1b_p.res$par,formula=Phoneme.Inventory.Size ~ L1_pop,data=A1b,spmatrix=A1b_spmatrix,phylomatrix=A1b_phylomatrix)
# Save model fit
save(A1b_model, file = "output/A1b/A1b_model.RData") 

## --- A2 --- ####

# --- Ordinary Least Squares ---
A2_l1 <- lm(formula = Phoneme.Inventory.Size ~ L1_pop, data = A2)
print(summary(A2_l1))
print(BIC(A2_l1))
print(logLik(A2_l1))

# --- Model predictor combinations ---
A2_predictors <- c("Island.Endemic", "Range.Size..km2.",
                "altitude_range", "bordering_language_richness",
                "L1_pop")
response <- "Phoneme.Inventory.Size"
n_pred <- length(A2_predictors)
pred_combs <- sapply(1:n_pred, function(x) combn(A2_predictors, x))
A2_models <- c(Phoneme.Inventory.Size ~ 1) # List of models
A2_models_var <- list() # List of predictors for each model
A2_models_var[[1]] <- character(0) # Intercept model
k <- 2
for (i in 1:n_pred) {
  for (j in 1:dim(pred_combs[[i]])[2]) {
    combination <- pred_combs[[i]][,j]
    A2_models_var[[k]] <- combination
    model <- (paste(response, paste(combination, collapse="+"), sep="~"))
    model <- as.formula(model)
    A2_models <- c(A2_models,model)
    k <- k + 1
  }
}
unlist(A2_models)

# --- Linear regression fits for all models (no correction for autocorrelation) ---
A2_fits_nc <- mclapply(
  X = A2_models,
  FUN = function(f){lm(data = A2, formula = f)},
  mc.cores = 8
)
# Save models
saveRDS(A2_fits_nc, "output/A2/A2_fits_nc.RDS") 

# --- Model Comparison (for OLS) ---
# Make empty NA list of coefficients and their respective p values in tuple form
n <- length(A2_fits_nc)
A2_summaries_nc <- list("Island.Endemic" = as.list(rep(NA,n)),
                     "Range.Size..km2." = as.list(rep(NA,n)),
                     "altitude_range" = as.list(rep(NA,n)), 
                     "bordering_language_richness" = as.list(rep(NA,n)),
                     "L1_pop" = as.list(rep(NA,n)))
ml_logLik <- c()
ml_BIC <- c()
ml_AIC <- c()
ml_MSE <- c()

for (i in 1:n){
  # get current model summary
  ml_summary <- summary(A2_fits_nc[[i]])
  # if intercept only model just get BIC, logLik and MLE
  if(dim(ml_summary$coefficients)[1] == 1) {
    ml_logLik <- c(ml_logLik, logLik(A2_fits_nc[[i]]) %>% as.numeric())
    ml_BIC <- c(ml_BIC, BIC(A2_fits_nc[[i]]))
    ml_AIC <- c(ml_AIC, AIC(A2_fits_nc[[i]]))
    ml_MSE <- c(ml_MSE, sum(A2_fits_nc[[i]]$residuals**2)/dim(A2)[1])
    next
  }
  # Get coefficient values
  coef_val <- ml_summary$coefficients[,1]
  # Get corresponding p-values
  coef_p <- ml_summary$coefficients[,4]
  # Get coefficient names
  coef_names <- names(ml_summary$coefficients[,4])
  # Get number of coefficients
  n_coef <- length(coef_names)
  # Add coef info into model_summaries
  # for (j in 2:n_coef){
  #   if(coef_p[j]<=0.05){
  #     A2_summaries[[coef_names[j]]][[i]] <- paste0(round(coef_val[j], 5),"***")
  #   } else {
  #     A2_summaries[[coef_names[j]]][[i]] <- paste0(round(coef_val[j], 5))
  #   }
  # }
  for (j in 2:n_coef){
    A2_summaries_nc[[coef_names[j]]][[i]] <- paste0(round(coef_val[j],4),",",round(coef_p[j],4))
  }
  # Add BIC, logLik and MLE
  ml_logLik <- c(ml_logLik, logLik(A2_fits_nc[[i]]) %>% as.numeric())
  ml_BIC <- c(ml_BIC, BIC(A2_fits_nc[[i]]))
  ml_AIC <- c(ml_AIC, AIC(A2_fits_nc[[i]]))
  ml_MSE <- c(ml_MSE, sum(A2_fits_nc[[i]]$residuals**2)/dim(A2)[1])
}

# Create new data frame from lists
A2_ml_data_nc <- data.frame(BIC = ml_BIC,
                         AIC = ml_AIC,
                         logLik = ml_logLik,
                         MSE = ml_MSE,
                         unlist(A2_summaries_nc$Island.Endemic),
                         unlist(A2_summaries_nc$Range.Size..km2.),
                         unlist(A2_summaries_nc$altitude_range),
                         unlist(A2_summaries_nc$bordering_language_richness),
                         unlist(A2_summaries_nc$L1_pop),
                         stringsAsFactors = FALSE)
colnames(A2_ml_data_nc) <- c("BIC", "AIC", "logLik", "MSE", A2_predictors)
# Sort based on increasing BIC
A2_ml_data_BIC_nc <- A2_ml_data_nc[order(A2_ml_data_nc$BIC), ]
# Add delta BIC column (min_BIC - current_BIC)
min_BIC <-  min(A2_ml_data_BIC_nc$BIC)
A2_ml_data_BIC_nc$deltaBIC <- A2_ml_data_BIC_nc$BIC - min_BIC
# Save output
write.csv(A2_ml_data_BIC_nc, file = "output/A2/A2_ml_data_BIC_nc.csv", row.names=FALSE)
# Sort based on increasing AIC
A2_ml_data_AIC_nc <- A2_ml_data_nc[order(A2_ml_data_nc$AIC), ]
# Add delta AIC column (min_AIC - current_AIC)
min_AIC <-  min(A2_ml_data_AIC_nc$AIC)
A2_ml_data_AIC_nc$deltaAIC <- A2_ml_data_AIC_nc$AIC - min_AIC
# Save output
write.csv(A2_ml_data_AIC_nc, file = "output/A2/A2_ml_data_AIC_nc.csv", row.names=FALSE)

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

# --- BMA and PIP for all models ---
logp <- sapply(A2_fits,function (i) -BIC(i)/2)
pp <- exp(logp-max(logp)) 
pp <- pp/sum(pp) 
beta <- lapply(1:length(A2_fits),function (j) {
  coef <- A2_fits[[j]]$coefficients[-1] #assuming intercept is the first element in coefficients
  names(coef) <- A2_models_var[[j]]
  coef
})
A2_bma <- A2_pip <- rep(NA,length(A2_predictors))
names(A2_bma) <- names(A2_pip) <- A2_predictors
for (var in A2_predictors) {
  idx <- sapply(1:length(A2_fits),function (i) is.element(var,A2_models_var[[i]]))
  beta_var <- sapply(which(idx),function (i) beta[[i]][names(beta[[i]])==var])
  A2_bma[var] <- sum(beta_var*pp[idx])
  A2_pip[var] <- sum(pp[idx])
}

# --- Model Comparison ---
# Make empty NA list of coefficients and their respective p values in tuple form
n <- length(A2_fits)
A2_summaries <- list("Island.Endemic" = as.list(rep(NA,n)),
                     "Range.Size..km2." = as.list(rep(NA,n)),
                     "altitude_range" = as.list(rep(NA,n)), 
                     "bordering_language_richness" = as.list(rep(NA,n)),
                     "L1_pop" = as.list(rep(NA,n)))
ml_logLik <- c()
ml_BIC <- c()
ml_AIC <- c()
ml_MSE <- c()

for (i in 1:n){
  # get current model summary
  ml_summary <- summary(A2_fits[[i]])
  # if intercept only model just get BIC, logLik and MLE
  if(dim(ml_summary$tTable)[1] == 1) {
    ml_logLik <- c(ml_logLik, ml_summary$logLik)
    ml_BIC <- c(ml_BIC, ml_summary$BIC)
    ml_AIC <- c(ml_AIC, ml_summary$AIC)
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
  # for (j in 2:n_coef){
  #   if(coef_p[j]<=0.05){
  #     A2_summaries[[coef_names[j]]][[i]] <- paste0(round(coef_val[j], 5),"***")
  #   } else {
  #     A2_summaries[[coef_names[j]]][[i]] <- paste0(round(coef_val[j], 5))
  #   }
  # }
  for (j in 2:n_coef){
    A2_summaries[[coef_names[j]]][[i]] <- paste0(round(coef_val[j],4),",",round(coef_p[j],4))
  }
  # Add BIC, logLik and MLE
  ml_logLik <- c(ml_logLik, ml_summary$logLik)
  ml_BIC <- c(ml_BIC, ml_summary$BIC)
  ml_AIC <- c(ml_AIC, ml_summary$AIC)
  ml_MSE <- c(ml_MSE, sum(A2_fits[[i]]$residuals**2)/dim(A2)[1])
}

# Create new data frame from lists
A2_ml_data <- data.frame(BIC = ml_BIC,
                      AIC = ml_AIC,
                      logLik = ml_logLik,
                      MSE = ml_MSE,
                      unlist(A2_summaries$Island.Endemic),
                      unlist(A2_summaries$Range.Size..km2.),
                      unlist(A2_summaries$altitude_range),
                      unlist(A2_summaries$bordering_language_richness),
                      unlist(A2_summaries$L1_pop),
                      stringsAsFactors = FALSE)
colnames(A2_ml_data) <- c("BIC", "AIC", "logLik", "MSE", A2_predictors)
# Sort based on increasing BIC
A2_ml_data_BIC <- A2_ml_data[order(A2_ml_data$BIC), ]
# Add delta BIC column (min_BIC - current_BIC)
min_BIC <-  min(A2_ml_data_BIC$BIC)
A2_ml_data_BIC$deltaBIC <- A2_ml_data_BIC$BIC - min_BIC
# Save output
write.csv(A2_ml_data_BIC, file = "output/A2/A2_ml_data_BIC.csv", row.names=FALSE)
# Sort based on increasing AIC
A2_ml_data_AIC <- A2_ml_data[order(A2_ml_data$AIC), ]
# Add delta AIC column (min_AIC - current_AIC)
min_AIC <-  min(A2_ml_data_AIC$AIC)
A2_ml_data_AIC$deltaAIC <- A2_ml_data_AIC$AIC - min_AIC
# Save output
write.csv(A2_ml_data_AIC, file = "output/A2/A2_ml_data_AIC.csv", row.names=FALSE)

## --- A3a --- ####

# --- Get PISc ~ L1_pop OLS ---
A3a_l1 <- lm(formula = Phoneme.Inventory.Size ~ L1_pop, data = A3a)

# --- Get PISc ~ Area OLS ---
A3a_PIS_Range_l1 = lm(data=A3a, Phoneme.Inventory.Size ~ Range.Size..km2.)
summary(A3a_PIS_Range_l1)
# For each region
A3a_PIS_Range_Region_lm <- lapply(region_order, function(x) {
  summary(lm(data=A3a, A3a$Phoneme.Inventory.Size[A3a$region==x] ~ A3a$Range.Size..km2.[A3a$region==x]))
})
# For each region with only non-island languages
A3a_PIS_NI_Range_Region_lm <- lapply(region_order, function(x) {
  summary(lm(data=A3a,
             A3a$Phoneme.Inventory.Size[A3a$region==x & A3a$Island.Endemic == 0] ~ A3a$Range.Size..km2.[A3a$region==x  & A3a$Island.Endemic == 0]))
})
# For only Oceania region island languages
A3a_PIS_Range_l2 = lm(data=A3a[A3a$region == "Oceania" & A3a$Island.Endemic == 1,], Phoneme.Inventory.Size ~ Range.Size..km2.)

# --- Model predictor combinations ---
A3a_predictors <- c("Island.Endemic" ,"Range.Size..km2.",
                "Distance.to.Mainland", "Distance.to.Continent",
                "L1_pop")
response <- "Phoneme.Inventory.Size"
n_pred <- length(A3a_predictors)
pred_combs <- sapply(1:n_pred, function(x) combn(A3a_predictors, x))
A3a_models <- c(Phoneme.Inventory.Size~1)
A3a_models_var <- list() # List of predictors for each model
A3a_models_var[[1]] <- character(0) # Intercept model
k <- 2
for (i in 1:n_pred) {
  #print(pred_combs[[i]])
  for (j in 1:dim(pred_combs[[i]])[2]) {
    combination <- pred_combs[[i]][,j]
    A3a_models_var[[k]] <- combination
    model <- (paste(response, paste(combination, collapse="+"), sep="~"))
    model <- as.formula(model)
    A3a_models <- c(A3a_models,model)
    k <- k + 1
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

# --- BMA and PIP for all models ---
logp <- sapply(A3a_fits,function (i) -BIC(i)/2)
pp <- exp(logp-max(logp)) 
pp <- pp/sum(pp) 
beta <- lapply(1:length(A3a_fits),function (j) {
  coef <- A3a_fits[[j]]$coefficients[-1] #assuming intercept is the first element in coefficients
  names(coef) <- A3a_models_var[[j]]
  coef
})
A3a_bma <- A3a_pip <- rep(NA,length(A2_predictors))
names(A3a_bma) <- names(A3a_pip) <- A3a_predictors
for (var in A3a_predictors) {
  idx <- sapply(1:length(A3a_fits),function (i) is.element(var,A3a_models_var[[i]]))
  beta_var <- sapply(which(idx),function (i) beta[[i]][names(beta[[i]])==var])
  A3a_bma[var] <- sum(beta_var*pp[idx])
  A3a_pip[var] <- sum(pp[idx])
}

# --- Model Comparison ---
# Make empty NA list of coefficients and their respective p values in tuple form
n <- length(A3a_fits)
A3a_summaries <- list("Island.Endemic" = as.list(rep(NA,n)),
                     "Range.Size..km2." = as.list(rep(NA,n)),
                     "Distance.to.Mainland" = as.list(rep(NA,n)), 
                     "Distance.to.Continent" = as.list(rep(NA,n)),
                     "L1_pop" = as.list(rep(NA,n)))
ml_logLik <- c()
ml_BIC <- c()
ml_AIC <- c()
ml_MSE <- c()

for (i in 1:n){
  # get current model summary
  ml_summary <- summary(A3a_fits[[i]])
  # if intercept only model just get BIC, logLik and MLE
  if(dim(ml_summary$tTable)[1] == 1) {
    ml_logLik <- c(ml_logLik, ml_summary$logLik)
    ml_BIC <- c(ml_BIC, ml_summary$BIC)
    ml_AIC <- c(ml_AIC, ml_summary$AIC)
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
  # for (j in 2:n_coef){
  #   if(coef_p[j]<=0.05){
  #     A3a_summaries[[coef_names[j]]][[i]] <- paste0(round(coef_val[j], 5),"***")
  #   } else {
  #     A3a_summaries[[coef_names[j]]][[i]] <- paste0(round(coef_val[j], 5))
  #   }
  # }
  for (j in 2:n_coef){
    A3a_summaries[[coef_names[j]]][[i]] <- paste0(round(coef_val[j],4),",",round(coef_p[j],4))
  }
  # Add BIC, logLik and MLE
  ml_logLik <- c(ml_logLik, ml_summary$logLik)
  ml_BIC <- c(ml_BIC, ml_summary$BIC)
  ml_AIC <- c(ml_AIC, ml_summary$AIC)
  ml_MSE <- c(ml_MSE, sum(A3a_fits[[i]]$residuals**2)/dim(A3a)[1])
}

# Create new data frame from lists
A3a_ml_data <- data.frame(BIC = ml_BIC,
                      AIC = ml_AIC,
                      logLik = ml_logLik,
                      MSE = ml_MSE,
                      unlist(A3a_summaries$Island.Endemic),
                      unlist(A3a_summaries$Range.Size..km2.),
                      unlist(A3a_summaries$Distance.to.Mainland),
                      unlist(A3a_summaries$Distance.to.Continent),
                      unlist(A3a_summaries$L1_pop),
                      stringsAsFactors = FALSE)
colnames(A3a_ml_data) <- c("BIC", "AIC", "logLik", "MSE", A3a_predictors)
# Sort based on increasing BIC
A3a_ml_data_BIC <- A3a_ml_data[order(A3a_ml_data$BIC), ]
# Add delta BIC column (min_BIC - current_BIC)
min_BIC <-  min(A3a_ml_data_BIC$BIC)
A3a_ml_data_BIC$deltaBIC <- A3a_ml_data_BIC$BIC - min_BIC
# Save output
write.csv(A3a_ml_data_BIC, file = "output/A3a/A3a_ml_data_BIC.csv", row.names=FALSE)
# Sort based on increasing AIC
A3a_ml_data_AIC <- A3a_ml_data[order(A3a_ml_data$AIC), ]
# Add delta AIC column (min_AIC - current_AIC)
min_AIC <-  min(A3a_ml_data$AIC)
A3a_ml_data_AIC$deltaAIC <- A3a_ml_data_AIC$AIC - min_AIC
# Save output
write.csv(A3a_ml_data_AIC, file = "output/A3a/A3a_ml_data_AIC.csv", row.names=FALSE)

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
      A3b_summaries[[coef_names[j]]][[i]] <- paste0(round(coef_val[j], 4),"***")
    } else {
      A3b_summaries[[coef_names[j]]][[i]] <- paste0(round(coef_val[j], 4))
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

## --- A5a --- ####

# --- Ordinary Least Squares A5a ---
A5a_l1 <- lm(formula = Phoneme.Inventory.Size ~ L1_pop, data = A5a)
print(summary(A5a_l1))
print(BIC(A5a_l1))
print(logLik(A5a_l1))

# --- Best p values for all models A5a ---
A5a_p.res <- sbplx(c(0.5, 0.5, 0.5),
                   best_p,
                   formula=Phoneme.Inventory.Size ~ L1_pop,
                   data=A5a,
                   spmatrix=A5a_spmatrix,phylomatrix=A5a_phylomatrix,
                   lower=c(0,0,0),upper=c(1,1,1),
                   nl.info = TRUE)
# Save p.res
save(A5a_p.res, file = "output/A5a/A5a_p_res.RData")

# --- Maximum likelihood fits for all models A5a ---
A5a_model <- ml_fit(p=A5a_p.res$par,formula=Phoneme.Inventory.Size ~ L1_pop,data=A5a,spmatrix=A5a_spmatrix,phylomatrix=A5a_phylomatrix)
# Save model fit
save(A5a_model, file = "output/A5a/A5a_model.RData")

## --- A5b --- ####

# --- Ordinary Least Squares ---
A5b_l1 <- lm(formula = Phoneme.Inventory.Size ~ L1_pop, data = A5b)
print(summary(A5b_l1))
print(BIC(A5b_l1))
print(logLik(A5b_l1))

# --- Model predictor combinations ---
A5b_predictors <- c("Island.Endemic", "Range.Size..km2.",
                   "altitude_range", "bordering_language_richness",
                   "L1_pop")
response <- "Phoneme.Inventory.Size"
n_pred <- length(A5b_predictors)
pred_combs <- sapply(1:n_pred, function(x) combn(A5b_predictors, x))
A5b_models <- c(Phoneme.Inventory.Size ~ 1)
for (i in 1:n_pred) {
  for (j in 1:dim(pred_combs[[i]])[2]) {
    combination <- pred_combs[[i]][,j]
    model <- (paste(response, paste(combination, collapse="+"), sep="~"))
    model <- as.formula(model)
    A5b_models <- c(A5b_models,model)
  }
}
unlist(A5b_models)

# --- Best p values for all models ---
A5b_p.res <- sbplx(c(0.5, 0.5, 0.5),
                  best_p,
                  formula=Phoneme.Inventory.Size~1,
                  data=A5b,
                  spmatrix=A5b_spmatrix,phylomatrix=A5b_phylomatrix,
                  lower=c(0,0,0),upper=c(1,1,1),
                  nl.info = TRUE)
# Save p.res and covariance matrix
save(A5b_p.res, file = "output/A5b/A5b_p_res.RData")
spmatrix_temp <- A5b_spmatrix/max(A5b_spmatrix)
spmatrix_temp <- exp(-(spmatrix_temp/A5b_p.res$par[2])^2)
A5b_mat <- as.matrix((A5b_p.res$par[3]*(1-A5b_p.res$par[1])*spmatrix_temp+(1-A5b_p.res$par[1])*(1-A5b_p.res$par[3])*A5b_phylomatrix+A5b_p.res$par[1]*diag(dim(A5b_phylomatrix)[1])))
save(A5b_mat, file = "output/A5b/A5b_mat.RData")

# --- Maximum likelihood fits for all models ---
A5b_fits <- mclapply(
  X = A5b_models,
  FUN = function(f) ml_fit(p=A5b_p.res$par,formula=f,data=A5b,spmatrix=A5b_spmatrix,phylomatrix=A5b_phylomatrix),
  mc.cores = 8
)
# Save models
saveRDS(A5b_fits, "output/A5b/A5b_fits.RDS") 

# --- Model Comparison ---
# Make empty NA list of coefficients and their respective p values in tuple form
n <- length(A5b_fits)
A5b_summaries <- list("Island.Endemic" = as.list(rep(NA,n)),
                     "Range.Size..km2." = as.list(rep(NA,n)),
                     "altitude_range" = as.list(rep(NA,n)), 
                     "bordering_language_richness" = as.list(rep(NA,n)),
                     "L1_pop" = as.list(rep(NA,n)))
ml_logLik <- c()
ml_BIC <- c()
ml_AIC <- c()
ml_MSE <- c()

for (i in 1:n){
  # get current model summary
  ml_summary <- summary(A5b_fits[[i]])
  # if intercept only model just get BIC, logLik and MLE
  if(dim(ml_summary$tTable)[1] == 1) {
    ml_logLik <- c(ml_logLik, ml_summary$logLik)
    ml_BIC <- c(ml_BIC, ml_summary$BIC)
    ml_AIC <- c(ml_AIC, ml_summary$AIC)
    ml_MSE <- c(ml_MSE, sum(A5b_fits[[i]]$residuals**2)/dim(A5b)[1])
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
  # for (j in 2:n_coef){
  #   if(coef_p[j]<=0.05){
  #     A5b_summaries[[coef_names[j]]][[i]] <- paste0(round(coef_val[j], 5),"***")
  #   } else {
  #     A5b_summaries[[coef_names[j]]][[i]] <- paste0(round(coef_val[j], 5))
  #   }
  # }
  for (j in 2:n_coef){
    A5b_summaries[[coef_names[j]]][[i]] <- paste0(round(coef_val[j],4),",",round(coef_p[j],4))
  }
  # Add BIC, logLik and MLE
  ml_logLik <- c(ml_logLik, ml_summary$logLik)
  ml_BIC <- c(ml_BIC, ml_summary$BIC)
  ml_AIC <- c(ml_AIC, ml_summary$AIC)
  ml_MSE <- c(ml_MSE, sum(A5b_fits[[i]]$residuals**2)/dim(A5b)[1])
}

# Create new data frame from lists
A5b_ml_data <- data.frame(BIC = ml_BIC,
                         AIC = ml_AIC,
                         logLik = ml_logLik,
                         MSE = ml_MSE,
                         unlist(A5b_summaries$Island.Endemic),
                         unlist(A5b_summaries$Range.Size..km2.),
                         unlist(A5b_summaries$altitude_range),
                         unlist(A5b_summaries$bordering_language_richness),
                         unlist(A5b_summaries$L1_pop),
                         stringsAsFactors = FALSE)
colnames(A5b_ml_data) <- c("BIC", "AIC", "logLik", "MSE", A5b_predictors)
# Sort based on increasing BIC
A5b_ml_data_BIC <- A5b_ml_data[order(A5b_ml_data$BIC), ]
# Add delta BIC column (min_BIC - current_BIC)
min_BIC <-  min(A5b_ml_data_BIC$BIC)
A5b_ml_data_BIC$deltaBIC <- A5b_ml_data_BIC$BIC - min_BIC
# Save output
write.csv(A5b_ml_data_BIC, file = "output/A5b/A5b_ml_data_BIC.csv", row.names=FALSE)
# Sort based on increasing AIC
A5b_ml_data_AIC <- A5b_ml_data[order(A5b_ml_data$AIC), ]
# Add delta AIC column (min_AIC - current_AIC)
min_AIC <-  min(A5b_ml_data_AIC$AIC)
A5b_ml_data_AIC$deltaAIC <- A5b_ml_data_AIC$AIC - min_AIC
# Save output
write.csv(A5b_ml_data_AIC, file = "output/A5b/A5b_ml_data_AIC.csv", row.names=FALSE)

## --- A5c --- ####

# --- Get PISc ~ Area OLS ---
A5c_PIS_Range_l1 = lm(data=A5c, Phoneme.Inventory.Size ~ Range.Size..km2.)
summary(A5c_PIS_Range_l1)
# For each region
A5c_PIS_Range_Region_lm <- lapply(region_order, function(x) {
  summary(lm(data=A5c, A5c$Phoneme.Inventory.Size[A5c$region==x] ~ A5c$Range.Size..km2.[A5c$region==x]))
})
# For each region with only non-island languages
A5c_PIS_NI_Range_Region_lm <- lapply(region_order, function(x) {
  summary(lm(data=A5c,
             A5c$Phoneme.Inventory.Size[A5c$region==x & A5c$Island.Endemic == 0] ~ A5c$Range.Size..km2.[A5c$region==x  & A5c$Island.Endemic == 0]))
})
# For only Oceania region island languages
A5c_PIS_Range_l2 = lm(data=A5c[A5c$region == "Oceania" & A5c$Island.Endemic == 1,], Phoneme.Inventory.Size ~ Range.Size..km2.)

# --- Model predictor combinations ---
A5c_predictors <- c("Island.Endemic" ,"Range.Size..km2.",
                    "Distance.to.Mainland", "Distance.to.Continent",
                    "L1_pop")
response <- "Phoneme.Inventory.Size"
n_pred <- length(A5c_predictors)
pred_combs <- sapply(1:n_pred, function(x) combn(A5c_predictors, x))
A5c_models <- c(Phoneme.Inventory.Size~1)
for (i in 1:n_pred) {
  #print(pred_combs[[i]])
  for (j in 1:dim(pred_combs[[i]])[2]) {
    model <- (paste(response, paste(pred_combs[[i]][,j], collapse="+"), sep="~"))
    model <- as.formula(model)
    A5c_models <- c(A5c_models,model)
  }
}

# --- Best p values for all models ---
A5c_p.res <- sbplx(c(0.5, 0.5, 0.5),
                   best_p,
                   formula=Phoneme.Inventory.Size~1,
                   data=A5c,
                   spmatrix=A5c_spmatrix,phylomatrix=A5c_phylomatrix,
                   lower=c(0,0,0),upper=c(1,1,1),
                   nl.info = TRUE)
# Save p.res and covariance matrix
save(A5c_p.res, file = "output/A5c/A5c_p_res.RData")
spmatrix_temp <- A5c_spmatrix/max(A5c_spmatrix)
spmatrix_temp <- exp(-(spmatrix_temp/A5c_p.res$par[2])^2)
A5c_mat <- as.matrix((A5c_p.res$par[3]*(1-A5c_p.res$par[1])*spmatrix_temp+(1-A5c_p.res$par[1])*(1-A5c_p.res$par[3])*A5c_phylomatrix+A5c_p.res$par[1]*diag(dim(A5c_phylomatrix)[1])))
save(A5c_mat, file = "output/A5c/A5c_mat.RData")

# --- Maximum likelihood fits for all models ---
A5c_fits <- mclapply(
  X = A5c_models,
  FUN = function(f) ml_fit(p=A5c_p.res$par,formula=f,data=A5c,spmatrix=A5c_spmatrix,phylomatrix=A5c_phylomatrix),
  mc.cores = 8
)
# Save models
saveRDS(A5c_fits, "output/A5c/A5c_fits.RDS") 

# --- Model Comparison ---
# Make empty NA list of coefficients and their respective p values in tuple form
n <- length(A5c_fits)
A5c_summaries <- list("Island.Endemic" = as.list(rep(NA,n)),
                      "Range.Size..km2." = as.list(rep(NA,n)),
                      "Distance.to.Mainland" = as.list(rep(NA,n)), 
                      "Distance.to.Continent" = as.list(rep(NA,n)),
                      "L1_pop" = as.list(rep(NA,n)))
ml_logLik <- c()
ml_BIC <- c()
ml_AIC <- c()
ml_MSE <- c()

for (i in 1:n){
  # get current model summary
  ml_summary <- summary(A5c_fits[[i]])
  # if intercept only model just get BIC, logLik and MLE
  if(dim(ml_summary$tTable)[1] == 1) {
    ml_logLik <- c(ml_logLik, ml_summary$logLik)
    ml_BIC <- c(ml_BIC, ml_summary$BIC)
    ml_AIC <- c(ml_AIC, ml_summary$AIC)
    ml_MSE <- c(ml_MSE, sum(A5c_fits[[i]]$residuals**2)/dim(A5c)[1])
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
  # for (j in 2:n_coef){
  #   if(coef_p[j]<=0.05){
  #     A5c_summaries[[coef_names[j]]][[i]] <- paste0(round(coef_val[j], 5),"***")
  #   } else {
  #     A5c_summaries[[coef_names[j]]][[i]] <- paste0(round(coef_val[j], 5))
  #   }
  # }
  for (j in 2:n_coef){
    A5c_summaries[[coef_names[j]]][[i]] <- paste0(round(coef_val[j],4),",",round(coef_p[j],4))
  }
  # Add BIC, logLik and MLE
  ml_logLik <- c(ml_logLik, ml_summary$logLik)
  ml_BIC <- c(ml_BIC, ml_summary$BIC)
  ml_AIC <- c(ml_AIC, ml_summary$AIC)
  ml_MSE <- c(ml_MSE, sum(A5c_fits[[i]]$residuals**2)/dim(A5c)[1])
}

# Create new data frame from lists
A5c_ml_data <- data.frame(BIC = ml_BIC,
                          AIC = ml_AIC,
                          logLik = ml_logLik,
                          MSE = ml_MSE,
                          unlist(A5c_summaries$Island.Endemic),
                          unlist(A5c_summaries$Range.Size..km2.),
                          unlist(A5c_summaries$Distance.to.Mainland),
                          unlist(A5c_summaries$Distance.to.Continent),
                          unlist(A5c_summaries$L1_pop),
                          stringsAsFactors = FALSE)
colnames(A5c_ml_data) <- c("BIC", "AIC", "logLik", "MSE", A5c_predictors)
# Sort based on increasing BIC
A5c_ml_data_BIC <- A5c_ml_data[order(A5c_ml_data$BIC), ]
# Add delta BIC column (min_BIC - current_BIC)
min_BIC <-  min(A5c_ml_data_BIC$BIC)
A5c_ml_data_BIC$deltaBIC <- A5c_ml_data_BIC$BIC - min_BIC
# Save output
write.csv(A5c_ml_data_BIC, file = "output/A5c/A5c_ml_data_BIC.csv", row.names=FALSE)
# Sort based on increasing AIC
A5c_ml_data_AIC <- A5c_ml_data[order(A5c_ml_data$AIC), ]
# Add delta AIC column (min_AIC - current_AIC)
min_AIC <-  min(A5c_ml_data$AIC)
A5c_ml_data_AIC$deltaAIC <- A5c_ml_data_AIC$AIC - min_AIC
# Save output
write.csv(A5c_ml_data_AIC, file = "output/A5c/A5c_ml_data_AIC.csv", row.names=FALSE)

## --- A5d --- ####

# --- Model predictor combinations ---
A5d_model <- c(Phoneme.Inventory.Size ~ documentation)

# --- Best p values for all models ---
A5d_p.res <- sbplx(c(0.5, 0.5, 0.5),
                  best_p,
                  formula=A5d_model[[1]],
                  data=A5d,
                  spmatrix=A5d_spmatrix,phylomatrix=A5d_phylomatrix,
                  lower=c(0,0,0),upper=c(1,1,1),
                  nl.info = TRUE)
# Save p.res and covariance matrix
save(A5d_p.res, file = "output/A5d/A5d_p_res.RData")
spmatrix_temp <- A5d_spmatrix/max(A5d_spmatrix)
spmatrix_temp <- exp(-(spmatrix_temp/A5d_p.res$par[2])^2)
mat <- as.matrix((A5d_p.res$par[3]*(1-A5d_p.res$par[1])*spmatrix_temp+(1-A5d_p.res$par[1])*(1-A5d_p.res$par[3])*A5d_phylomatrix+A5d_p.res$par[1]*diag(dim(A5d_phylomatrix)[1])))

# --- Maximum likelihood fits for all models A1b ---
A5d_model <- ml_fit(p=A5d_p.res$par,formula=A5d_model[[1]],data=A5d,spmatrix=A5d_spmatrix,phylomatrix=A5d_phylomatrix)
# Save model fit
save(A5d_model, file = "output/A5d/A5d_model.RData") 

## --- A6 --- ####

# --- Model predictor combinations ---
A6_model <- c(L1_pop ~ documentation)

# --- OLS ---
A6_OLS <- lm(data=A6, formula = L1_pop ~ documentation)
summary(A6_OLS)

# --- Best p values for GLS models ---
A6_p.res <- sbplx(c(0.5, 0.5, 0.5),
                  best_p,
                  formula=L1_pop ~ documentation,
                  data=A6,
                  spmatrix=A6_spmatrix,phylomatrix=A6_phylomatrix,
                  lower=c(0,0,0),upper=c(1,1,1),
                  nl.info = TRUE)
# Save p.res and covariance matrix
save(A6_p.res, file = "output/A6/A6_p_res.RData")
spmatrix_temp <- A6_spmatrix/max(A6_spmatrix)
spmatrix_temp <- exp(-(spmatrix_temp/A6_p.res$par[2])^2)
mat <- as.matrix((A6_p.res$par[3]*(1-A6_p.res$par[1])*spmatrix_temp+(1-A6_p.res$par[1])*(1-A6_p.res$par[3])*A6_phylomatrix+A6_p.res$par[1]*diag(dim(A6_phylomatrix)[1])))

# --- Maximum likelihood fits for all models A1b ---
A6_model <- ml_fit(p=A6_p.res$par,formula=L1_pop ~ documentation,data=A6,spmatrix=A6_spmatrix,phylomatrix=A6_phylomatrix)
# Save model fit
save(A6_model, file = "output/A6/A6_model.RData") 
summary(A6_model)

# --- (8) RESULTS PLOTS --- ####

## --- A1a --- ####

# Show regression lines (with and without correction for autocorrelation)
PISa_L1_pop_wr = ggplot(A1a, aes(x = L1_pop, y = Phoneme.Inventory.Size, color = factor(region, levels = region_order))) + 
  geom_point() + 
  geom_abline(intercept = A1a_OLS$coefficients[1], slope = A1a_OLS$coefficients[2],
              color = "red", linewidth = 1) +
  geom_abline(intercept = A1a_GLS_model$coefficients[1], slope = A1a_GLS_model$coefficients[2],
              color = "blue", linewidth = 1) +
  labs(x = "L1_pop", y = "PISa", color = "Region") + 
  scale_color_manual(values = region_colours) +
  theme_minimal()
# Save graph
ggsave("output/A1a/PISa_L1_pop_wr.png", plot = PISa_L1_pop_wr, width = 15, height = 10, dpi = 100,scale = 0.5)
# Show graph
print(PISa_L1_pop_wr)

## --- A3a --- ####

# --- PISc ~ L1_pop ---
# Show regression lines (with and without correction for autocorrelation)
PISc_L1_pop_wr = ggplot(A3a, aes(x = L1_pop, y = Phoneme.Inventory.Size, color = factor(region, levels = region_order))) + 
  geom_point() + 
  geom_abline(intercept = A3a_l1$coefficients[1], slope = A3a_l1$coefficients[2],
              color = "red", linewidth = 1) +
  geom_abline(intercept = A3a_fits[[6]]$coefficients[1], slope = A3a_fits[[6]]$coefficients[2],
              color = "blue", linewidth = 1) +
  labs(x = "L1_pop", y = "PISc", color = "Region") + 
  scale_color_manual(values = region_colours) +
  theme_minimal()
# Save graph
ggsave("output/A3a/PISc_L1_pop_wr.png", plot = PISc_L1_pop_wr, width = 15, height = 10, dpi = 100,scale = 0.5)
# Show graph
print(PISc_L1_pop_wr)

# --- PISc ~ Area ---
# Show regression lines (with and without correction for autocorrelation)
PISc_Area_wr = ggplot(A3a, aes(x = Range.Size..km2., y = Phoneme.Inventory.Size, color = factor(region, levels = region_order))) + 
  geom_point() + 
  geom_abline(intercept = A3a_PIS_Range_l1$coefficients[1], slope = A3a_PIS_Range_l1$coefficients[2],
              color = "red", linewidth = 1) +
  geom_abline(intercept = A3a_fits[[3]]$coefficients[1], slope = A3a_fits[[3]]$coefficients[2],
              color = "blue", linewidth = 1) +
  labs(x = "Area", y = "PISc", color = "Region") + 
  scale_color_manual(values = region_colours) +
  theme_minimal()
# Save graph
ggsave("output/A3a/PISc_Area_wr.png", plot = PISc_Area_wr, width = 15, height = 10, dpi = 100,scale = 0.5)
# Show graph
print(PISc_Area_wr)


# --- (9) REFERENCES --- ####
#
# (1) Bromham, L., Yaxley, K.J. & Cardillo, M. Islands are engines of language diversity. Nat Ecol Evol 8, 1991–2002 (2024). https://doi.org/10.1038/s41559-024-02488-4
# 
# (2) Hua, X., Greenhill, S.J., Cardillo, M. et al. The ecological drivers of variation in global language diversity. Nat Commun 10, 2047 (2019). https://doi.org/10.1038/s41467-019-09842-2
# 
# (3) Bromham, L., Dinnage, R., Skirgård, H. et al. Global predictors of language endangerment and the future of linguistic diversity. Nat Ecol Evol 6, 163–173 (2022). https://doi.org/10.1038/s41559-021-01604-y
# 
# (4) Cormac Anderson, Tiago Tresoldi, Simon J Greenhill, Robert Forkel, Russell Gray, Johann-Mattis List, Variation in phoneme inventories: quantifying the problem and improving comparability, Journal of Language Evolution, Volume 8, Issue 2, July 2023, Pages 149–168, https://doi.org/10.1093/jole/lzad011
# 