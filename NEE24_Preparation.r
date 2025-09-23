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

# --- Loading NEE24 (1) ---
nee24 <- read.csv("data/Islands are engines of language diversity data.csv",
                  stringsAsFactors = FALSE)

# Remove all NA entries
paste("Size of initial dataset:", dim(nee24)[1])
paste("#Island_Endemic: ", sum(nee24$Island.Endemic))
paste("#Island_Endemic to Total Ratio: ", sum(nee24$Island.Endemic) / dim(nee24)[1])
old_nee24 <- nee24
nee24 <- na.omit(nee24) 
paste("Number of incomplete (NA) language entries removed:", 
      dim(old_nee24)[1] - dim(nee24)[1])
paste("Size of dataset with NA remove:", dim(nee24)[1])
paste("#Island_Endemic: ", sum(nee24$Island.Endemic))
paste("#Island_Endemic to Total Ratio: ", sum(nee24$Island.Endemic) / dim(nee24)[1])

# Convert variable types
nee24$Island.Endemic <- as.numeric(nee24$Island.Endemic)

# --- Adjusting NEE24 data with matrices ---
common_ids <- Reduce(intersect, list(
  nee24$ISO,
  rownames(phylomatrix),
  rownames(spmatrix),
  rownames(Wnb)
))
nee24 <- nee24[nee24$ISO %in% common_ids, ]
phylomatrix <- phylomatrix[common_ids, common_ids]
spmatrix <- spmatrix[common_ids, common_ids]
Wnb <- Wnb[common_ids, common_ids]

# Preview and save data
paste("Size of dataset adjusted:", dim(nee24)[1])
paste("#Island_Endemic: ", sum(nee24$Island.Endemic))
paste("#Island_Endemic to Total Ratio: ", sum(nee24$Island.Endemic) / dim(nee24)[1])
# save data 
write.csv(nee24, file = "data/nee24_adjusted.csv", row.names=FALSE)

# Preview and save matrices
print(dim(phylomatrix))
print(dim(spmatrix))
print(dim(Wnb))
print(phylomatrix[1:10,1:10])
print(spmatrix[1:10,1:10])
print(Wnb[1:10,1:10])
save(phylomatrix, file = "data/phylomatrix.RData") 
save(spmatrix, file = "data/spmatrix.RData")
save(Wnb, file = "data/Wnb.Rdata")