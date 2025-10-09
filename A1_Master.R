library(nloptr)
library(nlme)
library(parallel)
library(ggplot2)

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

### Preview Matrices
print(phylomatrix[1:5,1:5])
print(spmatrix[1:5,1:5])

# --- Loading NEE22 (3) ---
nee22 <- read.csv("data/Global predictors of language endangerment and the future of linguistic diversity Data 2.csv",
                  stringsAsFactors = FALSE)

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
print(unique(phoible$Glottocode[glotto_duplicates]))
phoible <- phoible[!unlist(glotto_duplicates),]
print("Size of dataset after duplicates removed:")
print(dim(phoible)[1])

# rename Sounds column to Phoneme.Inventory.Size
colnames(phoible)[colnames(phoible) == "Sounds"] <- "Phoneme.Inventory.Size"

# --- Preview data ---
head(phoible)
summary(phoible)

## Merging Datasets
data <- merge(
  nee22[, 
        c("ISO", "L1_pop", "region")],
  phoible[, c("ISO639P3code", "Glottocode", "Latitude", "Longitude", "Phoneme.Inventory.Size")], by.x = "ISO", by.y = "ISO639P3code", all.x = TRUE)

# --- Data Cleanup ---
# Remove all NA entries
print("Size of initial dataset:")
print(dim(data)[1])
old_data <- data
data <- na.omit(data) 
print("Size of dataset with NA removed:")
print(dim(data)[1])
# Remove iso duplicate entries
print("Duplicate iso entries:")
iso_duplicates <- which(table(data$ISO) != 1)
print(which(table(data$ISO) != 1))
data <- data[!(data$ISO %in% names(which(table(data$ISO) != 1))),]
print("Size of dataset after duplicates removed:")
print(dim(data[1]))

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
write.csv(data, file = "output/A1/A1_adjusted.csv", row.names=FALSE)

# --- Preview and save matrices ---
print(dim(phylomatrix))
print(dim(spmatrix))
print(phylomatrix[1:10,1:10])
print(spmatrix[1:10,1:10])
save(phylomatrix, file = "output/A1/phylomatrix.RData") 
save(spmatrix, file = "output/A1/spmatrix.RData")

# --------
# ANALYSIS
# --------

# --- Get raw A1 adjusted data and matrices ---
raw_A1 <- read.csv("output/A1/A1_adjusted.csv")
load("output/A1/phylomatrix.RData")
load("output/A1/spmatrix.RData")

# --- Preview raw data ---
head(raw_A1)
summary(raw_A1)
# Phoneme Inventory Size distribution
PIS_hist <- ggplot(raw_A1, aes(x = Phoneme.Inventory.Size)) +
  geom_histogram() +
  labs(
    y = "Number of Languages Observed",
    x = "Phoneme.Inventory.Size",
  ) + 
  scale_x_continuous(n.breaks = 15) +
  scale_y_continuous(n.breaks = 15)
  
ggsave(
  filename = "output/A1/PIS_hist.png",
  plot = PIS_hist,
  scale = 1
)
print(PIS_hist)

# --- Log transform variables ---
A1 <- raw_A1
A1$Phoneme.Inventory.Size <- log(A1$Phoneme.Inventory.Size) 
A1$L1_pop <- log(A1$L1_pop + 0.5)

# --- Preview transformed data ---
head(A1)
summary(A1)

# --- Preview correlations ---
plot(Phoneme.Inventory.Size ~ L1_pop,data=A1)

# --- Ordinary Least Squares ---
l1 <- lm(formula = Phoneme.Inventory.Size ~ L1_pop, data = A1)
print(summary(l1))
print(BIC(l1))
print(logLik(l1))

# --- Testing GLS (2) ---
# p = c(0.5,0.5,0.5)
# spmatrix_test <- spmatrix/max(spmatrix)
# spmatrix_test <- exp(-(spmatrix_test/p[2])^2)
# mat <- as.matrix((p[3]*(1-p[1])*spmatrix_test+(1-p[1])*(1-p[3])*phylomatrix+p[1]*diag(dim(phylomatrix)[1])))
# res <- gls(model=Phoneme.Inventory.Size~1,data=A1,correlation=corSymm(mat[lower.tri(mat)],fixed=T),method="ML")
# print(-res$logLik)

# --- GLS (2) Setup ---
best_p <- function (p,formula,data,spmatrix,phylomatrix) {
  spmatrix <- spmatrix/max(spmatrix)
  spmatrix <- exp(-(spmatrix/p[2])^2)
  mat <- as.matrix((p[3]*(1-p[1])*spmatrix+(1-p[1])*(1-p[3])*phylomatrix+p[1]*diag(dim(phylomatrix)[1])))
  res <- try(gls(model=formula,data=data,correlation=corSymm(mat[lower.tri(mat)],fixed=T),method="ML"),silent=T)
  if (inherits(res,"try-error")) {
    out <- -10000
  } else {
    out <- res$logLik
    if (res$logLik>0) {out <- -10000}
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

# --- Model predictor combinations ---
predictors <- c("L1_pop")
models <- c(Phoneme.Inventory.Size ~ L1_pop)

# --- Best p values for all models ---
p.res <- sbplx(c(0.5, 0.5, 0.5),
               best_p,
               formula=models[[1]],
               data=A1,
               spmatrix=spmatrix,phylomatrix=phylomatrix,
               lower=c(0,0,0),upper=c(1,1,1),
               nl.info = TRUE)
# Save p.res
save(p.res, file = "output/A1/ml_p_res.RData")

# --- Maximum likelihood fits for all models ---
ml_model <- ml_fit(p=p.res$par,formula=models[[1]],data=A1,spmatrix=spmatrix,phylomatrix=phylomatrix)
# Save model fit
save(ml_model, file = "output/A1/ml_model.RData") 

# --- References ---
#
# (1) Bromham, L., Yaxley, K.J. & Cardillo, M. Islands are engines of language diversity. Nat Ecol Evol 8, 1991–2002 (2024). https://doi.org/10.1038/s41559-024-02488-4
# 
# (2) Hua, X., Greenhill, S.J., Cardillo, M. et al. The ecological drivers of variation in global language diversity. Nat Commun 10, 2047 (2019). https://doi.org/10.1038/s41467-019-09842-2
# 
# (3) Bromham, L., Dinnage, R., Skirgård, H. et al. Global predictors of language endangerment and the future of linguistic diversity. Nat Ecol Evol 6, 163–173 (2022). https://doi.org/10.1038/s41559-021-01604-y
# 
# (4) Cormac Anderson, Tiago Tresoldi, Simon J Greenhill, Robert Forkel, Russell Gray, Johann-Mattis List, Variation in phoneme inventories: quantifying the problem and improving comparability, Journal of Language Evolution, Volume 8, Issue 2, July 2023, Pages 149–168, https://doi.org/10.1093/jole/lzad011
# 
