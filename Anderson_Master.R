library(nloptr)
library(nlme)
library(parallel)

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
write.csv(data, file = "output/anderson/anderson_adjusted.csv", row.names=FALSE)

### Preview and save matrices
print(dim(phylomatrix))
print(dim(spmatrix))
print(phylomatrix[1:10,1:10])
print(spmatrix[1:10,1:10])
save(phylomatrix, file = "output/anderson/phylomatrix.RData") 
save(spmatrix, file = "output/anderson/spmatrix.RData")
save(Wnb, file = "output/anderson/Wnb.Rdata")

# --------
# ANALYSIS
# --------

# Get raw anderson adjusted data and matrices
raw_anderson <- read.csv("data/anderson_adjusted.csv")
load("data/anderson/phylomatrix.RData")
load("data/anderson/spmatrix.RData")

# --- Log transform variables ---
anderson <- raw_anderson
anderson$Phoneme.Inventory.Size <- log(anderson$Phoneme.Inventory.Size) 
anderson$Range.Size..km2. <- log(anderson$Range.Size..km2.)
anderson$L1_pop <- log(anderson$L1_pop + 0.5)
# Do not transform zero entries for distance variables
anderson$Distance.to.Mainland[anderson$Distance.to.Mainland != 0] <- log(
  anderson$Distance.to.Mainland[anderson$Distance.to.Mainland != 0])
anderson$Distance.to.Continent[anderson$Distance.to.Continent != 0] <- log(
  anderson$Distance.to.Continent[anderson$Distance.to.Continent != 0])

# Preview transformed data
head(anderson)
summary(anderson)

# --- Preview correlations ---
pairs(anderson[,
               c("Phoneme.Inventory.Size", "Range.Size..km2.",
                 "L1_pop", "Distance.to.Mainland", "Distance.to.Continent",
                 "altitude_range", "bordering_language_richness")],
      panel = function(x, y, ...) {
        points(x, y, cex = 0.1, ...) # cex = 0.5 makes points half the default size
      })
cor(anderson[,
             c("Phoneme.Inventory.Size", "Range.Size..km2.",
               "L1_pop", "Distance.to.Mainland", "Distance.to.Continent",
               "altitude_range", "bordering_language_richness")])
cor.test(anderson$Phoneme.Inventory.Size,anderson$Island.Endemic)

# --- Testing GLS (2) ---
p = c(0.5,0.5,0.5)
spmatrix_test <- spmatrix/max(spmatrix)
spmatrix_test <- exp(-(spmatrix_test/p[2])^2)
mat <- as.matrix((p[3]*(1-p[1])*spmatrix_test+(1-p[1])*(1-p[3])*phylomatrix+p[1]*diag(dim(phylomatrix)[1])))
res <- gls(model=Phoneme.Inventory.Size~1,data=anderson,correlation=corSymm(mat[lower.tri(mat)],fixed=T),method="ML")
print(-res$logLik)

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
predictors <- c("Island.Endemic", "Range.Size..km2.",
                "Distance.to.Mainland", "Distance.to.Continent",
                "L1_pop", "altitude_range", 
                "bordering_language_richness")
response <- "Phoneme.Inventory.Size"
n_pred <- length(predictors)
pred_combs <- sapply(1:n_pred, function(x) combn(predictors, x))
models <- c()
for (i in 1:n_pred) {
  #print(pred_combs[[i]])
  for (j in 1:dim(pred_combs[[i]])[2]) {
    model <- (paste(response, paste(pred_combs[[i]][,j], collapse="+"), sep="~"))
    model <- as.formula(model)
    models <- c(models,model)
  }
}
unlist(models)

# --- Best p values for all models ---
p.res <- sbplx(c(0.5, 0.5, 0.5),
               best_p,
               formula=Phoneme.Inventory.Size~1,
               data=anderson,
               spmatrix=spmatrix,phylomatrix=phylomatrix,
               lower=c(0,0,0),upper=c(1,1,1),
               nl.info = TRUE)

# --- Maximum likelihood fits for all models ---
ml_fits <- mclapply(
  X = models,
  FUN = function(f) ml_fit(p=p.res$par,formula=f,data=anderson,spmatrix=spmatrix,phylomatrix=phylomatrix),
  mc.cores = 8
)
# Save models
saveRDS(ml_fits, "output/anderson/ml_fits.RData") 

# --- Model Comparison ---
# Make empty NA list of coefficients and their respective p values in tuple form
n <- length(ml_fits)
ml_summaries <- list("Island.Endemic" = as.list(rep(NA,n)),
                     "Range.Size..km2." = as.list(rep(NA,n)),
                     "Distance.to.Mainland" = as.list(rep(NA,n)), 
                     "Distance.to.Continent" = as.list(rep(NA,n)),
                     "L1_pop" = as.list(rep(NA,n)), 
                     "altitude_range" = as.list(rep(NA,n)), 
                     "bordering_language_richness" = as.list(rep(NA,n)))
ml_logLik <- c()
ml_BIC <- c()
ml_MSE <- c()

for (i in 1:n){
  # get current model summary
  ml_summary <- summary(ml_fits[[i]])
  # if intercept only model just get BIC, logLik and MLE
  if(dim(ml_summary$tTable)[1] == 1) {
    ml_logLik <- c(ml_logLik, ml_summary$logLik)
    ml_BIC <- c(ml_BIC, ml_summary$BIC)
    ml_MSE <- c(ml_MSE, sum(ml_fits[[i]]$residuals**2)/dim(anderson)[1])
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
      ml_summaries[[coef_names[j]]][[i]] <- paste0(round(coef_val[j], 5),"***")
    } else {
      ml_summaries[[coef_names[j]]][[i]] <- paste0(round(coef_val[j], 5))
    }
  }
  # Add BIC, logLik and MLE
  ml_logLik <- c(ml_logLik, ml_summary$logLik)
  ml_BIC <- c(ml_BIC, ml_summary$BIC)
  ml_MSE <- c(ml_MSE, sum(ml_fits[[i]]$residuals**2)/dim(anderson)[1])
}

# Create new data frame from lists
ml_data <- data.frame(BIC = ml_BIC,
                      logLik = ml_logLik,
                      MSE = ml_MSE,
                      unlist(ml_summaries$Island.Endemic),
                      unlist(ml_summaries$Range.Size..km2.),
                      unlist(ml_summaries$Distance.to.Mainland),
                      unlist(ml_summaries$Distance.to.Continent),
                      unlist(ml_summaries$L1_pop),
                      unlist(ml_summaries$altitude_range),
                      unlist(ml_summaries$bordering_language_richness),
                      stringsAsFactors = FALSE)
colnames(ml_data) <- c("BIC", "logLik", "MSE", predictors)
# Sort based on increasing BIC
ml_data <- ml_data[order(ml_data$BIC), ]
# Add delta BIC column (min_BIC - current_BIC)
min_BIC <-  min(ml_data$BIC)
ml_data$deltaBIC <- ml_data$BIC - min_BIC
# Save output
write.csv(ml_data, file = "output/anderson/ml_data.csv", row.names=FALSE)

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

