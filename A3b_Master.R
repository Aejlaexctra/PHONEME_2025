library(nloptr)
library(nlme)
library(parallel)

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

# --- Loading NEE24 (1) ---
nee24 <- read.csv("data/Islands are engines of language diversity data.csv",
                  stringsAsFactors = FALSE)

# --- Loading NEE22 (3) ---
nee22 <- read.csv("data/Global predictors of language endangerment and the future of linguistic diversity Data 2.csv",
                  stringsAsFactors = FALSE)

# --- Preview NEE24 ---
paste("Size of initial dataset:", dim(nee24)[1])
paste("#Island_Endemic: ", sum(nee24$Island.Endemic))
paste("#Island_Endemic to Total Ratio: ", sum(nee24$Island.Endemic) / dim(nee24)[1])

# --- Prepare dataset ---
# Merge datasets
A3b <- merge(nee22[, c("ISO", "region", "bordering_language_richness")], 
            nee24[, c("ISO693.3", "L1.Population","Phoneme.Inventory.Size", "Island.Endemic", 
                      "Distance.to.Mainland", "Distance.to.Continent", "Range.Size..km2.")], 
              by.x = "ISO", by.y = "ISO693.3", all.x = TRUE)
colnames(A3b)[colnames(A3b) == "L1.Population"] <- "L1_pop"
# Remove all NA entries
A3b <- na.omit(A3b) 
paste("Size of dataset with NA remove:", dim(A3b)[1])
paste("#Island_Endemic: ", sum(A3b$Island.Endemic))
paste("#Island_Endemic to Total Ratio: ", sum(A3b$Island.Endemic) / dim(A3b)[1])
# Convert variable types
A3b$Island.Endemic <- as.numeric(A3b$Island.Endemic)
# Rename predictors
colnames(A3b)[colnames(A3b) == "L1.Population"] <- "L1_pop"
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
phylomatrix <- phylomatrix[common_ids, common_ids]
spmatrix <- spmatrix[common_ids, common_ids]

# Preview and save data
paste("Size of dataset adjusted:", dim(A3b)[1])
paste("#Island_Endemic: ", sum(A3b$Island.Endemic))
paste("#Island_Endemic to Total Ratio: ", sum(A3b$Island.Endemic) / dim(A3b)[1])
# save data 
write.csv(A3b, file = "output/A3b/A3b_adjusted.csv", row.names=FALSE)

# Preview and save matrices
print(dim(phylomatrix))
print(dim(spmatrix))
print(phylomatrix[1:10,1:10])
print(spmatrix[1:10,1:10])
save(phylomatrix, file = "output/A3b/phylomatrix.RData") 
save(spmatrix, file = "output/A3b/spmatrix.RData")

# --------
# ANALYSIS
# --------

# Get raw A3b adjusted data and matrices
raw_A3b <- read.csv("output/A3b/A3b_adjusted.csv")
load("output/A3b/phylomatrix.RData")
load("output/A3b/spmatrix.RData")

# --- Log transform variables ---
A3b <- raw_A3b
A3b$Phoneme.Inventory.Size <- log(A3b$Phoneme.Inventory.Size) 
A3b$Range.Size..km2. <- log(A3b$Range.Size..km2.)
A3b$L1_pop <- log(A3b$L1_pop + 0.5)
A3b$bordering_language_richness <- log(A3b$bordering_language_richness + 0.5)
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

# --- Testing GLS (2) ---
# p = c(0.5,0.5,0.5)
# spmatrix_test <- spmatrix/max(spmatrix)
# spmatrix_test <- exp(-(spmatrix_test/p[2])^2)
# mat <- as.matrix((p[3]*(1-p[1])*spmatrix_test+(1-p[1])*(1-p[3])*phylomatrix+p[1]*diag(dim(phylomatrix)[1])))
# res <- gls(model=Phoneme.Inventory.Size~1,data=A3b,correlation=corSymm(mat[lower.tri(mat)],fixed=T),method="ML")
# print(-res$logLik)

# --- Preview correlations ---
pairs(A3b[,
         c("Range.Size..km2.",
           "Distance.to.Mainland", "Distance.to.Continent",
           "L1_pop", "bordering_language_richness")],
      panel = function(x, y, ...) {
        points(x, y, cex = 0.1, ...) # cex = 0.5 makes points half the default size
      })
cor(A3b[,
       c("Range.Size..km2.",
         "Distance.to.Mainland", "Distance.to.Continent",
         "L1_pop", "bordering_language_richness")])

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

# --- Model predictor combinations ---
predictors <- c("Range.Size..km2.",
                "Distance.to.Mainland", "Distance.to.Continent",
                "L1_pop", "bordering_language_richness")
response <- "Phoneme.Inventory.Size"
n_pred <- length(predictors)
pred_combs <- sapply(1:n_pred, function(x) combn(predictors, x))
models <- c(Phoneme.Inventory.Size~1)
for (i in 1:n_pred) {
  #print(pred_combs[[i]])
  for (j in 1:dim(pred_combs[[i]])[2]) {
    model <- (paste(response, paste(pred_combs[[i]][,j], collapse="+"), sep="~"))
    model <- as.formula(model)
    models <- c(models,model)
  }
}

# --- Best p values for all models ---
p.res <- sbplx(c(0.5, 0.5, 0.5),
               best_p,
               formula=Phoneme.Inventory.Size~1,
               data=A3b,
               spmatrix=spmatrix,phylomatrix=phylomatrix,
               lower=c(0,0,0),upper=c(1,1,1),
               nl.info = TRUE)
# Save p.res and covariance matrix
save(p.res, file = "output/A3b/ml_p_res.RData")
spmatrix_temp <- spmatrix/max(spmatrix)
spmatrix_temp <- exp(-(spmatrix_temp/p.res$par[2])^2)
mat <- as.matrix((p.res$par[3]*(1-p.res$par[1])*spmatrix_temp+(1-p.res$par[1])*(1-p.res$par[3])*phylomatrix+p.res$par[1]*diag(dim(phylomatrix)[1])))

# --- Maximum likelihood fits for all models ---
ml_fits <- mclapply(
  X = models,
  FUN = function(f) ml_fit(p=p.res$par,formula=f,data=A3b,spmatrix=spmatrix,phylomatrix=phylomatrix),
  mc.cores = 8
)
# Save models
saveRDS(ml_fits, "output/A3b/ml_fits.RDS") 

# --- Model Comparison ---
# Make empty NA list of coefficients and their respective p values in tuple form
n <- length(ml_fits)
ml_summaries <- list("bordering_language_richness" = as.list(rep(NA,n)),
                     "Range.Size..km2." = as.list(rep(NA,n)),
                     "Distance.to.Mainland" = as.list(rep(NA,n)), 
                     "Distance.to.Continent" = as.list(rep(NA,n)),
                     "L1_pop" = as.list(rep(NA,n)))
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
    ml_MSE <- c(ml_MSE, sum(ml_fits[[i]]$residuals**2)/dim(A3b)[1])
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
  ml_MSE <- c(ml_MSE, sum(ml_fits[[i]]$residuals**2)/dim(A3b)[1])
}

# Create new data frame from lists
ml_data <- data.frame(BIC = ml_BIC,
                      logLik = ml_logLik,
                      MSE = ml_MSE,
                      unlist(ml_summaries$bordering_language_richness),
                      unlist(ml_summaries$Range.Size..km2.),
                      unlist(ml_summaries$Distance.to.Mainland),
                      unlist(ml_summaries$Distance.to.Continent),
                      unlist(ml_summaries$L1_pop),
                      stringsAsFactors = FALSE)
colnames(ml_data) <- c("BIC", "logLik", "MSE", predictors)
# Sort based on increasing BIC
ml_data <- ml_data[order(ml_data$BIC), ]
# Add delta BIC column (min_BIC - current_BIC)
min_BIC <-  min(ml_data$BIC)
ml_data$deltaBIC <- ml_data$BIC - min_BIC
# Save output
write.csv(ml_data, file = "output/A3b/ml_data.csv", row.names=FALSE)

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

