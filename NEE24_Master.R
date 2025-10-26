library(nloptr)
library(ggplot2)
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
# Rename predictors
colnames(nee24)[colnames(nee24) == "L1.Population"] <- "L1_pop"

# --- Adjusting nee24 data with matrices ---
common_ids <- Reduce(intersect, list(
  nee24$ISO,
  rownames(phylomatrix),
  rownames(spmatrix)
))
nee24 <- nee24[nee24$ISO %in% common_ids, ]
phylomatrix <- phylomatrix[common_ids, common_ids]
spmatrix <- spmatrix[common_ids, common_ids]

# Preview and save data
paste("Size of dataset adjusted:", dim(nee24)[1])
paste("#Island_Endemic: ", sum(nee24$Island.Endemic))
paste("#Island_Endemic to Total Ratio: ", sum(nee24$Island.Endemic) / dim(nee24)[1])
# save data 
write.csv(nee24, file = "output/NEE24/nee24_adjusted.csv", row.names=FALSE)

# Preview and save matrices
print(dim(phylomatrix))
print(dim(spmatrix))
print(phylomatrix[1:10,1:10])
print(spmatrix[1:10,1:10])
save(phylomatrix, file = "output/NEE24/phylomatrix.RData") 
save(spmatrix, file = "output/NEE24/spmatrix.RData")

# --------
# ANALYSIS
# --------

# Get raw NEE24 adjusted data and matrices
raw_NEE24 <- read.csv("output/NEE24/NEE24_adjusted.csv")
load("output/NEE24/phylomatrix.RData")
load("output/NEE24/spmatrix.RData")

# --- Resizing matrices and dataset to remove too-similar languages ---
# phylo_dist <- dist(phylomatrix, method = "euclidean")
# phylo_dist <- as.matrix(phylo_dist)
# threshold <- 1.3
# omit_lang <- which(phylo_dist < threshold & phylo_dist > 0, arr.ind = T)[,1]
# raw_NEE24 <- raw_NEE24[-omit_lang,]
# spmatrix <- spmatrix[-omit_lang,
#                      -omit_lang]
# phylomatrix <- phylomatrix[-omit_lang,
#                            -omit_lang]

# --- Log transform variables ---
NEE24 <- raw_NEE24
NEE24$Phoneme.Inventory.Size <- log(NEE24$Phoneme.Inventory.Size) 
NEE24$Range.Size..km2. <- log(NEE24$Range.Size..km2.)
NEE24$L1_pop <- log(NEE24$L1_pop + 0.5)
# Do not transform zero entries for distance variables
NEE24$Distance.to.Mainland[NEE24$Distance.to.Mainland != 0] <- log(
  NEE24$Distance.to.Mainland[NEE24$Distance.to.Mainland != 0])
NEE24$Distance.to.Continent[NEE24$Distance.to.Continent != 0] <- log(
  NEE24$Distance.to.Continent[NEE24$Distance.to.Continent != 0])

# --- Preview transformed data ---
head(NEE24)
summary(NEE24)
# PIS ~ Range_Size, Grouping by Island Endemics
PIS_Range <- ggplot(data = NEE24,aes(x = Range.Size..km2., y = Phoneme.Inventory.Size, colour = as.logical(Island.Endemic))) +
  geom_point() + 
  guides(colour = guide_legend(title = "Island Endemic"))
print(PIS_Range)
ggsave(
  filename = "output/NEE24//PIS_Range.png",
  plot = PIS_Range,
  scale = 1,
  width=7,
  height=5
)

# # --- Testing GLS (2) ---
# p = c(0.5,0.5,0.5)
# spmatrix_test <- spmatrix/max(spmatrix)
# spmatrix_test <- exp(-(spmatrix_test/p[2])^2)
# mat <- as.matrix((p[3]*(1-p[1])*spmatrix_test+(1-p[1])*(1-p[3])*phylomatrix+p[1]*diag(dim(phylomatrix)[1])))
# res <- gls(model=Phoneme.Inventory.Size~1,data=NEE24,correlation=corSymm(mat[lower.tri(mat)],fixed=T),method="ML")
# print(-res$logLik)

# # --- Preview correlations ---
# pairs(NEE24[,
#             c("Phoneme.Inventory.Size", "Range.Size..km2.",
#               "L1.Population", "Distance.to.Mainland", "Distance.to.Continent")],
#       panel = function(x, y, ...) {
#         points(x, y, cex = 0.1, ...) # cex = 0.5 makes points half the default size
#       })
# cor(NEE24[,
#           c("Phoneme.Inventory.Size", "Range.Size..km2.",
#             "L1.Population", "Distance.to.Mainland", "Distance.to.Continent")])
# cor.test(NEE24$Phoneme.Inventory.Size,NEE24$Island.Endemic)

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
predictors <- c("Island.Endemic", "Range.Size..km2.",
                "Distance.to.Mainland", "Distance.to.Continent",
                "L1_pop")
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
unlist(models)

# --- Best p values for all models ---
p.res <- sbplx(c(0.5, 0.5, 0.5),
               best_p,
               formula=Phoneme.Inventory.Size~1,
               data=NEE24,
               spmatrix=spmatrix,phylomatrix=phylomatrix,
               lower=c(0,0,0),upper=c(1,1,1),
               nl.info = TRUE)
# Save p.res and covariance matrix
save(p.res, file = "output/NEE24/ml_p_res.RData")
spmatrix_temp <- spmatrix/max(spmatrix)
spmatrix_temp <- exp(-(spmatrix_temp/p.res$par[2])^2)
mat <- as.matrix((p.res$par[3]*(1-p.res$par[1])*spmatrix_temp+(1-p.res$par[1])*(1-p.res$par[3])*phylomatrix+p.res$par[1]*diag(dim(phylomatrix)[1])))

# --- Maximum likelihood fits for all models ---
ml_fits <- mclapply(
  X = models,
  FUN = function(f) ml_fit(p=p.res$par,formula=f,data=NEE24,spmatrix=spmatrix,phylomatrix=phylomatrix),
  mc.cores = 8
)
# Save models
saveRDS(ml_fits, "output/NEE24/ml_fits.RDS") 

# --- Model Comparison ---
# Make empty NA list of coefficients and their respective p values in tuple form
n <- length(ml_fits)
ml_summaries <- list("Island.Endemic" = as.list(rep(NA,n)),
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
    ml_MSE <- c(ml_MSE, sum(ml_fits[[i]]$residuals**2)/dim(NEE24)[1])
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
  ml_MSE <- c(ml_MSE, sum(ml_fits[[i]]$residuals**2)/dim(NEE24)[1])
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
                      stringsAsFactors = FALSE)
colnames(ml_data) <- c("BIC", "logLik", "MSE", predictors)
# Sort based on increasing BIC
ml_data <- ml_data[order(ml_data$BIC), ]
# Add delta BIC column (min_BIC - current_BIC)
min_BIC <-  min(ml_data$BIC)
ml_data$deltaBIC <- ml_data$BIC - min_BIC
# Save output
write.csv(ml_data, file = "output/NEE24/ml_data.csv", row.names=FALSE)

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

