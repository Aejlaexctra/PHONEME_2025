library(nloptr)
library(nlme)
library(parallel)

# Set random seed
set.seed(2077)

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
# Save output
write.csv(ml_data, file = "output/anderson/ml_data.csv", row.names=FALSE)
