library(nloptr)
library(nlme)
library(parallel)

# Set random seed
set.seed(2077)

# Get raw A2 adjusted data and matrices
raw_A2 <- read.csv("output/A2/A2_adjusted.csv")
load("output/A2/phylomatrix.RData")
load("output/A2/spmatrix.RData")

# --- Log transform variables ---
A2 <- raw_A2
A2$Phoneme.Inventory.Size <- log(A2$Phoneme.Inventory.Size) 
A2$area <- log(A2$area)
A2$altitude_range <- log(A2$altitude_range + 0.5)
A2$L1_pop <- log(A2$L1_pop + 0.5)
A2$bordering_language_richness <- log(A2$bordering_language_richness + 0.5)

# Preview transformed data
head(A2)
summary(A2)

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

# # --- Testing GLS (2) ---
# p = c(0.5,0.5,0.5)
# spmatrix_test <- spmatrix/max(spmatrix)
# spmatrix_test <- exp(-(spmatrix_test/p[2])^2)
# mat <- as.matrix((p[3]*(1-p[1])*spmatrix_test+(1-p[1])*(1-p[3])*phylomatrix+p[1]*diag(dim(phylomatrix)[1])))
# res <- gls(model=Phoneme.Inventory.Size~1,data=A2,correlation=corSymm(mat[lower.tri(mat)],fixed=T),method="ML")
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
predictors <- c("Island.Endemic", "area",
                "altitude_range", "bordering_language_richness",
                "L1_pop")
response <- "Phoneme.Inventory.Size"
n_pred <- length(predictors)
pred_combs <- sapply(1:n_pred, function(x) combn(predictors, x))
models <- c()
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
    models <- c(models,model)
  }
}
unlist(models)

# --- Best p values for all models ---
p.res <- sbplx(c(0.5, 0.5, 0.5),
               best_p,
               formula=Phoneme.Inventory.Size~1,
               data=A2,
               spmatrix=spmatrix,phylomatrix=phylomatrix,
               lower=c(0,0,0),upper=c(1,1,1),
               nl.info = TRUE)
# Save p.res and covariance matrix
save(p.res, file = "output/A2/ml_p_res.RData")
spmatrix_temp <- spmatrix/max(spmatrix)
spmatrix_temp <- exp(-(spmatrix_temp/p.res$par[2])^2)
mat <- as.matrix((p.res$par[3]*(1-p.res$par[1])*spmatrix_temp+(1-p.res$par[1])*(1-p.res$par[3])*phylomatrix+p.res$par[1]*diag(dim(phylomatrix)[1])))
save(mat, file = "output/A2/mat.RData")

# --- Maximum likelihood fits for all models ---
ml_fits <- mclapply(
  X = models,
  FUN = function(f) ml_fit(p=p.res$par,formula=f,data=A2,spmatrix=spmatrix,phylomatrix=phylomatrix),
  mc.cores = 8
)
# Save models
saveRDS(ml_fits, "output/A2/ml_fits.RDS") 

# --- Model Comparison ---
# Make empty NA list of coefficients and their respective p values in tuple form
n <- length(ml_fits)
ml_summaries <- list("Island.Endemic" = as.list(rep(NA,n)),
                     "area" = as.list(rep(NA,n)),
                     "altitude_range" = as.list(rep(NA,n)), 
                     "bordering_language_richness" = as.list(rep(NA,n)),
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
    ml_MSE <- c(ml_MSE, sum(ml_fits[[i]]$residuals**2)/dim(A2)[1])
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
  ml_MSE <- c(ml_MSE, sum(ml_fits[[i]]$residuals**2)/dim(A2)[1])
}

# Create new data frame from lists
ml_data <- data.frame(BIC = ml_BIC,
                      logLik = ml_logLik,
                      MSE = ml_MSE,
                      unlist(ml_summaries$Island.Endemic),
                      unlist(ml_summaries$area),
                      unlist(ml_summaries$altitude_range),
                      unlist(ml_summaries$bordering_language_richness),
                      unlist(ml_summaries$L1_pop),
                      stringsAsFactors = FALSE)
colnames(ml_data) <- c("BIC", "logLik", "MSE", predictors)
# Sort based on increasing BIC
ml_data <- ml_data[order(ml_data$BIC), ]
# Add delta BIC column (min_BIC - current_BIC)
min_BIC <-  min(ml_data$BIC)
ml_data$deltaBIC <- ml_data$BIC - min_BIC
# Save output
write.csv(ml_data, file = "output/A2/ml_data.csv", row.names=FALSE)


