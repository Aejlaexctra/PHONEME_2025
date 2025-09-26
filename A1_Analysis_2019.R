library(nloptr)
library(nlme)
library(parallel)

# Set random seed
set.seed(2077)

# Get raw A1 adjusted data and matrices
raw_A1 <- read.csv("output/A1/A1_adjusted.csv")
load("output/A1/phylomatrix.RData")
load("output/A1/spmatrix.RData")

# --- Log transform variables ---
A1 <- raw_A1
A1$Phoneme.Inventory.Size <- log(A1$Phoneme.Inventory.Size) 
A1$L1_pop <- log(A1$L1_pop + 0.5)

# Preview transformed data
head(A1)
summary(A1)

# --- Preview correlations ---
plot(Phoneme.Inventory.Size ~ L1_pop,data=A1)

# --- Testing GLS (2) ---
p = c(0.5,0.5,0.5)
spmatrix_test <- spmatrix/max(spmatrix)
spmatrix_test <- exp(-(spmatrix_test/p[2])^2)
mat <- as.matrix((p[3]*(1-p[1])*spmatrix_test+(1-p[1])*(1-p[3])*phylomatrix+p[1]*diag(dim(phylomatrix)[1])))
res <- gls(model=Phoneme.Inventory.Size~1,data=A1,correlation=corSymm(mat[lower.tri(mat)],fixed=T),method="ML")
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

## Model predictor combinations
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
