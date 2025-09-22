# Set random seed
set.seed(2077)

# Get raw nee24 adjusted data and matrices
raw_nee24 <- read.csv("data/nee24_adjusted.csv")
load("data/phylomatrix.RData")
load("data/spmatrix.RData")

# --- Log transform variables ---
nee24 <- raw_nee24
nee24$Phoneme.Inventory.Size <- log(nee24$Phoneme.Inventory.Size) 
nee24$Range.Size..km2. <- log(nee24$Range.Size..km2.)
nee24$L1.Population <- log(nee24$L1.Population + 0.5)
nee24$Distance.to.Mainland <- log(nee24$Distance.to.Mainland + 0.5)
nee24$Distance.to.Continent <- log(nee24$Distance.to.Continent + 0.5)

# Preview transformed data
head(nee24)
summary(nee24)

# --- Preview correlations ---
pairs(nee24[,
            c("Phoneme.Inventory.Size", "Range.Size..km2.",
              "L1.Population", "Distance.to.Mainland", "Distance.to.Continent")],
      panel = function(x, y, ...) {
        points(x, y, cex = 0.1, ...) # cex = 0.5 makes points half the default size
      })
cor(nee24[,
          c("Phoneme.Inventory.Size", "Range.Size..km2.",
            "L1.Population", "Distance.to.Mainland", "Distance.to.Continent")])
cor.test(nee24$Phoneme.Inventory.Size,nee24$Island.Endemic)

# --- Testing GLS (2) ---
p = c(0.5,0.5,0.5)
spmatrix_test <- spmatrix/max(spmatrix)
spmatrix_test <- exp(-(spmatrix_test/p[2])^2)
mat <- as.matrix((p[3]*(1-p[1])*spmatrix_test+(1-p[1])*(1-p[3])*phylomatrix+p[1]*diag(dim(phylomatrix)[1])))
res <- gls(model=Phoneme.Inventory.Size~1,data=nee24,correlation=corSymm(mat[lower.tri(mat)],fixed=T),method="ML")
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
