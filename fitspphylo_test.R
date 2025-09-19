set.seed(42)

# Prepare scripts, models and data
load("Data/models/phoneme_glm_models.Rdata")
source("fitspphylo.R") # (1) The ecological drivers of variation in global language diversity
data <- read.csv("Data/phoneme_count/phoneme_data.csv")
load("Data/phoneme_count/phoneme_spmatrix.Rdata")
load("Data/phoneme_count/phoneme_phylomatrix.Rdata")

fitspphylo <- function (formula,data,spmatrix,phylomatrix,p=c(0.5, 0.5, 0.5)) {
  cal <- function (p,formula,data,spmatrix,phylomatrix) {
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
    # print(out)
    -out
  }
  cal2 <- function (p,formula,data,spmatrix,phylomatrix) {
    spmatrix <- spmatrix/max(spmatrix)
    spmatrix <- exp(-(spmatrix/p[2])^2)
    mat <- as.matrix((p[3]*(1-p[1])*spmatrix+(1-p[1])*(1-p[3])*phylomatrix+p[1]*diag(dim(phylomatrix)[1])))
    res <- try(gls(model=formula,data=data,correlation=corSymm(mat[lower.tri(mat)],fixed=T),method="ML"),silent=T)
    res
  }
  library(nloptr)
  library(nlme)
  p.res <- sbplx(p,cal,formula=formula,data=data,spmatrix=spmatrix,phylomatrix=phylomatrix,lower=c(0,0,0),upper=c(1,1,1),nl.info = TRUE, 
                 nl.opts(list(xtol_rel = 1e-6, maxeval = 2000)))
  # p.res <- optim(par=p,fn=cal, method = "L-BFGS-B",formula=formula,data=data,spmatrix=spmatrix,phylomatrix=phylomatrix,lower=c(0,0,0),upper=c(1,1,1))
  lm.res <- cal2(p=p.res$par,formula=formula,data=data,spmatrix=spmatrix,phylomatrix=phylomatrix)
  list(p.res=p.res,lm.res=lm.res)
}

# Previous function
autoglm <- function (a, y, X, Wsp, Wphy, fml, step) {
  # W <- a[2]*(a[1]*Wsp+(1-a[1])*Wnb)+(1-a[2])*Wphy
  W <- a*Wsp+(1-a)*Wphy
  Wy <- W %*% y  # calculates a weighted sum of y using W
  X2 <- W %*% X  # calculates a weighted sum of X using W
  X2 <- cbind(X, X2)  # concatenates X and X2 as columns in a new matrix
  res <- lm(Wy~X2)$residuals  # fits a linear model of Wy on X2 and extracts the residuals
  X1 <- cbind(y,Wy, X, res)  # concatenates Wy, X, and the residuals as columns in a new matrix
  colnames(X1)[1:2] <- c('LOG_Sounds', 'weighted_sums')
  # print(head(X1))
  
  if(step == 1){
    out2 <- glm(formula = fml, data = as.data.frame(X1), family = gaussian())
    # summary(out2)
    #-out$loglik  # returns the negative log-likelihood of the model
    out2 <- -logLik(out2)[[1]]
  } else if(step == 2){
    out2 <- glm(formula = fml, data = as.data.frame(X1), family=gaussian())
  } else {
    out2 <- glm(formula = fml, data = as.data.frame(X1), family=gaussian())
    out2 <- BIC(out2)
  }
  return(out2)
}

# Predictors matrix
X <- data.frame(
  LOG_L1_pop = data$LOG_L1_pop,
  LOG_Area =data$LOG_Area,
  LOG_Altitude_Range = data$LOG_Altitude_Range,
  LOG_Distance.to.Continent = data$LOG_Distance.to.Continent,
  LOG_Distance.to.Mainland = data$LOG_Distance.to.Mainland,
  LOG_Roughness = data$LOG_Roughness,
  LOG_Bordering_Language_Richness = data$LOG_Bordering_Language_Richness,
  documentation = data$documentation,
  Island.Endemic = data$Island.Endemic
)

X <- as.matrix(X)
# Response variable
y <- data$LOG_Sounds

# Define GLM dataframe
glm_data <- data.frame(Model_ID = integer(),
                       Coefficients = list(),
                       BIC = numeric(),
                       AIC = numeric(),
                       logLik = numeric(),
                       Fitted = list(),
                       stringsAsFactors = FALSE)

# Define columns as lists
Model_IDs <- vector("list", length(phoneme_glm_models))
Coefficients <- vector("list", length(phoneme_glm_models))
BICs <- vector("list", length(phoneme_glm_models))
AICs <- vector("list", length(phoneme_glm_models))
Fitted <- vector("list", length(phoneme_glm_models))

# Load and add data from each model
for (i in 1:length(phoneme_glm_models))
{
  # Get current model
  formula = phoneme_glm_models[[i]]
  
  # Get GLM fit
  best_a <- try(optim((1), autoglm, method="Brent", lower=c(0), upper=c(1), y=y, X=X, Wsp=phoneme_spmatrix, Wphy=phoneme_phylomatrix, fml = formula, step = 1))
  
  print(best_a$par)
  
  glm_model <- autoglm(a = best_a$par, y=y, X=X, Wsp=phoneme_spmatrix, Wphy=phoneme_phylomatrix, fml = formula, step = 2)
  Model_IDs[[i]] <- i
  BICs[[i]] <- BIC(glm_model)
  AICs[[i]] <- AIC(glm_model)
  # Get summary
  # Extract coefficients with p-values
  Coefficients[[i]] <- summary(glm_model)$coefficients
  Fitted[[i]] <- glm_model$fitted.values
}

# Create new data frame from lists
new_glm <- data.frame(Model_ID = unlist(Model_IDs),
                      BIC = unlist(BICs),
                      AIC = unlist(AICs),
                      stringsAsFactors = FALSE)

# Add coefficients, model fits and true values as list-columns
new_glm$Coefficients <- Coefficients
new_glm$Fitted <- Fitted

# Combine with existing data
glm_data <- rbind(glm_data, new_glm)

# Sort entries based on ascending BICs
glm_data <- glm_data[order(glm_data$BIC), ] 

# Show first 5 entries
head(glm_data[,c("Model_ID","BIC")])

# Show top BIC ranked coefficients
print("BIC Ranked")
glm_data$Coefficients[1]

