# Set seed
set.seed(9291)

## read in data
rawdata <- read.csv('data/nee24_adjusted.csv',header=T)

#run eigenvector spatial filtering
library(spfilteR)
library(parallel)

#define matrices for pairwise spatial similarity and phylogenetic similarity
#by definition, the diagonal elements are zero and each row is normalized to 1
#use the matrices derived in Hua et al. 2022
# load("data/Wnb.Rdata") #if two langauges overlap
# load("Wphy.Rdata") #similarity based on phylogenetic correlation
# load("Wsp.Rdata") #similarity based on spatial distance-decay in correlation

# Load matrices
load("data/phylomatrix.RData")
load("data/spmatrix.RData")
load("data/Wnb.Rdata")
Wphy <- phylomatrix 
Wsp <- spmatrix

# Set phylomatrix diagonal to 0
diag(Wphy) <- 0

# Find languages which have zero rows in matrices (Wphy and Wnb)
length(which(rowSums(Wnb) == 0))
length(which(rowSums(Wphy) == 0))
length(which(rowSums(Wsp) == 0))
zero_row_languages <- union(which(rowSums(Wnb) == 0),which(rowSums(Wphy) == 0))
paste("Number of languages to be removed:", length(zero_row_languages))

# Remove zero row languages from data and matrices
rawdata <- rawdata[-zero_row_languages,]
Wphy <- Wphy[-zero_row_languages,-zero_row_languages]
Wsp <- Wsp[-zero_row_languages,-zero_row_languages]
Wnb <- Wnb[-zero_row_languages,-zero_row_languages]
paste("Size of dataset after zero row languages removed:", dim(rawdata)[1])

# Row standardise matrices used 
Wphy <- Wphy / rowSums(Wphy)
Wsp <- Wsp / rowSums(Wsp)
Wnb <- Wnb / rowSums(Wnb)

# # Ensure that NA valeus in matrices are replaced by zero (due to zero rows during row-standardisation)
# Wphy[is.na(Wphy)] <- 0
# Wsp[is.na(Wsp)] <- 0
# Wnb[is.na(Wnb)] <- 0

#extract eigenvectors
#using SAR model
EVs.nb <- getEVs(W=Wnb,covars=NULL)
EVs.phy <- getEVs(W=Wphy,covars=NULL)
EVs.sp <- getEVs(W=Wsp,covars=NULL)
Ec.nb <- EVs.nb$moran/max(EVs.nb$moran)>=.25 #the threshold 0.25 is suggested by Griffith (2003)
Ec.phy <- EVs.phy$moran/max(EVs.phy$moran)>=.25
Ec.sp <- EVs.sp$moran/max(EVs.sp$moran)>=.25
evecs <- cbind(EVs.nb$vectors[,Ec.nb],EVs.phy$vectors[,Ec.phy],EVs.sp$vectors[,Ec.sp])

#following Koplenig 2024, using intercept-only model
#eigenvectors are selected by the objective: eigenvectors are added to the model until no more 
#	eigenvectors can further reduce the residual autocorrelation. 

########################
#polysynthesis analysis
########################
#defind variables
variables <- c("Phoneme.Inventory.Size", "L1.Population", "Range.Size..km2.", 
               "Island.Endemic", "Distance.to.Mainland", "Distance.to.Continent"
)

data <- rawdata[,variables]

#log transform highly right skewed predictors
data$Phoneme.Inventory.Size <- log(data$Phoneme.Inventory.Size) 
data$Range.Size..km2. <- log(data$Range.Size..km2.)
data$L1.Population <- log(data$L1.Population + 0.5)
data$Distance.to.Mainland <- log(data$Distance.to.Mainland + 0.5)
data$Distance.to.Continent <- log(data$Distance.to.Continent + 0.5)

#scaling continuous data to help interpreting results
data$L1.Population<- scale(data$L1.Population)
data$Range.Size..km2. <- scale(data$Range.Size..km2.)
data$Distance.to.Mainland <- scale(data$Distance.to.Mainland)
data$Distance.to.Continent <- scale(data$Distance.to.Continent)
data$Phoneme.Inventory.Size <- scale(data$Phoneme.Inventory.Size)

#define response variable
y <- data$Phoneme.Inventory.Size

#define similarity matrix to calculate Moran's I
W <- list(Wnb,Wphy,Wsp)

#define intercept
n <- length(y)
x <- as.matrix(rep(1, n))

#define Moran's I^2 test, using the general formula of MC with extension for multiple matrices
MI.resid <- function(resid, W, alternative = "greater") {
  
  n <- length(resid)
  sig2 <- t(resid) %*% resid / n
  
  if (length(W)==1) {
    W <- W[[1]]
    I1 <- t(resid) %*% W %*% resid / sig2
    I2 <- sqrt(sum(diag((t(W)+W)%*%W)))
    Isq <- (I1/I2)^2
    p <- 1-pchisq(Isq,df=1)
  } else if (length(W)>1) {
    q <- length(W)
    I1 <- numeric(q)
    I2 <- matrix(,q,q)
    for (i in 1:q) {
      I1[i] <- t(resid) %*% W[[i]] %*% resid / sig2
      for (j in 1:q) {
        I2[i,j] <- sum(diag((W[[i]]+t(W[[i]])) %*% (W[[j]]+t(W[[j]])))) / 2
      }
    }
    I2 <- solve(I2)
    Isq <- t(I1) %*% I2 %*% I1
    p <- 1-pchisq(Isq,df=q)
  }
  list(Isq=Isq,p=p)
}

#update objective function with general formula of MC with extension for multiple matrices
objfunc <-function(y, xe, W) {
  MI.resid <- function(resid, W, alternative = "greater") {
    
    n <- length(resid)
    sig2 <- t(resid) %*% resid / n
    
    if (length(W)==1) {
      W <- W[[1]]
      I1 <- t(resid) %*% W %*% resid / sig2
      I2 <- sqrt(sum(diag((t(W)+W)%*%W)))
      Isq <- (I1/I2)^2
      p <- 1-pchisq(Isq,df=1)
    } else if (length(W)>1) {
      q <- length(W)
      I1 <- numeric(q)
      I2 <- matrix(,q,q)
      for (i in 1:q) {
        I1[i] <- t(resid) %*% W[[i]] %*% resid / sig2
        for (j in 1:q) {
          I2[i,j] <- sum(diag((W[[i]]+t(W[[i]])) %*% (W[[j]]+t(W[[j]])))) / 2
        }
      }
      I2 <- solve(I2)
      Isq <- t(I1) %*% I2 %*% I1
    }
    return(Isq)
  }
  o <- glm(y~xe-1, family=binomial(link=probit), maxit = 100)
  resid <- (y - o$fitted.values) / sqrt(o$fitted.values * (1 - o$fitted.values)) 
  test <- MI.resid(resid=resid, W=W)
  return(test)
}

#calculate initial AIC and Moran's I^2=
init <- glm(y~x-1, family=gaussian, maxit = 100)
AIC_init <- init$aic
Rsquared_init <- 1-(init$deviance/init$null.deviance)
resid_init <- (y - init$fitted.values) / sqrt(init$fitted.values * (1 - init$fitted.values)) #pearson residual
MI_init <- MI.resid(resid=resid_init,W=W)

#search algorithm based on model fit, using AIC with minimum reduction =0.01 as criteria
min.reduction <- 0.01
IC <- AIC_init
mindiff <- abs(IC * min.reduction)
sel <- as.logical(rep(1,dim(evecs)[2]))
sel_id <- NULL
selset <- which(sel)

for (i in which(sel)) {
  
  ref <- Inf
  sid <- NULL
  
  for (j in selset) {
    xe <- cbind(x, evecs[, sel_id], evecs[, j])
    test <- glm(y~xe-1, family=gaussian, maxit = 1000)
    test <- test$aic
    if (test < ref) {
      sid <- j
      ref <- test
    }
  }
  
  if (ref < IC & abs(IC - ref) >= mindiff) {
    IC <- ref
    sel_id <- c(sel_id, sid)
    mindiff <- abs(IC * min.reduction)
  } else {
    break
  }
  
  # remove selected eigenvectors from candidate set
  selset <- selset[!(selset %in% sel_id)]
}
sel_id
# end selection

#search algorithm based on residual autocorrelation.
oldpMI <- MI_init$p
oldzMI <- MI_init$Isq
sel <- as.logical(rep(1,dim(evecs)[2]))
sel_id <- NULL
selset <- which(sel)

for (i in which(sel)) {
  
  ref <- Inf
  sid <- NULL
  
  f <- function (j,y,x,W,evecs,sel_id) {
    o <- try(objfunc(y = y, xe = cbind(x, evecs[, sel_id], evecs[, j]), W = W),silent=T)
    if (class(o)=="try-error" || is.na(o)) {
      o <- Inf
    }
    o
  }
  
  test <- unlist(mclapply(selset, f, y, x, W, evecs, sel_id, mc.cores=10))
  
  sid <- which(test==min(test))
  ref <- test[sid]
  
  # stopping rules
  pMI <- 1-pchisq(ref,df=3)
  if (pMI > oldpMI || ref < oldzMI) {
    oldzMI <- ref
    oldpMI <- pMI
    sel_id <- c(sel_id, selset[sid])
  } else {
    break #adding more eigenvectors does not further reduce Moran's I
  }
  
  # remove selected eigenvectors from candidate set
  selset <- selset[!(selset %in% sel_id)]
}

# filtered regression
xev <- cbind(x,evecs[, sel_id])
out <- glm(y~xev-1, family=binomial(link=probit), maxit = 100)
AIC_out <- out$aic
Rsquared_out <- 1-(out$deviance/out$null.deviance)
resid_out <- (y - out$fitted.values) / sqrt(out$fitted.values * (1 - out$fitted.values)) #pearson residual
MI_out <- MI.resid(resid=resid_out,W=W)

#using rjMCMC to estimate PIP
#generate datasets for gregl
data <- cbind(data,ev=evecs[, sel_id])
write.csv(data,file="data/NEE24_polydata_parma.csv")

#install gretl
#open the csv file in gretl
