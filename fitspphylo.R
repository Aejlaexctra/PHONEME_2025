#functions to correct for both phylo and spatial autocorrelation
	#formula: regression formula: response variable ~ predictor1 + predictor2 + ...
	#data: each row is the response variable and the independent variables of a lanugage
	#spmatrix: diagonal is 0 and the ij-th offdiagonal is geographic distance between language i and language j, assuming decrease in spatial autocorrelation using Gaussian:exp(-scale factor*distance))^2; 
	#phylomatrix: phylogenetic correlation matrix
	#starting value for maximization, p[1] is the total contribution of autocorrelation, p[3] is the relative contribution of spatial versus phylogenetic autocorrelation, p[2] is the scaling factor in the Gaussian function to model spatial autocorrelation, the smaller p[2] is , the faster spatial autocorrelation decreases with spatial distance, the default starting value is c(0.5, 0.5, 0.5). If the model has convergence problem, then you probably need to try different p iniital values. For example, for nested model, you can use the resulting p value of a simpler model as the initial p value for a more complex model.
#The outputs are:
	#p.res: the maximum likelihood estimates of p value
	#lm.res: the fitting result of the regression model under the maximum likelihood estimates of p value

# Including p values and log likelihood history
fitspphylo <- function (formula,data,spmatrix,phylomatrix,p=c(0.5, 0.5, 0.5)) {
  
  p1_history <<- c()
  p2_history <<- c()
  p3_history <<- c()
  logLik_history <<- c()
  
  # Store optimisation history
  opt_data <- data.frame(p1 = numeric(),
                         p2 = numeric(),
                         p3 = numeric(),
                         logLik = numeric(),
                         stringsAsFactors = FALSE)
  
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
    p1_history <<- c(p1_history,p[1])
    p2_history <<- c(p2_history,p[2])
    p3_history <<- c(p3_history,p[3])
    logLik_history <<- c(logLik_history, out)
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
  p.res <- sbplx(p,cal,formula=formula,data=data,spmatrix=spmatrix,phylomatrix=phylomatrix,lower=c(0,0,0),upper=c(1,1,1),nl.info = TRUE)
  lm.res <- cal2(p=p.res$par,formula=formula,data=data,spmatrix=spmatrix,phylomatrix=phylomatrix)
  
  # Create new data frame from lists
  new_opt <- data.frame(p1 = p1_history,
                        p2 = p2_history,
                        p3 = p3_history,
                        logLik = logLik_history,
                        stringsAsFactors = FALSE)
  # Combine with existing data
  opt_data <- rbind(opt_data, new_opt)
  
  list(p.res=p.res,lm.res=lm.res,opt_data=opt_data)
}