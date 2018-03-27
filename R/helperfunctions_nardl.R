# Helper functions for nardl

#-------------------------------------------------------------------------------
# Function trimr
trimr <- function(x,rb,re) {
  x <- cbind(x)
  n <- nrow(x)
  if ((rb+re) >= n) {
    stop('Attempting to trim too much')
  }
  z <- x[(rb+1):(n-re),]
  return(z)
}

#-------------------------------------------------------------------------------
# Function lagm
lagm <- function(m, nLags) {
  # JM: This code is redundant. More than 2 arguments will cause an error
  # anyway, less than 2 arguments will cause an error about missing
  # arguments that's more informative than yours.

  # nargin <- length(as.list(match.call())) - 1
  # if (nargin != 2) {
  #   stop('Check function inputs')
  # }

  lagM <- c()
  for(i in seq(nLags)) {
    for(j in seq(ncol(m))) {
      tmp <- c(rep(NA, i), trimr(m[,j], 0, i))
      #colnames(tmp)<-colnames(m)

      lagM <- cbind(lagM, tmp)

    }
  }
  if(ncol(m)>=2){
    colnames(lagM)<-paste(colnames(m)[1:ncol(m)],rep(seq(nLags), each=ncol(m)),sep = "_")
  }
  else{
    colnames(lagM)<-paste(colnames(m),c(1:ncol(lagM)),sep = "_")
  }
  #seq(nLags*ncol(m))
  return(as.matrix(lagM))
}

#-------------------------------------------------------------------------------
# Function ArchTest

ArchTest <- function (x, lags=12, demean = FALSE)
{
  # Capture name of x for documentation in the output
  xName <- deparse(substitute(x))
  #
  x <- as.vector(x)
  if(demean) x <- scale(x, center = TRUE, scale = FALSE)
  #
  lags <- lags + 1
  mat <- embed(x^2, lags)
  arch.lm <- summary(lm(mat[, 1] ~ mat[, -1]))
  STATISTIC <- arch.lm$r.squared * length(resid(arch.lm))
  names(STATISTIC) <- "Chi-squared"
  PARAMETER <- lags - 1
  names(PARAMETER) <- "df"
  PVAL <- 1 - pchisq(STATISTIC, df = PARAMETER)
  METHOD <- "ARCH LM-test;  Null hypothesis:  no ARCH effects"
  result <- list(statistic = STATISTIC, parameter = PARAMETER,
                 p.value = PVAL, method = METHOD, data.name =
                   xName)
  class(result) <- "htest"
  return(result)
}

#-------------------------------------------------------------------------------
# Function ginv

# This function is not used anywhere
ginv <- function(x, tol = sqrt(.Machine$double.eps))
{
  ## Generalized Inverse of a Matrix
  dnx <- dimnames(x)
  if(is.null(dnx)) dnx <- vector("list", 2)
  s <- svd(x)
  nz <- s$d > tol * s$d[1]
  structure(
    if(any(nz)) s$v[, nz] %*% (t(s$u[, nz])/s$d[nz]) else x,
    dimnames = dnx[2:1])
}

#-------------------------------------------------------------------------------
# Function bp2
bp2<-function(object,nlags,fill=NULL,type=c("F","Chi2")){
  e<-as.matrix(object$residuals)
  n<-nrow(e)
  x<-as.matrix(object$model[,-1])
  xx<-cbind(x,lagm(e,nlags))
  if(fill==0){
    xx<-na.replace(xx,0)
    u<-lm(e~xx)
    r2<-summary(u)$r.squared
    LM<-n*r2
    pv<-pchisq(LM,nlags,lower.tail = FALSE)

  }
  else{
    u<-lm(e~xx)
    r2<-summary(u)$r.squared
    LM<-(n-nlags)*r2
  }
  if(type=="Chi2"){pv<-pchisq(LM,nlags,lower.tail = FALSE)}
  if(type=="F"){pv<-pf(LM,nlags,1,lower.tail = FALSE)}

  res<-cbind(LM,pv,nlags)
  colnames(res)<-c("LM","P-value","lags")
  return(res)
}


