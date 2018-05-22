# Helper functions for nardl

#-------------------------------------------------------------------------------
# Function trimr
#
# Became redundant after rewriting lagm
# trimr <- function(x,rb,re) {
#   # Performance improvement over cbind
#   # x <- cbind(x)
#   if(!is.matrix(x)) x <- matrix(x, ncol=1)
#
#   n <- nrow(x)
#   if ((rb+re) >= n) {
#     stop('Attempting to trim too much')
#   }
#   z <- x[(rb+1):(n-re),]
#   return(z)
# }

#-------------------------------------------------------------------------------
# Function lagm
# INPUT: m is a matrix, nLags the number of lags
lagm <- function(m, nLags) {
  # JM: This code is redundant. More than 2 arguments will cause an error
  # anyway, less than 2 arguments will cause an error about missing
  # arguments that's more informative than yours.

  # nargin <- length(as.list(match.call())) - 1
  # if (nargin != 2) {
  #   stop('Check function inputs')
  # }

  if(!is.matrix(m))
    stop("Trying to lag something that's not a matrix")

  d <- dim(m)

  #Add column names if they don't exist yet
  if(is.null(colnames(m)))
    colnames(m) <- as.character(seq_len(d[2]))

  #Check
  if(nLags > d[1])
    stop(sprintf("You try to create %d lags but there's only %d rows in m.",
                 nLags, d[1]))

  lagM <- matrix(NA,nrow=d[1], ncol = d[2]*nLags)

  for(i in seq_len(nLags)){
    #Make ids for the columns in result
    cid <- seq(1:d[2]) + d[2]*(i-1)

    lagM[(i+1):d[1],cid] <- m[1:(d[1]-i),]
  }

  cnames <- outer(colnames(m),seq_len(nLags), FUN = paste, sep = "_")

  colnames(lagM) <- c(cnames)

  return(lagM)

}

#-------------------------------------------------------------------------------
#' ARCH test
#'
#'Computes the Lagrange multiplier test for conditional heteroscedasticity of Engle (1982), as described by Tsay (2005, pp. 101-102).
#'
#'@param x numeric vector
#'@param lags positive integer number of lags
#'@param demean logical: If TRUE, remove the mean before computing the test statistic.
#'@importFrom stats pchisq embed lm
#'@examples
#'
#'reg<-nardl(food~inf,fod,ic="aic",maxlags = TRUE,graph = TRUE,case=3)
#'x<-reg$selresidu
#'nlag<-reg$np
#'ArchTest(x,lags=nlag)
#'
#'@export

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
#' LM test for serial correlation
#'
#'
#'@param object fitted lm model
#'@param nlags positive integer number of lags
#'@param fill starting values for the lagged residuals in the auxiliary regression. By default 0.
#'@param type Fisher or Chisquare statistics
#'@importFrom stats lm pchisq
#'@importFrom gtools na.replace
#'@examples
#'
#'reg<-nardl(food~inf,fod,ic="aic",maxlags = TRUE,graph = TRUE,case=3)
#'lm2<-bp2(reg$fit,reg$np,fill=0,type="F")
#'
#'@export

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



#-------------------------------------------------------------------------------
# Function seqa
seqa<-function(a,b,c){
  #seq=(a:b:(a+b*(c-1)))';
  se<-seq(a,(a+b*(c-1)),by=b)
  return(t(se))
}

