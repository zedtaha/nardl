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

  cnames <- paste(colnames(m),seq_len(nLags), sep = "_")
  if(d[2] > 1) cnames <- rep(cnames, each = d[2])
  colnames(lagM) <- cnames

  return(lagM)

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

#-------------------------------------------------------------------------------
# Function cumsq
cumsq<-function(e,k,n){
  w<-as.matrix(na.omit(e))
  w=cumsum(w^2)/sum(w^2)
  m=abs(0.5*(n-k)-1) #abs to avoid negative log
  c=0.74191-0.17459*log(m)-0.26526*(1/m)+0.0029985*m-0.000010943*m^2
  w2=c+(seqa(k,1,length(w))-k)/(n-k)
  w1=-c+(seqa(k,1,length(w))-k)/(n-k)
  x<-seqa(k,1,length(w))
  w1<-matrix(w1,length(w1),1)
  w<-matrix(w,length(w),1)
  w2<-matrix(w2,length(w2),1)
  grange<-range(w1,w2)
  par(mar = c(5,4,4,8))
  plot(x,w,main="CUSUM of Squares Test",type = "l",ylim = grange,xlab ="",ylab = "Emperical fluctuation process",col="blue")
  lines(x,w1,col="red")
  lines(x,w2,col="red")
  abline(h=0,lty=2)
  legend(par("usr")[2],par("usr")[4],
         xpd = TRUE ,
         bty = "n",
         c("CUSUM of squares","5% significance"),
         lty = c(1, 1),
         cex=0.6,
         col=c("blue","red") )

}

#-------------------------------------------------------------------------------
# Function cusum
cusum<-function(e,k,n){
  w<-as.matrix(na.omit(e))
  #n<-length(e)
  w=cumsum(w/apply(w, 2, sd))
  c=0.984
  w2=seqa(c*sqrt(n-k),(2*c*sqrt(n-k))/length(w),length(w))
  w1=seqa(-c*sqrt(n-k),(-2*c*sqrt(n-k))/length(w),length(w))
  x<-seqa(k,1,length(w))
  w1<-matrix(w1,length(w1),1)
  w<-matrix(w,length(w),1)
  w2<-matrix(w2,length(w2),1)
  grange<-range(w1,w2)
  par(mar = c(5,4,4,8))
  plot(x,w,main="CUSUM Test",type = "l",ylim = grange,xlab = "",ylab = "Emperical fluctuation process",col="blue")
  lines(x,w1,col="red")
  lines(x,w2,col="red")
  abline(h=0,lty=2)
  legend(par("usr")[2],par("usr")[4],
         xpd = TRUE ,
         bty = "n",
         c("CUSUM ","5% significance"),
         lty = c(1, 1),
         cex=0.6,
         col=c("blue","red") )

}

#-------------------------------------------------------------------------------
# Function seqa
seqa<-function(a,b,c){
  #seq=(a:b:(a+b*(c-1)))';
  se<-seq(a,(a+b*(c-1)),by=b)
  return(t(se))
}

