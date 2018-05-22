#' Function cusum
#'
#' This function creates two plots for the cusum test.
#'
#' @param x is a object of the class \code{\link{nardl}} or a vector
#' with the recursive errors
#' @param k is the estimated coefficients length
#' @param n is the recursive errors length
#' @param ... arguments passed to other methods. Currently ignored.
#' @importFrom graphics abline legend lines par plot
#' @examples
#'
#' reg<-nardl(food~inf,fod,ic="aic",maxlags = TRUE,graph = TRUE,case=3)
#' cusum(reg)
#'
#' # Or by hand
#' e<-reg$rece
#' k<-reg$k
#' n<-reg$n
#' cusum(x=e,k=k,n=n)
#'
#'@export
cusum <- function(x, ...){
  UseMethod("cusum")
}

#' @rdname cusum
#' @export
cusum.nardl <- function(x, ...){
  e <- x$rece
  k <- x$k
  n <- x$n
  cusum.default(e,k,n)
}

#' @rdname cusum
#' @export
cusum.default <- function(x,k,n, ...){
  w<-as.matrix(na.omit(x))
  #n<-length(x)
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
