#-------------------------------------------------------------------------------
#' Function cumsq
#'
#' @inheritParams cusum
#' @importFrom graphics abline legend lines par plot
#' @examples
#'
#' reg<-nardl(food~inf,fod,ic="aic",maxlags = TRUE,graph = TRUE,case=3)
#' cumsq(reg)
#'
#' # Or manually
#' e<-reg$rece
#' k<-reg$k
#' n<-reg$n
#' cumsq(x=e,k=k,n=n)
#'
#' @export
cumsq <- function(x, ...){
  UseMethod("cumsq")
}

#' @rdname cumsq
#' @export
cumsq.nardl <- function(x, ...){
  e <- x$rece
  k <- x$k
  n <- x$n
  cumsq(e, k, n)
}

#' @rdname cumsq
#' @export
cumsq.default <- function(x,k,n, ...){

  w<-as.matrix(na.omit(x))
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
