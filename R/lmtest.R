#' LM test for serial correlation
#'
#' This function performs a LM test for serial correlation. The function
#' \code{bp2()} is the workhorse, but the function \code{LMTest()} is
#' preferred as that function does the proper dispatching of methods.
#'
#' @param x a \code{\link{nardl}} or \code{\link{summary.nardl}} object.
#' @param object fitted lm model
#' @param nlags positive integer number of lags
#' @param fill starting values for the lagged residuals in the
#' auxiliary regression. By default 0.
#' @param type Fisher or Chisquare statistics
#' @param ... arguments passed from and to other methods
#'
#' @return a matrix with 1 row and 3 named columns, containing the
#' value for LM, the p value and the number of lags used.
#'
#' @importFrom stats lm pchisq
#' @importFrom gtools na.replace
#' @examples
#'
#' reg<-nardl(food~inf,fod,ic="aic",maxlags = TRUE,graph = TRUE,case=3)
#' LMTest(reg)
#' lm2<-bp2(reg$fit,reg$np,fill=0,type="F")
#'
#'@export
LMTest <- function(x, ...){
  UseMethod("LMTest")
}

#' @rdname LMTest
#' @export
LMTest.nardl <- function(x,
                         fill = NULL,
                         type = c("F", "Chi2"),
                         ...){

  if(missing(fill) && missing(type)){
    x$lm2
  } else {
    fit <- x$fit
    nlag <- x$np
    bp2(fit, nlag, fill = fill, type = type)
  }
}

#' @rdname LMTest
#' @export
LMTest.summary.nardl <- function(x,
                                 ...){
  matrix(
    c(x$statistic, x$p.value, x$np),
    nrow = 1,
    dimnames = list(NULL,
                    c("LM","P-value","lags"))
  )

}
#' @rdname LMTest
#' @aliases bp2
#' @export
bp2<-function(object,nlags,fill=NULL,type=c("F","Chi2")){

  # Needed to avoid errors about condition having length > 1
  type <- match.arg(type)
  if(is.null(fill))
    stop("You need to specify fill.")

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
