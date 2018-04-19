#' Summary of a nardl model
#'
#' \code{summary} method for a \code{\link{nardl}} model.
#'
#' @param object is the object of the function
#' @param ... not used
#'
#' @return an object of the S3 class \code{summary.nardl} with the
#' following components:
#'
#' @importFrom stats printCoefmat
#' @rdname summary.nardl
#' @name summary.nardl
#' @export
summary.nardl<-function(object,...)
{

  cat("==============================================================\n")
  cat("\n Lag and lead selection:\n")
  icres<-cbind(object$AK,object$SC,object$ll,object$R2)
  colnames(icres)<-c("AIC","BIC","ll","R2")
  print(icres)
  cat("\n NARDL model:\n")
  print(object$sel)
  cat("---------------------------------\n")
  cat("\n model diagnostic tests:\n----------\n")
  cat(" JB test:\n","JB:",object$jbtest$statistic[[1]],"Pvalue",object$jbtest$p.value[[1]],"\n----------\n")
  cat(" LM test for serial correlation:\n","LM(",object$np,"):",object$lm2[1],"Pvalue",object$lm2[2],"\n----------\n")
  cat(" ARCH test:\n","ARCH(",object$np,"):",object$arch$statistic[[1]],"Pvalue",object$arch$p.value[[1]],"\n----------\n")
  cat("==============================================================\n")
  cat("\nLong-run coefficients\n")
  printCoefmat(object$lres,has.Pvalue = TRUE,signif.stars = TRUE)
  cat("==============================================================\n")
  cat(" Long Run Asymmety test\n","F-stat:",object$tasym[2],"Pvalue:",object$pasym[2],"\n")
  cat("==============================================================\n")

}
