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
  cat("\n NARDL model:\n")
  print(object$sels)
  cat("\n model diagnostic tests:\n")
  diagnos=cbind(object$jbtest$statistic[[1]],object$lm2[1],object$arch$statistic[[1]])
  pval=cbind(object$jbtest$p.value[[1]],object$lm2[2],object$arch$p.value[[1]])
  nlag=cbind(0,object$nl,object$nl)
  resdiag=rbind(diagnos,pval,nlag)
  colnames(resdiag)=c(" JB test","LM test","ARCH test")
  rownames(resdiag)=c("Stat","Pvalue","lags")
  print(resdiag)
  cat("==============================================================\n")
  cat(" Short Run Asymmety test\n","W-stat:",object$wldsr[1],"Pvalue:",object$wldsr[2],"\n")
  cat("==============================================================\n")
  pssbounds(case=3,fstat= object$fstat,obs=object$Nobs,k=2)
  cat("==============================================================\n")
  cat("\nLong-run coefficients\n")
  printCoefmat(object$lres,has.Pvalue = TRUE,signif.stars = TRUE)
  cat("==============================================================\n")
  cat(" Long Run Asymmety test\n","W-stat:",object$wldq[1],"Pvalue:",object$wldq[2],"\n")
  cat("==============================================================\n")

}
