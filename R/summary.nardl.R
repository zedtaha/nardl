#' Summary of a nardl model
#'
#' \code{summary} method for a \code{\link{nardl}} model.
#'
#' @param object an object of class \code{\link{nardl}}
#' @param ... not used
#'
#' @return an object of the S3 class \code{summary.nardl} with the
#' following components:
#'
#'  \itemize{
#'    \item icres: a matrix with the AIC, BIC, LL and R2 values
#'    \item sel: a \code{\link[stats]{summary.lm}} object with the actual fit
#'    \item np: the number of lags
#'    \item jbtest: a list with the results of the JB test
#'    \item lmtest: a list with the statistic nd p value related to
#'    the LM test.
#'    \item archtest: a list with the results of the arch test
#'    \item lres: a matrix with the long run coefficients and tests
#'    \item asymtest: a list with the statistic and p value of the
#'    long run asymmetry test.
#'
#'  }
#'
#' @examples
#' data(fod)
#' reg<-nardl(food~inf,fod,ic="aic",maxlags = TRUE,graph = FALSE,case=3)
#'
#' summ.reg <- summary(reg)
#' summ.reg
#'
#' @importFrom stats printCoefmat
#' @rdname summary.nardl
#' @name summary.nardl
#' @export
summary.nardl<-function(object,...)
{

  icres<-cbind(object$AK,object$SC,object$ll,object$R2)
  colnames(icres)<-c("AIC","BIC","ll","R2")

  out <- list(
    icres = icres,
    sel = object$sel,
    jbtest = object$jbtest,
    np = object$np,
    lmtest = list(statistic = object$lm2[1],
                  p.value = object$lm2[2]),
    archtest = object$arch,
    lres = object$lres,
    asymtest = list(statistic = object$tasym[2],
                    p.value = object$pasym[2])
  )
  class(out) <- "summary.nardl"
  return(out)

}

# print method for summary.nardl class
#' @export
print.summary.nardl<-function(x,...)
{

  cat("==============================================================\n")
  cat("\n Lag and lead selection:\n")

  print(x$icres)
  cat("\n NARDL model:\n")
  print(x$sel)
  cat("---------------------------------\n")
  cat("\n model diagnostic tests:\n----------\n")
  cat(" JB test:\n","JB:",x$jbtest$statistic[[1]],"Pvalue",x$jbtest$p.value[[1]],"\n----------\n")
  cat(" LM test for serial correlation:\n","LM(",x$np,"):",
      x$lmtest$statistic,
      "Pvalue",x$lmtest$p.value,"\n----------\n")
  cat(" ARCH test:\n","ARCH(",x$np,"):",x$arch$statistic[[1]],"Pvalue",x$arch$p.value[[1]],"\n----------\n")
  cat("==============================================================\n")
  cat("\nLong-run coefficients\n")
  printCoefmat(x$lres,has.Pvalue = TRUE,signif.stars = TRUE)
  cat("==============================================================\n")
  cat(" Long Run Asymmety test\n","F-stat:",x$asymtest$statistic,
      "Pvalue:",x$asymtest$p.value,"\n")
  cat("==============================================================\n")

}
