#' ARCH test
#'
#' Computes the Lagrange multiplier test for conditional
#' heteroscedasticity of Engle (1982), as described by
#' Tsay (2005, pp. 101-102).
#'
#' @param x numeric vector, \code{\link{nardl}} or
#' \code{\link{summary.nardl}} object.
#' @param lags positive integer number of lags. This will be ignored
#' when a \code{nardl} or \code{summary.nardl} is passed as \code{x}.
#' @param demean logical: If TRUE, remove the mean before computing
#' the test statistic. This will be ignored when
#' @param ... arguments passed from or to other methods
#'
#' @note This function is modeled after the \code{htest} class functions
#' like \code{\link[stats]{t.test}} and others. This also means that
#' when using the method on a \code{nardl} or \code{summary.nardl}
#' object, the name of the data displayed will always be ''residu''
#'
#' @return a \code{htest} object.
#' @importFrom stats pchisq embed lm
#' @examples
#'
#' reg<-nardl(food~inf,fod,ic="aic",maxlags = TRUE,graph = TRUE,case=3)
#'
#' ArchTest(reg)
#' ArchTest(reg, demean = TRUE)
#'
#' reg.summ <- summary(reg)
#' ArchTest(reg)
#'
#' # Manual calculation
#' x<-reg$selresidu
#' nlag<-reg$np
#' ArchTest(x,lags=nlag)
#'
#' @export
ArchTest <- function(x, ...){
  UseMethod("ArchTest")
}

#' @rdname ArchTest
#' @export
ArchTest.nardl <- function(x, demean = FALSE, ...){

  if(demean){
    residu <- x$selresidu
    nlag <- x$np
    ArchTest.default(residu, lags = nlag, demean = demean)
  } else {
    # extract from nardl object if demean is FALSE
    x$arch
  }
}

#' @rdname ArchTest
#' @export
ArchTest.summary.nardl <- function(x, ...){

  x$archtest

}

#' @rdname ArchTest
#' @export
ArchTest.default <- function (x, lags=12, demean = FALSE , ...)
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
