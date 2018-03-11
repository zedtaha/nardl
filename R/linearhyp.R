makeHypothesis <- function(cnames, hypothesis, rhs = NULL){
  parseTerms <- function(terms){
    component <- gsub("^[-\\ 0-9\\.]+", "", terms)
    component <- gsub(" ", "", component, fixed=TRUE)
    component
  }
  stripchars <- function(x) {
    x <- gsub("\\n", " ", x)
    x <- gsub("\\t", " ", x)
    x <- gsub(" ", "", x, fixed = TRUE)
    x <- gsub("*", "", x, fixed = TRUE)
    x <- gsub("-", "+-", x, fixed = TRUE)
    x <- strsplit(x, "+", fixed = TRUE)[[1]]
    x <- x[x!=""]
    x
  }
  char2num <- function(x) {
    x[x == ""] <- "1"
    x[x == "-"] <- "-1"
    as.numeric(x)
  }
  constants <- function(x, y) {
    with.coef <- unique(unlist(sapply(y,
                                      function(z) which(z == parseTerms(x)))))
    if (length(with.coef) > 0) x <- x[-with.coef]
    x <- if (is.null(x)) 0 else sum(as.numeric(x))
    if (any(is.na(x)))
      stop('The hypothesis "', hypothesis,
           '" is not well formed: contains bad coefficient/variable names.')
    x
  }
  coefvector <- function(x, y) {
    rv <- gsub(" ", "", x, fixed=TRUE) ==
      parseTerms(y)
    if (!any(rv)) return(0)
    if (sum(rv) > 1) stop('The hypothesis "', hypothesis,
                          '" is not well formed.')
    rv <- sum(char2num(unlist(strsplit(y[rv], x, fixed=TRUE))))
    if (is.na(rv))
      stop('The hypothesis "', hypothesis,
           '" is not well formed: contains non-numeric coefficients.')
    rv
  }

  if (!is.null(rhs)) rhs <- rep(rhs, length.out = length(hypothesis))
  if (length(hypothesis) > 1)
    return(rbind(Recall(cnames, hypothesis[1], rhs[1]),
                 Recall(cnames, hypothesis[-1], rhs[-1])))

  cnames_symb <- sapply(c("@", "#", "~"), function(x) length(grep(x, cnames)) < 1)

  if(any(cnames_symb)) {
    cnames_symb <- head(c("@", "#", "~")[cnames_symb], 1)
    cnames_symb <- paste(cnames_symb, seq_along(cnames), cnames_symb, sep = "")
    hypothesis_symb <- hypothesis
    for(i in order(nchar(cnames), decreasing = TRUE))
      hypothesis_symb <- gsub(cnames[i], cnames_symb[i], hypothesis_symb, fixed = TRUE)
  } else {
    stop('The hypothesis "', hypothesis,
         '" is not well formed: contains non-standard coefficient names.')
  }

  lhs <- strsplit(hypothesis_symb, "=", fixed=TRUE)[[1]]
  if (is.null(rhs)) {
    if (length(lhs) < 2) rhs <- "0"
    else if (length(lhs) == 2) {
      rhs <- lhs[2]
      lhs <- lhs[1]
    }
    else stop('The hypothesis "', hypothesis,
              '" is not well formed: contains more than one = sign.')
  }
  else {
    if (length(lhs) < 2) as.character(rhs)
    else stop('The hypothesis "', hypothesis,
              '" is not well formed: contains a = sign although rhs was specified.')
  }
  lhs <- stripchars(lhs)
  rhs <- stripchars(rhs)
  rval <- sapply(cnames_symb, coefvector, y = lhs) - sapply(cnames_symb, coefvector, y = rhs)
  rval <- c(rval, constants(rhs, cnames_symb) - constants(lhs, cnames_symb))
  names(rval) <- c(cnames, "*rhs*")
  rval
}

printHypothesis <- function(L, rhs, cnames){
  hyp <- rep("", nrow(L))
  for (i in 1:nrow(L)){
    sel <- L[i,] != 0
    h <- L[i, sel]
    h <- ifelse(h < 0, as.character(h), paste("+", h, sep=""))
    nms <- cnames[sel]
    h <- paste(h, nms)
    h <- gsub("-", " - ", h)
    h <- gsub("+", "  + ", h, fixed=TRUE)
    h <- paste(h, collapse="")
    h <- gsub("  ", " ", h, fixed=TRUE)
    h <- sub("^\\ \\+", "", h)
    h <- sub("^\\ ", "", h)
    h <- sub("^-\\ ", "-", h)
    h <- paste(" ", h, sep="")
    h <- paste(h, "=", rhs[i])
    h <- gsub(" 1([^[:alnum:]_.]+)[ *]*", "",
              gsub("-1([^[:alnum:]_.]+)[ *]*", "-",
                   gsub("- +1 +", "-1 ", h)))
    h <- sub("Intercept)", "(Intercept)", h)
    h <- gsub("-", " - ", h)
    h <- gsub("+", "  + ", h, fixed=TRUE)
    h <- gsub("  ", " ", h, fixed=TRUE)
    h <- sub("^ *", "", h)
    hyp[i] <- h
  }
  hyp
}



linearHypothesis <- function(model, hypothesis.matrix, rhs=NULL,
                                     test=c("Chisq", "F"), vcov.=NULL, singular.ok=FALSE, verbose=FALSE,
                                     coef. = coef(model), ...){
  df <- df.residual(model)
  if (is.null(df)) df <- Inf ## if no residual df available
  if (df == 0) stop("residual df = 0")
  V <- if (is.null(vcov.)) vcov(model, complete=FALSE)
  else if (is.function(vcov.)) vcov.(model) else vcov.
  b <- coef.
  if (any(aliased <- is.na(b)) && !singular.ok)
    stop("there are aliased coefficients in the model")
  b <- b[!aliased]
  if (is.null(b)) stop(paste("there is no coef() method for models of class",
                             paste(class(model), collapse=", ")))
  if (is.character(hypothesis.matrix)) {
    L <- makeHypothesis(names(b), hypothesis.matrix, rhs)
    if (is.null(dim(L))) L <- t(L)
    rhs <- L[, NCOL(L)]
    L <- L[, -NCOL(L), drop = FALSE]
    rownames(L) <- hypothesis.matrix
  }
  else {
    L <- if (is.null(dim(hypothesis.matrix))) t(hypothesis.matrix)
    else hypothesis.matrix
    if (is.null(rhs)) rhs <- rep(0, nrow(L))
  }
  q <- NROW(L)
  value.hyp <- L %*% b - rhs
  vcov.hyp <- L %*% V %*% t(L)
  if (verbose){
    cat("\nHypothesis matrix:\n")
    print(L)
    cat("\nRight-hand-side vector:\n")
    print(rhs)
    cat("\nEstimated linear function (hypothesis.matrix %*% coef - rhs)\n")
    print(drop(value.hyp))
    cat("\n")
    if (length(vcov.hyp) == 1) cat("\nEstimated variance of linear function\n")
    else cat("\nEstimated variance/covariance matrix for linear function\n")
    print(drop(vcov.hyp))
    cat("\n")
  }
  SSH <- as.vector(t(value.hyp) %*% solve(vcov.hyp) %*% value.hyp)
  test <- match.arg(test)
  if (!(is.finite(df) && df > 0)) test <- "Chisq"
  name <- try(formula(model), silent = TRUE)
  if (inherits(name, "try-error")) name <- substitute(model)
  title <- "Linear hypothesis test\n\nHypothesis:"
  topnote <- paste("Model 1: restricted model","\n", "Model 2: ",
                   paste(deparse(name), collapse = "\n"), sep = "")
  note <- if (is.null(vcov.)) ""
  else "\nNote: Coefficient covariance matrix supplied.\n"
  rval <- matrix(rep(NA, 8), ncol = 4)
  colnames(rval) <- c("Res.Df", "Df", test, paste("Pr(>", test, ")", sep = ""))
  rownames(rval) <- 1:2
  rval[,1] <- c(df+q, df)
  if (test == "F") {
    f <- SSH/q
    p <- pf(f, q, df, lower.tail = FALSE)
    rval[2, 2:4] <- c(q, f, p)
  }
  else {
    p <- pchisq(SSH, q, lower.tail = FALSE)
    rval[2, 2:4] <- c(q, SSH, p)
  }
  if (!(is.finite(df) && df > 0)) rval <- rval[,-1]
  result <- structure(as.data.frame(rval),
                      heading = c(title, printHypothesis(L, rhs, names(b)), "", topnote, note),
                      class = c("anova", "data.frame"))
  attr(result, "value") <- value.hyp
  attr(result, "vcov") <- vcov.hyp
  result
}
