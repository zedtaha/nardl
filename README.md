---
output: github_document
---
[![](http://www.r-pkg.org/badges/version/nardl)]()
[![](http://cranlogs.r-pkg.org/badges/nardl)](http://cran.rstudio.com/web/packages/nardl/index.html)
# nardl
nardl:An R package to estimate the nonlinear cointegrating autoregressive distributed lag model

In this package, we apply the ordinary least squares method to estimate the cointegrating nonlinear ARDL (NARDL) model in which short and long-run nonlinearities are introduced via positive and negative partial sum decompositions of the explanatory variables.Besides, we provide the CUSUM, CUSUMSQ model stability tests, model selection via aic and bic.

# Example 
 To dwonload the package you can just type install.packages("nardl") in R or Rstudio.
 
```R
# load the package 
 library(nardl)
```
The data used in this example is the Indian yearly data of inflation rate and percentage food import to total import.
```{r}
# formula food~inf or food~inf|I(inf^2)
# data the dataframe
# ic : c("aic","bic") criteria model selection
# maxlag the maximum lag number
# graph TRUE to show stability tests plot
# case case number 3 for (unrestricted intercert, no trend) and 5 (unrestricted intercept, unrestricted trend), 1 2 and 4 not supported

############################################
# Fit the nonlinear cointegrating autoregressive distributed lag model
############################################
# Load data
data(fod)
reg<-nardl(food~inf,fod,ic="aic",maxlag = 4,graph = TRUE,case=3)
summary(reg)

```

# Cointegration bounds test
```{r}
reg<-nardl(food~inf,fod,ic="aic",maxlag = 4,graph = TRUE,case=3)
pssbounds(case=reg$case,fstat=reg$fstat,obs=reg$obs,k=reg$k)
```

# License
All code is licensed [![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
