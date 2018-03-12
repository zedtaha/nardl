---
output: github_document
---

[![](http://www.r-pkg.org/badges/version/nardl)]()
[![](http://cranlogs.r-pkg.org/badges/nardl)](http://cran.rstudio.com/web/packages/nardl/index.html)
# nardl
nardl:An R package to estimate the nonlinear cointegrating autoregressive distributed lag model

In this package, we apply the ordinary least squares method to estimate the cointegrating nonlinear ARDL (NARDL) model in which short and long-run nonlinearities are introduced via positive and negative partial sum decompositions of the explanatory variables.Besides, we provide the CUSUM, CUSUMSQ model stability tests, model selection via aic, bic and rsqaured criteria and the dynamic multipliers plot.

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
# p  lags of dependent variable
# q  lags of independent variables
# ic : c("aic","bic","ll","R2") criteria model selection
# maxlags if TRUE auto lags selection
# graph TRUE to show stability tests plot
# case case number 3 for (unrestricted intercert, no trend) and 5 (unrestricted intercept, unrestricted trend), 1 2 and 4 not supported

############################################
# Fit the nonlinear cointegrating autoregressive distributed lag model
############################################
# Load data
data(fod)
############################################
#example 1: nardl with fixed p and q lags
############################################
reg<-nardl(food~inf,p=4,q=4,fod,ic="aic",maxlags = FALSE,graph = FALSE,case=3)
summary(reg)

############################################
# example 2:auto selected lags (maxlags=TRUE)
############################################
reg<-nardl(food~inf,fod,ic="aic",maxlags = TRUE,graph = FALSE,case=3)
summary(reg)

# ############################################
# example 3: Cusum and CusumQ plot (graph=TRUE)
# ############################################
reg<-nardl(food~inf,fod,ic="aic",maxlags = TRUE,graph = TRUE,case=3)

```
# Dynamic multipliers plot
```{r}
############################
# Dynamic multipliers plot
############################
# Load data
data(fod)
reg<-nardl(food~inf,p=4,q=4,fod,ic="aic",maxlags = FALSE,graph = TRUE,case=3)
plotmplier(reg,reg$np,1,10)
```
# License
All code is licensed [![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
