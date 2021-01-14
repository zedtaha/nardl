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
 ==============================================================

 NARDL model: 
 Call:
 lm(formula = dy ~ lay + lxp + lxn, na.action = na.exclude)
 
 Residuals:
    Min      1Q  Median      3Q     Max 
 -4.4288 -1.6970 -0.0655  1.1475 10.2080 
 
 Coefficients:
         Estimate Std. Error t value Pr(>|t|)   
 Const    5.30911    3.50852   1.513  0.13790   
 food_1  -0.35328    0.11640  -3.035  0.00417 **
 inf_p   -0.09718    0.15203  -0.639  0.52622   
 inf_p_1  0.50629    0.21717   2.331  0.02473 * 
 inf_p_2 -0.10542    0.20334  -0.518  0.60694   
 inf_p_3  0.10229    0.19521   0.524  0.60312   
 inf_p_4 -0.25546    0.14409  -1.773  0.08368 . 
 inf_n    0.20342    0.12783   1.591  0.11921   
 ---
 Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
 
 Residual standard error: 3.025 on 41 degrees of freedom
   (4 observations deleted due to missingness)
 Multiple R-squared:  0.396,  Adjusted R-squared:  0.2929 
 F-statistic:  3.84 on 7 and 41 DF,  p-value: 0.002689

 
  model diagnostic tests:
            JB test   LM test  ARCH test
 Stat   0.898800930 0.1081374 0.03920182
 Pvalue 0.000502159 0.7977438 0.84304937
 lags   0.000000000 1.0000000 1.00000000
 ==============================================================
  Short Run Asymmety test
  W-stat: 2.392574 Pvalue: 0.3023147 
 ==============================================================
 
  PESARAN, SHIN AND SMITH (2001) COINTEGRATION TEST 
 
  Observations: 49 
  Number of Regressors (k): 2 
  Case: 3 
 
  ------------------------------------------------------ 
  -                       F-test                       - 
  ------------------------------------------------------ 
                  <------- I(0) ------------ I(1) -----> 
  10% critical value       3.333            4.313 
  5% critical value        4.07            5.19 
  1% critical value        5.817            7.303 
  
 
  F-statistic = 3.84006390330111 
  
  ------------------------------------------------------ 
  
  
 ==============================================================
 
 Long-run coefficients
         Estimate Std. Error t value Pr(>|t|)  
 inf_p   -0.27508    0.41952 -0.6557  0.51201  
 inf_p_1  1.43310    0.80035  1.7906  0.07336 .
 inf_p_2 -0.29839    0.57873 -0.5156  0.60613  
 inf_p_3  0.28953    0.54957  0.5268  0.59831  
 inf_p_4 -0.72311    0.41602 -1.7382  0.08218 .
 inf_n    0.57579    0.46801  1.2303  0.21859  
 ---
 Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
 ==============================================================
  Long Run Asymmety test
  W-stat: 19.16983 Pvalue: 6.875809e-05 
 ==============================================================
# Cointegration bounds test
```{r}
reg<-nardl(food~inf,fod,ic="aic",maxlag = 4,graph = TRUE,case=3)
pssbounds(case=reg$case,fstat=reg$fstat,obs=reg$obs,k=reg$k)
```

# License
All code is licensed [![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
