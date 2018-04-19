---
title: 'nardl:An R package to estimate the nonlinear cointegrating autoregressive distributed lag model'
tags:
    - Autoregressive distributed lag model
    - Nonlinear cointegrating
    - R
    - nardl
authors: 
    - name: Taha Zaghdoudi
      orcid: 0000-0002-4488-5774
      affiliation: 1   
affiliations:
    - name: Faculty of Law, Economics and Management of Jendouba, Tunisia
      index: 1    
date: 12 Mars 2018
bibliography: paper.bib
---

  
# Summary #
In the nardl package [@zaghdoudi], we apply the ordinary least squares method to estimate the cointegrating nonlinear ARDL (NARDL) 
model developed by [@shin2014modelling] in which short and long-run nonlinearities are introduced via positive and negative partial 
sum decompositions of the explanatory variables.Besides, we provide the CUSUM, CUSUMSQ model stability tests, model selection via aic, 
bic and rsqaured criteria and the dynamic multipliers plot.

# Example #

In this example we examine the impact of both long and short-run asymmetries effect of inflation on food price in India. we  specify the following asymmetric long-run equation of inflation:

$$Food_{t}= \alpha_{0}+\alpha_{1}infp_{t}+\alpha_{é}infn_{t}+\varepsilon_{t}$$

Where $Food_{t}$ refers to the food price and$\alpha=(\alpha_{0},\alpha_{1},\alpha_{2}$ is a $infp_{t}$ and$infn_{t}$ are partial sums of positive and negative changes in $inf_{t}$:


	
$$infp_{t}=\sum_{i=1}^{t}\Delta infp_{t}=\sum_{i=1}^{t} max(\Delta inf_{t},0)$$ 
		
	
and
	
$$infn_{t}=\sum_{i=1}^{t}\Delta infn_{t}=\sum_{i=1}^{t} min(\Delta inf_{t},0)$$ 


```{r}
# Loading package
library(nardl)
# Importing data from Excel sheet
data(fod)
# Estimate the NARDL model
# food: the dependent variable
# inf : the decomposed dependent variable in infp and infn
# (p,q)=(4,4): max lags for dependent and independent variables
# ic=c(aic,bic, R2): selection model criteria
# if graph=TRUE this make a stability tests plot 
reg<-nardl(food~inf,p=4,q=4,fod,ic="aic",maxlags = FALSE,graph = TRUE)
summary(reg)
```

# Cointegration bounds test
```{r}
pssbounds(case=reg$case,fstat=reg$fstat,obs=reg$obs,k=reg$k)
```
# Dynamic multipliers plot
```{r}
############################
# Dynamic multipliers plot
############################
plotmplier(reg,reg$np,1,10)
```


nardl is GPL-3 licensed and can be retrieved from GitHub [Repository](https://github.com/cran/nardl)

[nardl CRAN](https://cran.r-project.org/web/packages/nardl/index.html)

[Archive](https://zenodo.org/record/1221245#.WtidtS7wbIU)

# References #

