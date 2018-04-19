#'Nonlinear ARDL function
#'
#'@param formula food~inf or food~inf|I(inf^2)
#'@param data the dataframe
#'@param p  lags of dependent variable
#'@param q  lags of independent variables
#'@param ic : c("aic","bic","ll","R2") criteria model selection
#'@param maxlags if TRUE auto lags selection
#'@param graph TRUE to show stability tests plot
#'@param case case number 3 for (unrestricted intercert, no trend) and 5 (unrestricted intercept, unrestricted trend), 1 2 and 4 not supported
#'@importFrom stats lm AIC BIC pchisq as.formula model.frame model.matrix model.response na.omit sd update vcov embed resid coef logLik nobs pf pnorm df.residual formula
#'@importFrom graphics abline legend lines par plot
#'@importFrom strucchange recresid
#'@importFrom tseries jarque.bera.test
#'@importFrom gtools na.replace
#'@importFrom utils head
#'@import Formula
#'@examples
#'
#'############################################
#'# Fit the nonlinear cointegrating autoregressive distributed lag model
#'############################################
#'# Load data
#'data(fod)
#'############################################
#'#example 1: nardl with fixed p and q lags
#'############################################
#'reg<-nardl(food~inf,p=4,q=4,fod,ic="aic",maxlags = FALSE,graph = FALSE,case=3)
#'summary(reg)
#'
#'############################################
#'# example 2:auto selected lags (maxlags=TRUE)
#'############################################
#'reg<-nardl(food~inf,fod,ic="aic",maxlags = TRUE,graph = FALSE,case=3)
#'summary(reg)
#'
#'############################################
#'# example 3: Cusum and CusumQ plot (graph=TRUE)
#'############################################
#'reg<-nardl(food~inf,fod,ic="aic",maxlags = TRUE,graph = TRUE,case=3)
#'
#'@export
nardl<-function(formula,data,p=NULL,q=NULL,ic=c("aic","bic","ll","R2"),
                maxlags=TRUE,graph=FALSE,case= 3){

  f<-formula
  a<-unlist(strsplit(as.character(f),"[|]"))
  #core<-paste(un[[2]],"~",un[[3]],"+",un[[4]])
  if(length(a)==4){
    lhs   <- a[[2]]
    core  <- a[[3]]
    suffix <- a[[4]]
    fm <-as.formula(paste(lhs,"~",core,"|",suffix))
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0)
    mf <- mf[c(1, m)]
    ffm <- Formula::Formula(fm)
    mf[[1]] <- as.name("model.frame")
    mf$formula <- ffm
    fffm<-update(ffm, ~ . -1|~ . -1)
    mf <- model.frame(formula=fffm, data=data)
    x <- model.matrix(fffm, data = mf, rhs = 1)
    if(ncol(x)>=2) stop("nardl package accept only one decomposed variable")
    h <- model.matrix(fffm, data = mf, rhs = 2)
    y<-model.response(mf)
    y<-as.matrix(y)
    k<-ncol(x)+ncol(h)
    #colnames(h)<-suffix
    colnames(y)<-lhs
    dy<-diff(y)
    colnames(dy)<-paste("D",lhs,sep=".")
    dx<-diff(x)
    colnames(dx)<-paste("D",colnames(x),sep=".")
    dh<-diff(h)
    colnames(dh)<-paste("D",colnames(h),sep=".")
    n<-nrow(dx)
    cl<-ncol(dx)
    pos<-dx[,1:cl]>=0
    dxp<-as.matrix(as.numeric(pos)*dx[,1:cl])
    colnames(dxp)<-paste(colnames(dx),"p",sep="_")
    dxn<-as.matrix((1-as.numeric(pos))*dx[,1:cl])
    colnames(dxn)<-paste(colnames(dx),"n",sep="_")
      xp<-apply(dxp,2,cumsum)
      colnames(xp)<-paste(colnames(x),"p",sep="_")
      xn<-apply(dxn,2,cumsum)
      colnames(xn)<-paste(colnames(x),"n",sep="_")

    lagy<-lagm(as.matrix(y),1)
    lagh<-lagm(as.matrix(h),1)
    lxp<-lagm(as.matrix(xp),1)
    lxn<-lagm(as.matrix(xn),1)
    #fit<-lm(dy~na.omit(lagy)+lxp+lxn)

    #lag order selection
    R2 = NULL
    AK = NULL
    SC = NULL
    ll= NULL
    if(maxlags==TRUE){
      ordmax<-floor(nrow(dy)^(1/3))
    }
    if(is.null(p) && is.null(q)){

      ordmax<-floor(nrow(dy)^(1/3))

    }else{
      ordmax<-max(p,q)
      }

    #
    for (i in 1:ordmax){
      #lagmat = cbind(lagmat[-i,],x[(1):(l1-i)]) # lagged matrix
      # armod <- lm(x[(i+1):l1]~lagmat)
      ldy<-lagm(dy,i)
      ldh<-lagm(as.matrix(dh),i)
      ldxp<-lagm(as.matrix(dxp),i)
      ldxn<-lagm(as.matrix(dxn),i)
      lagy<-na.omit(lagy)
      lagh<-na.omit(lagh)
      fit<-lm(dy~lagy+lxp+lxn+ldy+ldxp+ldxn+lagh+ldh)
      R2[i] = summary(fit)$adj.r.squared
      AK[i] = AIC(fit)
      SC[i] = BIC(fit)
      ll[i] =as.numeric(logLik(fit))
    }
    np<-0
    #ores<-c(which.max(R2), which.min(AK),which.min(SC) )
    if(ic=="aic") np<-which.min(AK)
    if(ic=="bic") np<-which.min(SC)
    if(ic=="R2") np<-which.max(R2)
    if(ic=="ll") np<-which.max(ll)
    #nnp<-ores[duplicated(ores)]
    #nnp<-nnp[1]
    #np<-3
    ldy<-lagm(dy,np)
    ldh<-lagm(as.matrix(dh),np)
    ldxp<-lagm(as.matrix(dxp),np)
    ldxn<-lagm(as.matrix(dxn),np)
    lagy<-na.omit(lagy)
    lagh<-na.omit(lagh)
    lhnames<-colnames(dy)
    if(case==5){
      trend<-seq_along(dy)
      rhnams<-c("Const",colnames(lagy),colnames(lxp),colnames(lxn),
                colnames(ldy),colnames(ldxp),colnames(ldxn),
                colnames(lagh),colnames(ldh),"trend")
      fit<-lm(dy~lagy+lxp+lxn+ldy+ldxp+ldxn+lagh+ldh+trend)
    }
    if(case==3){
      rhnams<-c("Const",colnames(lagy),colnames(lxp),colnames(lxn),
                colnames(ldy),colnames(ldxp),colnames(ldxn),
                colnames(lagh),colnames(ldh))
      fit<-lm(dy~lagy+lxp+lxn+ldy+ldxp+ldxn+lagh+ldh)
    }
    names(fit$coefficients)<-rhnams
    fitcoef<-summary(fit)

  }else{
    lhs   <- a[[2]]
    core  <- a[[3]]
    fm <- paste(lhs,"~",core)
    ffm<-as.formula(fm)
    fffm<-update(ffm, ~ . -1)
    mf <- model.frame(formula=fffm, data=data)
    x <- model.matrix(attr(mf, "terms"), data=mf)
    k<-ncol(x)
    if(ncol(x)>=2) stop("nardl package accept only one decomposed variable")
    y <- model.response(mf)
    y <- model.response(mf)
    y<-as.matrix(y)
    colnames(y)<-lhs
    dy<-diff(y)
    colnames(dy)<-paste("D",lhs,sep=".")
    dx<-diff(x)
    colnames(dx)<-paste("D",core,sep=".")
    n<-nrow(dx)
    cl<-ncol(dx)
    pos<-dx[,1:cl]>=0
    dxp<-as.matrix(as.numeric(pos)*dx[,1:cl])
    colnames(dxp)<-paste(colnames(dx),"p",sep="_")
    dxn<-as.matrix((1-as.numeric(pos))*dx[,1:cl])
    colnames(dxn)<-paste(colnames(dx),"n",sep="_")
      xp<-apply(dxp,2,cumsum)
      colnames(xp)<-paste(colnames(x),"p",sep="_")
      xn<-apply(dxn,2,cumsum)
      colnames(xn)<-paste(colnames(x),"n",sep="_")
    lagy<-lagm(as.matrix(y),1)
    lxp<-lagm(as.matrix(xp),1)
    lxn<-lagm(as.matrix(xn),1)
    #lag order selection
    R2 = NULL
    AK = NULL
    SC = NULL
    ll= NULL
    if(maxlags){
      ordmax<-floor(nrow(dy)^(1/3))
    }
    if(is.null(p) && is.null(q)){

      ordmax<-floor(nrow(dy)^(1/3))

    }else{
      ordmax<-max(p,q)
    }

    for (i in 1:ordmax){
      ldy<-lagm(dy,i)
      ldxp<-lagm(as.matrix(dxp),i)
      ldxn<-lagm(as.matrix(dxn),i)
      lagy<-na.omit(lagy)
      fit<-lm(dy~lagy+lxp+lxn+ldy+ldxp+ldxn)
      R2[i] = summary(fit)$adj.r.squared
      AK[i] = AIC(fit)
      SC[i] = BIC(fit)
      ll[i] =as.numeric(logLik(fit))
    }
    np<-0
    if(ic=="aic") np<-which.min(AK)
    if(ic=="bic") np<-which.min(SC)
    if(ic=="R2") np<-which.max(R2)
    if(ic=="ll") np<-which.max(ll)
    ldy<-lagm(dy,np)
    ldxp<-lagm(as.matrix(dxp),np)
    ldxn<-lagm(as.matrix(dxn),np)
    lagy<-na.omit(lagy)


    # flm<-paste(lhnames,"~",paste(rhnams,collapse = " + "))
    if(case==5){
      trend<-seq_along(dy)
      rhnams<-c("Const",colnames(lagy),colnames(lxp),colnames(lxn),
                colnames(ldy),colnames(ldxp),colnames(ldxn),"trend")
      fit<-lm(dy~lagy+lxp+lxn+ldy+ldxp+ldxn+trend)
    }
    if(case==3){
      rhnams<-c("Const",colnames(lagy),colnames(lxp),colnames(lxn),
                colnames(ldy),colnames(ldxp),colnames(ldxn))
      fit<-lm(dy~lagy+lxp+lxn+ldy+ldxp+ldxn)
    }
    names(fit$coefficients)<-rhnams
    fitcoef<-summary(fit)
  }

  #long run coeff
  sel<-summary(fit)
  residu<-sel$residuals
  #jarque berra test for residual normallity
  jbtest<-tseries::jarque.bera.test(residu)
  #lmtest and ARCH tests
  #lmtest<-lmtest::bgtest(fit,order = np, type = "F", fill=NA)
  #lm2<-bptest(fit,np)
  lm2<-bp2(fit,np,fill=0,type="F")
  arch<-ArchTest(residu,np)

  coeff<-sel$coefficients
  nlvars<-length(coeff[,1])
  lvars<-coeff[3:nlvars,1]
  seldata<-data.matrix(coeff)
  coof<- -lvars/coeff[[2]]
  cof<- matrix(coof, length(lvars),1)
  #-----SE by delta method
  vc<-vcov(sel)
  vcc<-vc[2:nrow(vc),2:ncol(vc)]
  nsel<-length(sel$coefficients[,1])
  fb1<-lvars/coeff[[2]]^2
  fb2<-(-coeff[[2]])*diag(nrow(as.matrix(fb1)))
  fb<-cbind(as.matrix(fb1),fb2)
  lrse<-sqrt(diag(fb%*%vcc%*%t(fb)))
  lrt<-coof/lrse
  #wald test
 # lrse2<-lrse^2
  rcap<-matrix(0,1,2)
 rcap[1,1]<-1
  rcap[1,2]<--1
  ###############
  lrpv<-2 * pnorm(-abs(lrt))
  lres<-cbind(cof,lrse,lrt,lrpv)
  #lres
  rownames(lres)<-names(lvars)
  colnames(lres)<-c("Estimate", "Std. Error", "t value", "Pr(>|t|)" )
  #---------------------------------------------------
  #cointegration test
  l<-c(paste(names(coeff[-1,1]),"=0"))
  coin<-linearHypothesis(fit,l,test="F")
  fstat<-coin$F[2]
  ll2<-c(paste(names(coof[1]),"=",names(coof[2]) ))
  asym<-linearHypothesis(fit,ll2,test="F")
  tasym<-asym$F
  pasym<-asym$`Pr(>F)`
  #long run wald test
  #get recursive residuals
  rece<-strucchange::recresid(fit)
 # plot cumsum and cumsumq
  if(graph==TRUE){
    bbst<-sel$coefficients[,1]
    kt<-length(bbst[-1])
    n<-length(residu)
    cusum(rece,kt,n)
    cumsq(rece,kt,n)
  }

  obs<-nobs(fit)

 out<-list(fstat=fstat, fit=fit,fitcoef=fitcoef, sel=sel, cof=cof,coof=coof,
       k=k,case=case,lvars=lvars,selresidu=residu,
       tasym=tasym,pasym=pasym,jbtest=jbtest,arch=arch,
       np=np,rece=rece,obs=obs,AK=AK,R2=R2,SC=SC,ll=ll,vcc=vcc,fb=fb,
       fb1=fb1,fb2=fb2,lrse=lrse,lrt=lrt,lrpv=lrpv,lres=lres,lm2=lm2)

 class(out) <- "nardl"
 # Return results.
 return(out)



}




