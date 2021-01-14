#'Nonlinear ARDL function
#'
#'@param formula food~inf or food~inf|I(inf^2)
#'@param data the dataframe
#'@param ic : c("aic","bic") criteria model selection
#'@param maxlag maximum lag number
#'@param graph TRUE to show stability tests plot
#'@param case case number 3 for (unrestricted intercert, no trend) and 5 (unrestricted intercept, unrestricted trend), 1 2 and 4 not supported
#'@importFrom stats lm AIC BIC pchisq as.formula model.frame model.matrix model.response na.omit sd update vcov embed resid coef logLik nobs pf pnorm df.residual formula na.exclude shapiro.test
#'@importFrom strucchange recresid
#'@importFrom tseries jarque.bera.test
#'@importFrom gtools na.replace
#'@importFrom utils head
#'@importFrom car linearHypothesis
#'@import Formula
#'@examples
#'
#'############################################
#'# Fit the nonlinear cointegrating autoregressive distributed lag model
#'############################################
#'# Load data
#'data(fod)
#'############################################
#'# example 1:auto selected lags (maxlags=TRUE)
#'############################################
#'reg<-nardl(food~inf,fod,ic="aic",maxlag = 4,graph = FALSE,case=3)
#'summary(reg)
#'
#'############################################
#'# example 2: Cusum and CusumQ plot (graph=TRUE)
#'############################################
#'reg<-nardl(food~inf,fod,ic="aic",maxlag = 4,graph = TRUE,case=3)
#'
#'@export
nardl<-function(formula,data,ic=c("aic","bic"),maxlag=4,graph=FALSE,case= 3){
 #♣ formula=ly~lco2
  #data=china
  #ic="aic"
  #maxlag=4
 # case=3
  f<-formula
  a<-unlist(strsplit(as.character(f),"[|]"))
  if(length(a)==4){
    lhs   <- a[[2]]
    core  <- a[[3]]
    suffix <- a[[4]]
    fm <-as.formula(paste(lhs,"~",core,"-1","|",suffix,"-1"))
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0)
    mf <- mf[c(1, m)]
    ffm <- Formula::Formula(fm)
    mf[[1]] <- as.name("model.frame")
    mf$formula <- ffm
    mf <- model.frame(formula=ffm, data=data)
    x <- model.matrix(ffm, data = mf, rhs = 1)
    if(ncol(x)>=2) stop("nardl package accept only one decomposed variable")
    h <- model.matrix(ffm, data = mf, rhs = 2)
   # if(ncol(h)>=2) stop("nardl package accept only 2 exogenous variables")
    y<-model.response(mf)
    y<-as.matrix(y)
    k<-ncol(x)+ncol(h)
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

    #*************************************************
    # grid lag selection
    #************************************************

    l <- rep(list(0:maxlag),(k+2))
    grid0=expand.grid(l)
    b=which(grid0[1:nrow(grid0),]==0)
    grid=grid0[-b,]


    #************************************************************
    if(case==5){
      trend<-seq_along(dy)
      if(ncol(h)==2){

        ss=lapply(1:nrow(grid),function(x) lm(dy~nlagm1(y,grid[x,1])[-1,]+xn+nlagm1(xn,grid[x,3])+
                                                xp+nlagm1(xp,grid[x,2])
                                              +slag(h,grid[x,4:ncol(grid)])[-1,]+trend))
        if(ic=="aic"){
          ak=sapply(1:nrow(grid),function(x) AICC(ss[[x]]))
        }
        if(ic=="bic"){
          ak=sapply(1:nrow(grid),function(x) BICC(ss[[x]]))
        }
        cii=as.matrix(ak)
        np=which.min(cii)
        pg=grid[np,]
        lay1=nlagm1(y,pg[,1])
        lay=as.matrix(nlagm1(y,pg[,1])[-1,])
        colnames(lay)=colnames(lay1)
        lxp=nlagm(xp,pg[,2])
        lxn=nlagm(xn,pg[,3])
        lagh1=slag(h,pg[,4:ncol(pg)])[-1,,drop=FALSE]
        rhnams<-c("Const",colnames(lay),colnames(lxp),colnames(lxn),colnames(lagh1),"trend")
        fits= lm(dy~lay+lxp+lxn+lagh1+trend,na.action =na.exclude )
      }else{
        ss=lapply(1:nrow(grid),function(x) lm(dy~nlagm1(y,grid[x,1])[-1,]+xn+nlagm1(xn,grid[x,3])+
                                                    xp+nlagm1(xp,grid[x,2])+nlagm(h,grid[x,4])[-1,]+trend))


      if(ic=="aic"){
        ak=sapply(1:nrow(grid),function(x) AIC(ss[[x]]))
      }
      if(ic=="bic"){
        ak=sapply(1:nrow(grid),function(x) BIC(ss[[x]]))
      }
      cii=as.matrix(ak)
      np=which.min(cii)
      pg=grid[np,]
      lay1=nlagm1(y,pg[,1])
      lay=as.matrix(nlagm1(y,pg[,1])[-1,])
      colnames(lay)=colnames(lay1)
      lxp=nlagm(xp,pg[,2])
      lxn=nlagm(xn,pg[,3])
      lagh1=as.matrix(nlagm(h,pg[,4])[-1,,drop=FALSE])
      rhnams<-c("Const",colnames(lay),colnames(lxp),colnames(lxn),colnames(lagh1),"trend")
      fits= lm(dy~lay+lxp+lxn+lagh1+trend,na.action =na.exclude )

      }
    }
    if(case==3){
       if(ncol(h)==2){

    ss=lapply(1:nrow(grid),function(x) lm(dy~nlagm1(y,grid[x,1])[-1,]+xn+nlagm1(xn,grid[x,3])+
                                            xp+nlagm1(xp,grid[x,2])
                                             +slag(h,grid[x,4:ncol(grid)])[-1,]))
      if(ic=="aic"){
        ak=sapply(1:nrow(grid),function(x) AICC(ss[[x]]))
      }
      if(ic=="bic"){
        ak=sapply(1:nrow(grid),function(x) BICC(ss[[x]]))
      }
      cii=as.matrix(ak)
      np=which.min(cii)
      pg=grid[np,]
      lay1=nlagm1(y,pg[,1])
      lay=as.matrix(nlagm1(y,pg[,1])[-1,])
      colnames(lay)=colnames(lay1)
      lxp=nlagm(xp,pg[,2])
      lxn=nlagm(xn,pg[,3])
      lagh1=slag(h,pg[,4:ncol(pg)])[-1,,drop=FALSE]
      rhnams<-c("Const",colnames(lay),colnames(lxp),colnames(lxn),colnames(lagh1))
      fits= lm(dy~lay+lxp+lxn+lagh1,na.action =na.exclude )


       }
      else{ ss=lapply(1:nrow(grid),function(x) lm(dy~nlagm1(y,grid[x,1])[-1,]+xn+nlagm1(xn,grid[x,3])+
                                                    xp+nlagm1(xp,grid[x,2])+nlagm(h,grid[x,4])[-1,]))


      if(ic=="aic"){
        ak=sapply(1:nrow(grid),function(x) AIC(ss[[x]]))
      }
      if(ic=="bic"){
        ak=sapply(1:nrow(grid),function(x) BIC(ss[[x]]))
      }
      cii=as.matrix(ak)
      np=which.min(cii)
      pg=grid[np,]
      lay1=nlagm1(y,pg[,1])
      lay=as.matrix(nlagm1(y,pg[,1])[-1,])
      colnames(lay)=colnames(lay1)
      lxp=nlagm(xp,pg[,2])
      lxn=nlagm(xn,pg[,3])
      lagh1=as.matrix(nlagm(h,pg[,4])[-1,,drop=FALSE])
      rhnams<-c("Const",colnames(lay),colnames(lxp),colnames(lxn),colnames(lagh1))
      fits= lm(dy~lay+lxp+lxn+lagh1,na.action =na.exclude )
      }
      #summaries=lapply(ss,summary)


    }
  }else{
    aa<-unlist(as.character(f))
    lhs   <- aa[[2]]
    core  <- aa[[3]]
    fm <- paste(lhs,"~",core)
    ffm<-as.formula(fm)
    fffm<-update(ffm, ~ . -1)
    mf <- model.frame(formula=fffm, data=data)
    x <- model.matrix(attr(mf, "terms"), data=mf)
    k<-ncol(x)
    if(ncol(x)>=2) stop("nardl package accept only one decomposed variable")
    y <- model.response(mf)
    y<-as.matrix(y)
    colnames(y)<-lhs
    dy<-diff(y)
    colnames(dy)<-paste("D",lhs,sep=".")
    dx<-diff(x)
    core=strsplit(a[[3]],"[+]")
    core=unlist(core)
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
    #*************************************************
    # grid lag selection
    #**********************************************

    l <- rep(list(0:maxlag),(k+2))
    grid0=expand.grid(l)
    b=which(grid0[1:nrow(grid0),]==0)
    grid=grid0[-b,]

    #************************************************************
    if(case==5){
      trend<-seq_along(dy)
      ss=lapply(1:nrow(grid),function(x) lm(dy~nlagm1(y,grid[x,1])[-1,]+xn+nlagm1(xn,grid[x,3])+
                                              xp+nlagm1(xp,grid[x,2])+trend))
     # summaries=lapply(ss,summary)
      if(ic=="aic"){
        ak=sapply(1:nrow(grid),function(x) AICC(ss[[x]]))
      }
      if(ic=="bic"){
        ak=sapply(1:nrow(grid),function(x) BICC(ss[[x]]))
      }
      cii=as.matrix(ak)
      np=which.min(cii)
      pg=grid[np,]
      lay1=nlagm1(y,pg[,1])
      lay=as.matrix(nlagm1(y,pg[,1])[-1,])
      colnames(lay)=colnames(lay1)
      lxp=nlagm(xp,pg[,2])
      lxn=nlagm(xn,pg[,3])
      rhnams<-c("Const",colnames(lay),colnames(lxp),colnames(lxn),"trend")
      fits= lm(dy~lay+lxp+lxn+trend,na.action =na.exclude )
    }
    if(case==3){
      ss=lapply(1:nrow(grid),function(x) lm(dy~nlagm1(y,grid[x,1])[-1,]+xn+nlagm1(xn,grid[x,2])+
                                              xp+nlagm1(xp,grid[x,3])))

      if(ic=="aic"){
        ak=sapply(1:nrow(grid),function(x) AICC(ss[[x]]))
      }
      if(ic=="bic"){
        ak=sapply(1:nrow(grid),function(x) BIC(ss[[x]]))
      }
      cii=as.matrix(ak)
      np=which.min(cii)
      pg=grid[np,]
      lay1=nlagm1(y,pg[,1])
      lay=as.matrix(nlagm1(y,pg[,1])[-1,])
      colnames(lay)=colnames(lay1)
      lxp=nlagm(xp,pg[,2])
      lxn=nlagm(xn,pg[,3])
      rhnams<-c("Const",colnames(lay),colnames(lxp),colnames(lxn))
      fits= lm(dy~lay+lxp+lxn,na.action =na.exclude )

    }
  }
  names(fits$coefficients)<-rhnams
  sels=summary(fits)
  #**************************
  #diagnos
  #**************************
  errors=sels$residuals
  nl=pg[[1]]
  jbtest<-shapiro.test(errors)
  lm2<-bp2(fits,nl,fill=0,type="F")
  arch<-ArchTest(errors,nl)
  #**************************
  #long run estimation
  #**************************
  coeff<-sels$coefficients
  nlvars<-length(coeff[,1])
  lvars<-coeff[3:nlvars,1]
  seldata<-data.matrix(coeff)
  coof<- -lvars/coeff[[2]]
  cof<- matrix(coof, length(lvars),1)
  #**************************
  #SE by delta method
  #**************************
  vc<-vcov(fits)#vcov(sel)
  vcc<-vc[2:nrow(vc),2:ncol(vc)]

  nsel<-length(sels$coefficients[,1])
  fb1<-lvars/coeff[[2]]^2
  fb2<-(-1/coeff[[2]])*diag(nrow(as.matrix(fb1)))
  fb<-cbind(as.matrix(fb1),fb2)
  lrse<-sqrt(diag(fb%*%vcc%*%t(fb)))
  #lrse2=lrse^2
  lrse2=lrse
  lrt<-coof/lrse2
  lrpv<-2 * pnorm(-abs(lrt))
  lres<-cbind(cof,lrse2,lrt,lrpv)
  #lres
  rownames(lres)<-names(lvars)
  colnames(lres)<-c("Estimate", "Std. Error", "t value", "Pr(>|t|)" )
  #**************************
  #cointegration test
  #**************************
  l<-c(paste(names(coeff[-1,1]),"=0"))
  coin<-linearHypothesis(fits,l,test="F")
  fstat<-coin$F[2]
  tcoin=fits$coefficients[2]
  #**************************
  #wald test
  #**************************
  rcap=matrix(c(0),1,2)
  rcap[1,1]=1
  rcap[1,2]=-1
  rsml=matrix(c(0),1,1)

  bp=paste(colnames(x),"p",sep="_")
  bn=paste(colnames(x),"n",sep="_")
  vcc1=vcc[c(bp,bn),c(bp,bn)]
  #***************************
  # long run wlad test
  #***************************
  lbetas=coof[c(bp,bn)]
  wldq=wtest(rcap,lbetas,vcc1,rsml)
  #**************************
  # short run wald test
  #**************************
  sbeta=coeff[c(bp,bn),1]
  wldsr=wtest(rcap,sbeta,vcc1,rsml)
  Nobs=length(fits$residuals)
  #************************
  # plot
  #************************
  #get recursive residuals
  rece<-strucchange::recresid(fits)
  # plot cumsum and cumsumq
  if(graph==TRUE){
    if(length(rece)<=1) stop("length of recursive residuals is less then 1! ")
    #x11()
    bbst<-sels$coefficients[,1]
    k<-length(bbst[-1])
    n<-length(errors)
    cusum(rece,k,n)
    cumsq(rece,k,n)
  }

  out<-list(fstat=fstat,tcoin=tcoin, fits=fits, sels=sels, cof=cof,coof=coof,
            k=k,n=n,case=case,lvars=lvars,selresidu=errors,
            jbtest=jbtest,arch=arch,
            nl=nl,rece=rece,Nobs=Nobs,vcc=vcc,fb=fb,
            fb1=fb1,fb2=fb2,lrse=lrse,lrt=lrt,lrpv=lrpv,
            lres=lres,lm2=lm2,wldsr=wldsr,wldq=wldq,pg=pg,
            ak=ak)

  class(out) <- "nardl"

  return(out)

}

