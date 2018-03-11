
mplier<-function(alpha,beta,varphi,vpi,l,h,k){
  # k number of independent variables
  # h is the horizon over which multipliers will be computed
  # l max(p,q)
  vphi<-matlab::zeros(1,l)
  vphi[1] = 1 + alpha + varphi[1]
  i=2
  while(i<=l){
    vphi[i] = varphi[i] - varphi[i-1]
   i=i+1
  }

  vphi[l] = -varphi[l-1]

  mtheta = matlab::zeros(l+1,k)
  mtheta[1,] = vpi[1]
  mtheta[2:3,] = vpi[2] - vpi[1] + beta
  i=3
  while(i<=l+1){
    mtheta[i,] = vpi[i] - vpi[i-1]
    i=i+1
  }
  mtheta[l+1,] = -vpi[l]
  mpsi =matlab::zeros(h,k)
  #mpsi<-NULL
  mpsi[1,] = mtheta[1,]
  i=1
  while(i<=l){
    mpsi[i+1,] = vphi[1:i]%*%mpsi[i:1,] + mtheta[i+1,]

    i=i+1
  }

  i=l+1
  mps<-vector('list',length(mpsi))
  while(i<=h-1){
    mps[[i+1]]<-t(vphi)%*%mpsi[i:i-l+1]
    #mpsi<-c(mpsi[i+1],t(vphi)%*%mpsi[i-l+1:i])
    i=i+1
  }
  mps<-matrix(do.call(rbind,mps)[1:h,])
  cm<-apply(mps, 2, cumsum)
  return(cm)
}

#'Dynamic multiplier plot
#'
#'@param model the fitted model
#'@param np the selected number of lags
#'@param k number of decomposed independent variables
#'@param h is the horizon over which multipliers will be computed
#'@importFrom stats lm AIC BIC as.formula model.frame model.matrix model.response na.omit sd update vcov embed resid coef
#'@importFrom graphics abline legend lines par plot
#'@import matlab
#'@examples
#'
#' ############################
#' # Dynamic multipliers plot
#' ############################
#' # Load data
#'data(fod)
#' reg<-nardl(food~inf,p=4,q=4,fod,ic="aic",maxlags = FALSE,graph = TRUE,case=3)
#' plotmplier(reg,reg$np,1,10)
#'
#'@export
plotmplier<-function(model,np,k,h){

  if(np<=1) stop("number of lags must be > 1 ")
  #positive
  alpha<-as.matrix(coef(model$fit)[2:2])
  betap<-as.matrix(coef(model$fit))[3:3]
  betan<-as.matrix(coef(model$fit))[4:4]
  varphi<-as.matrix(coef(model$fit)[5:(4+np)])
  vpip<-as.matrix(coef(model$fit)[(5+np):((4+np)+np)])
  vpin<-as.matrix(coef(model$fit)[((5+np)+np):(((4+np)+np)+np)])
  mpp<-mplier(alpha,betap,varphi,vpip,np,h,k)
  mpn<-mplier(alpha,betan,varphi,vpin,np,h,k)
  x = seqa(0,1,(h+1))
  yp<-cbind(matlab::zeros(1,1),t(mpp))
  yn = cbind(matlab::zeros(1,1),t(mpn))
  dmp<-mpp-mpn
  yrange<-range(yp,yn,(yp-yn))
  par(mar = c(5,4,4,8))
  plot(x,yp,main="Dynamic multiplier",type = "l",ylim =yrange ,xlab ="",
       ylab = "", col="red")
  lines(x,yn,col="blue")
  lines(x,(yp-yn),lty=2,col="darkgreen")
  abline(h=0)
  legend(par("usr")[2],par("usr")[4],
         xpd = TRUE ,
         bty = "n",
   c(names(coef(model$fit)[3:3]),names(coef(model$fit)[4:4]),"Difference"),
  lty = c(1, 1, 2),
  cex=0.6,
  col=c("blue","red","darkgreen"))



}
