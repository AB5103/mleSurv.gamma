#' @title LikMarg
#'
#' @description This function calculates neg-loglikelihood for a proportional hazards model.
#'
#' @param y Vector of parameters in the form:
#' \deqn{y = (ln(a), \beta , ln(\sigma ^2))}
#' for the exponential baseline hazard function,
#' \deqn{y = (ln(a), \beta , ln(\sigma ^2),b)}
#' for the Weibull baseline hazard function and
#' \deqn{y = (ln(10^3a), \beta , ln(\sigma ^2), 10b)}
#' for the Gompertz baseline hazard function, where a and b are the slope and the shape parameters,
#' \eqn{\beta } are the Cox-regression parameters, and \eqn{\sigma ^2} is the variance of frailty.
#' @param D  A data.frame in which to interpret the variables named in the formula.
#' The data set includes the following fields:
#' \enumerate{
#' \item time-to-failure 'time';
#' \item censoring 'event' (censoring must be either 0 for no event or 1 for event);
#' \item  covariates (continuous or categorical, first 'nf' columns) used in a study (can be empty set);
#' \item cluster - the name of a cluster variable in data (is equal to NULL for the fixed-effect model).}
#' @param nf The number of covariates under study.
#' @param ncl The number of clusters.
#' @param dist Baseline hazard function ('Exponential', 'Weibull' or 'Gompertz').
#'
#' @details Three kinds of the baseline hazards are used in this function:
#' \enumerate{
#' \item Exponential with baseline cumulative hazard function \deqn{H_{0}(t;a)=at;}
#' \item Weibull with baseline cumulative hazard function \deqn{H_{0}(t;a,b)=at^b;}
#' \item Gompertz with baseline cumulative hazard function \deqn{H_{0}(t;a,b)=a(exp(bt)-1)).}
#' }
#'
#' @return Neg-loglikelihood
#'
#' @examples
#' \dontrun{
#' LikMarg(y,D,nf,ncl,dist)
#' }
#' @export
#'
LikMarg=function(y,D,nf,ncl,dist){
     k0=1
    if (dist=='Weibull'){
    lambda0=exp(y[1])
    k0=y[nf+3]} else if (dist=='Gompertz') {
      lambda0=1e-3*exp(y[1])
      k0=1e-1*y[nf+3]
    } else if (dist=='Exponential'){
      lambda0=exp(y[1])
    }
  if (k0>0){
  Cox=c(rep(0,nrow(D)))
  if (nf>0){
    for (i in 1:nf){
        Cox=Cox+D[,i]*y[1+i]
    }
  }
  LCox=Cox
  Cox=exp(Cox)
  if (ncl>0) {
    G2=exp(y[2+nf])
    list=unique(D$cluster)
    nl=length(list)
      }
  Cens=D$event
  x1=D$time

  if (dist=='Exponential'){
    if (G2>1e-8 & ncl>0){
      Lik=0
      for (i in 1:nl){
        ID=list[i]
        ind=which(D$cluster==ID)
        nn=length(ind)
        ic1=1*(Cens[ind]==1)
        nc1=sum(ic1)
        Hfull1= Cox[ind]*lambda0*x1[ind]
        Lmufull1=LCox[ind]+log(lambda0)
        if (nc1<=1) {
          ee1=0} else {
            ee1=sum(log(c(1:(nc1-1))*G2+1))
          }
        Lik=Lik+ee1+sum(Lmufull1*ic1)-(1/G2+nc1)*log(1+G2*sum(Hfull1))
      }
    } else {
      ic=1*(Cens==1)
      Hfull1= Cox*lambda0*x1
      Lmufull1=LCox+log(lambda0)
      Lik=sum(Lmufull1*ic)-sum(Hfull1)
    }
  }


  if (dist=='Weibull'){
    if (G2>1e-8 & ncl>0){
      Lik=0
      for (i in 1:nl){
        ID=list[i]
        ind=which(D$cluster==ID)
        nn=length(ind)
        ic1=1*(Cens[ind]==1)
        nc1=sum(ic1)
        Hfull1= Cox[ind]*lambda0*x1[ind]^k0
        Lmufull1=LCox[ind]+log(lambda0)+log(k0)+(k0-1)*log(x1[ind])
        if (nc1<=1) {
          ee1=0} else {
            ee1=sum(log(c(1:(nc1-1))*G2+1))
          }
        Lik=Lik+ee1+sum(Lmufull1*ic1)-(1/G2+nc1)*log(1+G2*sum(Hfull1))
      }
    } else {
      ic=1*(Cens==1)
      Hfull1= Cox*lambda0*x1^k0
      Lmufull1=LCox+log(lambda0)+log(k0)+(k0-1)*log(x1)
      Lik=sum(Lmufull1*ic)-sum(Hfull1)
    }
  }
  if (dist=='Gompertz'){
    if (G2>1e-8 & ncl>0){
      Lik=0
      for (i in 1:nl){
        ID=list[i]
        ind=which(D$cluster==ID)
        nn=length(ind)
        ic1=1*(Cens[ind]==1)
        nc1=sum(ic1)
        Hfull1= Cox[ind]*lambda0*(exp(k0*x1[ind])-1)
        Lmufull1=LCox[ind]+log(lambda0)+log(k0)+k0*x1[ind]
        if (nc1<=1) {
          ee1=0} else {
            ee1=sum(log(c(1:(nc1-1))*G2+1))
          }
        Lik=Lik+ee1+sum(Lmufull1*ic1)-(1/G2+nc1)*log(1+G2*sum(Hfull1))
      }
    } else {
      ic=1*(Cens==1)
      Hfull1= Cox*lambda0*(exp(k0*x1)-1)
      Lmufull1=LCox+log(lambda0)+log(k0)+k0*x1
      Lik=sum(Lmufull1*ic)-sum(Hfull1)

    }
  }
  Lik=-Lik
  if (is.na(Lik) | is.infinite(Lik)) Lik=1e+50
  } else {Lik=1e+50}
  return(Lik)
}
