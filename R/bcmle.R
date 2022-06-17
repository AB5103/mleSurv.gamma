#' @title bcmle
#'
#' @description   Calculates ML estimates and their bias corrected values for the proportional hazards model with gamma frailty.
#'
#' @param formula A formula object, with the response on the left of a ~ operator, and the terms on
#' the right. The response must be a survival object as returned by the Surv() function.
#' The status indicator corresponding to censoring in the Surv object must be 0 for censored objects and 1 for non-censored ones.
#' This object describes the effect of several factors in the Cox proportional hazards model.
#' @param data The data frame. Includes fields for:
#' \enumerate{
#' \item Time-to-failure;
#' \item Censoring;
#' \item Possible covariates (optional);
#' \item Cluster.
#' }
#' The names of these fields must be included in 'formula'.
#' @param cluster The name of a cluster variable in data.
#' @param dist Baseline hazard function ('Exponential', 'Weibull' or 'Gompertz').
#'
#' @return The list object which components are:
#' \enumerate{
#' \item baseline hazard function ('Exponential', 'Weibull' or 'Gompertz') - 'dist';
#' \item ML parameter estimates 'parEst';
#' \item standard errors for ML parameter estimates - 'stErr';
#' \item bias corrected parameter estimates  - 'parCor';
#' \item latex table for ML estimates, standard errors and bias corrected ML estimates - 'xTab';
#' \item times to failure in increasing order - 't';
#' \item empirical survivals - 'sEmp';
#' \item survivals based on the ML estimates - 'sEst';
#' \item survivals based on the bias corrected estimates - 'sCor'.
#' }
#' @details Bias correction is computed using the Cox-Snell method [1]-[2] involving calculation the expectations for
#' the second and third order derivatives of the log-likelihood function with respect to parameters of the model.
#' If the calculation of the bias fails return NAs.
#' \deqn{}
#' Remark. Ever if the calculation of the bias does not fail it does not guarantee that the calculated bias is reliable.
#' To test it it is recommended to produce the plot for empirical, estimated and corrected survivals using function 'plotSurv' and check
#' fit visually.
#'
#' @references
#'
#' 1. Cox D., Snell E. A general definition of residuals. \eqn{R Stat Soc series B} 1968; 30(2):248--275.
#'
#' 2. Cordeiro G, Klein R. Bias correction in ARMA models. \eqn{Stat Probab Lett} 1994; 19(3):169--176.
#'
#' @import survival
#' @import MASS
#' @import ucminf
#' @import stats
#' @import xtable
#'
#' @examples
#' \dontrun{
#' library(survival)
#' set.seed(1)
#' formula=as.formula("Surv(time, status) ~ Sex+Score")
#' dist="Exponential"
#' NN=100
#' ncl=5
#' npercl=floor(NN/ncl)
#' status=rep(1,NN)
#' Sex=rbinom(NN,1,0.5)
#' Score=rnorm(NN,0,0.5)
#' ID=rep(1:ncl,npercl)
#' event=rep(1,NN)
#' G2=1
#' scale=0.2
#' bSex=0.5
#' bScore=1
#' Cox=exp(bSex*Sex+bScore*Score)
#' Z=rep(rgamma(ncl,shape=1/G2,rate=1/G2),npercl)
#' time=-log(runif(NN))/(Z*scale*Cox)
#' D=data.frame(Sex,Score,time,event,ID)
#' colnames(D)=c("Sex","Score","time","status","id")
#' cluster="id"
#' R.model=bcmle(formula,D,cluster,dist)
#' plotSurv(R.model)
#' }
#' \dontrun{
#' library(survival)
#' set.seed(1)
#' formula=as.formula("Surv(time, status) ~ Sex+Score")
#' dist="Weibull"
#' NN=100
#' ncl=5
#' npercl=floor(NN/ncl)
#' status=rep(1,NN)
#' Sex=rbinom(NN,1,0.5)
#' Score=rnorm(NN,0,0.5)
#' ID=rep(1:ncl,npercl)
#' event=rep(1,NN)
#' G2=1
#' scale=0.05
#' shape=1.5
#' bSex=0.5
#' bScore=1
#' Cox=exp(bSex*Sex+bScore*Score)
#' Z=rep(rgamma(ncl,shape=1/G2,rate=1/G2),npercl)
#' time=(-log(runif(NN))/(Z*scale*Cox))^(1/shape)
#' D=data.frame(Sex,Score,time,event,ID)
#' colnames(D)=c("Sex","Score","time","status","id")
#' cluster="id"
#' R.model=bcmle(formula,D,cluster,dist)
#' plotSurv(R.model)
#' }
#'
#' \dontrun{
#' library(survival)
#' set.seed(1)
#' formula=as.formula("Surv(time, status) ~ Sex+Score")
#' dist="Gompertz"
#' NN=100
#' ncl=5
#' npercl=floor(NN/ncl)
#' status=rep(1,NN)
#' Sex=rbinom(NN,1,0.5)
#' Score=rnorm(NN,0,0.5)
#' ID=rep(1:ncl,npercl)
#' event=rep(1,NN)
#' G2=1
#' scale=1e-3
#' shape=0.1
#' bSex=0.5
#' bScore=1
#' Cox=exp(bSex*Sex+bScore*Score)
#' Z=rep(rgamma(ncl,shape=1/G2,rate=1/G2),npercl)
#' time=log((-log(runif(NN))/(Z*scale*Cox))+1)/shape
#' D=data.frame(Sex,Score,time,event,ID)
#' colnames(D)=c("Sex","Score","time","status","id")
#' cluster="id"
#' R.model=bcmle(formula,D,cluster,dist)
#' plotSurv(R.model)
#' }
#'
#' @export
#'
bcmle=function(formula,data,cluster,dist)
{
  nr=nrow(data)
  obsdata <- NULL
  obsdata$time <- eval(formula[[2]][[2]], envir = data)
  obsdata$event <- eval(formula[[2]][[3]], envir = data)
  obsdata$x <- as.data.frame(model.matrix.lm(formula, data = data,na.action='na.pass'))
  names(obsdata$x)=names(obsdata$x)

  ind.x=which(is.na(c(apply(obsdata$x,1,sum))))

  if (is.null(cluster)) {
    obsdata$ncl <- 0
    obsdata$di <- sum(obsdata$event)
    obsdata$cluster <- c(rep(0,nrow(data)))
    ind.cl <- as.numeric({})
    stop(paste0("object '", cluster, "' not found"))
  }   else {
    if (!cluster %in% names(data)) {
      stop(paste0("object '", cluster, "' not found"))
    }
    obsdata$cluster <- eval(as.name(cluster), envir = data)
    obsdata$ncl <- length(levels(as.factor(obsdata$cluster)))
    obsdata$di <- aggregate(obsdata$event, by = list(obsdata$cluster),FUN = sum)[, , drop = FALSE]
    cnames <- obsdata$di[, 1]
    obsdata$di <- as.vector(obsdata$di[, 2])
    names(obsdata$di) <- cnames
    ind.cl=which(is.na(obsdata$cl))
  }
  ind=sort(unique(c(ind.x,ind.cl)))
  nx=rep(0,nr)
  nx[ind]=1
  if (is.factor(obsdata$cluster)) obsdata$cluster=as.character(obsdata$cluster)
  obs=data.frame(obsdata$x[-1],obsdata$event,obsdata$time,obsdata$cluster)[nx!=1,]
  nr=nrow(obs)
  namesf=names(obsdata$x)[-1]
  names(obs)=c(names(obsdata$x)[-1],names(obsdata$event),"event","time","cluster")
  ncl=obsdata$ncl
  nf=length(namesf)
  Fact=matrix(unlist(obsdata$x[-1]),nrow(obsdata$x[-1]),nf)
  par0=c(0,0)
  if (nf>0)  par0=c(0,c(rep(0,nf)),0)
  if (!(dist=='Exponential')){
  par0=c(0,0,1)
  if (nf>0)  par0=c(0,c(rep(0,nf)),0,1)
  }
  ResultMarg=ucminf(par=par0,fn=LikMarg,gr=NULL,D=obs,nf=nf,ncl=ncl,dist=dist,hessian=3)
  para=ResultMarg$par
    if (any(!is.finite(as.matrix(ResultMarg$hessian))))
    stop("infinite or missing values in hessian. It is not possible to calculate the matrix of covariates. \n  Change the model and try again.")
    if (any(suppressWarnings(diag(ginv(ResultMarg$hessian)))<0))
    stop("hessian cannot be correctly calculated. \n  Change the model and try again.")
    se=sqrt(diag(ginv(ResultMarg$hessian)))
    Lik=-ResultMarg$value

  if (dist=='Exponential'){
    Bia=BiasCorrectionE(obs,nf,ncl,para)
    piar=para
    se=se} else if (dist=='Weibull'){
      Bia=BiasCorrectionW(obs,nf,ncl,para)
      piar=para
      se=se} else if (dist=='Gompertz') {
        Bia=BiasCorrectionG(obs,nf,ncl,para)
        piar=c(para[1]+log(1e-3),para[2:(nf+2)],0.1*para[nf+3])
        se=c(se[1]+log(1e-3),se[2:(nf+2)],0.1*se[nf+3])
      }
   parcor=piar-Bia

   if (dist=='Exponential' & nf==0) namVar=c("$\\ln(a)$","$\\ln(\\sigma ^2)$")
   if (dist!='Exponential' & nf==0) namVar=c("$\\ln(a)$","$log(\\sigma ^2)$","$b$")
   if (dist=='Exponential' & nf>0){
     namVar="$\\ln(a)$"
     for (it in 1:nf){
       namVar=c(namVar,paste0("$\\beta _{",namesf[it],"}$"))
     }
     namVar=c(namVar,"$\\ln(\\sigma ^2)$")
   }
   if (dist!='Exponential' & nf>0){
     namVar="$\\ln(a)$"
     for (it in 1:nf){
       namVar=c(namVar,paste0("$\\beta _{",namesf[it],"}$"))
     }
     namVar=c(namVar,"$\\ln(\\sigma ^2)$","$b$")
   }
   xTab=data.frame(piar,se,parcor)
   colnames(xTab)=c("MLE","SE","Corrected MLE")
   rownames(xTab)=namVar
   xTab=xtable(xTab)
   digits(xTab)=3
   align(xTab)=c("l",rep("c",3))
   xTab=print(xTab,sanitize.rownames.function=function(x){x})
  MeanCov=0
  if (nf>0){
  MeanCov=c(apply(Fact,2,mean))
  }
  Coxp=1
  if (nf>0) Coxp=exp(sum(MeanCov*para[2:(1+nf)]))
  km_fit <- survfit(Surv(time, event) ~ 1, data=obs)
  t=summary(km_fit)$time
  sEmp=summary(km_fit)$surv
  if (dist=='Exponential') {
    ae=exp(piar[1])*Coxp
    G2=exp(piar[nf+2])
    acor=exp(parcor[1])*Coxp
    G2cor=exp(parcor[nf+2])} else {
    ae=exp(piar[1])*Coxp
    G2=exp(piar[nf+2])
    be=piar[nf+3]
    acor=exp(parcor[1])*Coxp
    G2cor=exp(parcor[nf+2])
    bcor=parcor[nf+3]
    }
  sEst=(1+G2*ae*t)^(-1/G2)
  sCor=(1+G2cor*acor*t)^(-1/G2cor)
  if (dist=='Weibull') {
    sEst=(1+G2*ae*t^be)^(-1/G2)
    sCor=(1+G2cor*acor*t^bcor)^(-1/G2cor)
  } else if (dist=='Gompertz'){
    sEst=(1+G2*ae*(exp(be*t)-1))^(-1/G2)
    sCor=(1+G2cor*acor*(exp(bcor*t)-1))^(-1/G2cor)
  }
  return(list(dist=dist,parEst=piar,stErr=se,parCor=parcor,xTab=xTab,time=t,sEmp=sEmp,sEst=sEst,sCor=sCor))
}
