% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bcmle.R
\name{bcmle}
\alias{bcmle}
\title{bcmle}
\usage{
bcmle(formula, data, cluster, dist)
}
\arguments{
\item{formula}{A formula object, with the response on the left of a ~ operator, and the terms on
the right. The response must be a survival object as returned by the Surv() function.
The status indicator corresponding to censoring in the Surv object must be 0 for censored objects and 1 for non-censored ones.
This object describes the effect of several factors in the Cox proportional hazards model.}

\item{data}{The data frame. Includes fields for:
\enumerate{
\item Time-to-failure;
\item Censoring;
\item Possible covariates (optional);
\item Cluster.
}
The names of these fields must be included in 'formula'.}

\item{cluster}{The name of a cluster variable in data.}

\item{dist}{Baseline hazard function ('Exponential', 'Weibull' or 'Gompertz').}
}
\value{
The list object which components are:
\enumerate{
\item baseline hazard function ('Exponential', 'Weibull' or 'Gompertz') - 'dist';
\item ML parameter estimates 'parEst';
\item standard errors for ML parameter estimates - 'stErr';
\item bias corrected parameter estimates  - 'parCor';
\item latex table for ML estimates, standard errors and bias corrected ML estimates - 'xTab';
\item times to failure in increasing order - 't';
\item empirical survivals - 'sEmp';
\item survivals based on the ML estimates - 'sEst';
\item survivals based on the bias corrected estimates - 'sCor'.
}
}
\description{
Calculates ML estimates and their bias corrected values for the proportional hazards model with gamma frailty.
}
\details{
Bias correction is computed using the Cox-Snell method [1]-[2] involving calculation the expectations for
the second and third order derivatives of the log-likelihood function with respect to parameters of the model.
If the calculation of the bias fails return NAs.
\deqn{}
Remark. Ever if the calculation of the bias does not fail it does not guarantee that the calculated bias is reliable.
To test it it is recommended to produce the plot for empirical, estimated and corrected survivals using function 'plotSurv' and check
fit visually.
}
\examples{
\dontrun{
library(survival)
set.seed(1)
formula=as.formula("Surv(time, status) ~ Sex+Score")
dist="Exponential"
NN=100
ncl=5
npercl=floor(NN/ncl)
status=rep(1,NN)
Sex=rbinom(NN,1,0.5)
Score=rnorm(NN,0,0.5)
ID=rep(1:ncl,npercl)
event=rep(1,NN)
G2=1
scale=0.2
bSex=0.5
bScore=1
Cox=exp(bSex*Sex+bScore*Score)
Z=rep(rgamma(ncl,shape=1/G2,rate=1/G2),npercl)
time=-log(runif(NN))/(Z*scale*Cox)
D=data.frame(Sex,Score,time,event,ID)
colnames(D)=c("Sex","Score","time","status","id")
cluster="id"
R.model=bcmle(formula,D,cluster,dist)
plotSurv(R.model)
}
\dontrun{
library(survival)
set.seed(1)
formula=as.formula("Surv(time, status) ~ Sex+Score")
dist="Weibull"
NN=100
ncl=5
npercl=floor(NN/ncl)
status=rep(1,NN)
Sex=rbinom(NN,1,0.5)
Score=rnorm(NN,0,0.5)
ID=rep(1:ncl,npercl)
event=rep(1,NN)
G2=1
scale=0.05
shape=1.5
bSex=0.5
bScore=1
Cox=exp(bSex*Sex+bScore*Score)
Z=rep(rgamma(ncl,shape=1/G2,rate=1/G2),npercl)
time=(-log(runif(NN))/(Z*scale*Cox))^(1/shape)
D=data.frame(Sex,Score,time,event,ID)
colnames(D)=c("Sex","Score","time","status","id")
cluster="id"
R.model=bcmle(formula,D,cluster,dist)
plotSurv(R.model)
}

\dontrun{
library(survival)
set.seed(1)
formula=as.formula("Surv(time, status) ~ Sex+Score")
dist="Gompertz"
NN=100
ncl=5
npercl=floor(NN/ncl)
status=rep(1,NN)
Sex=rbinom(NN,1,0.5)
Score=rnorm(NN,0,0.5)
ID=rep(1:ncl,npercl)
event=rep(1,NN)
G2=1
scale=1e-3
shape=0.1
bSex=0.5
bScore=1
Cox=exp(bSex*Sex+bScore*Score)
Z=rep(rgamma(ncl,shape=1/G2,rate=1/G2),npercl)
time=log((-log(runif(NN))/(Z*scale*Cox))+1)/shape
D=data.frame(Sex,Score,time,event,ID)
colnames(D)=c("Sex","Score","time","status","id")
cluster="id"
R.model=bcmle(formula,D,cluster,dist)
plotSurv(R.model)
}

}
\references{
1. Cox D., Snell E. A general definition of residuals. \eqn{R Stat Soc series B} 1968; 30(2):248--275.

2. Cordeiro G, Klein R. Bias correction in ARMA models. \eqn{Stat Probab Lett} 1994; 19(3):169--176.
}
