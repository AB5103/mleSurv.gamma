#' @title BiasCorrectionG
#'
#' @description   Returns the bias correction for ML estimates.
#'
#' @param D The data frame. The first 'nf' columns are the covariates under study. Other columns must include fields
#' 'time' for time-to-failure, 'Cens' for censoring and 'cluster' for cluster identification. The number of rows is equal to the number of subjects.
#' @param nf The number of covariates under study.
#' @param ncl The number of clusters.
#' @param para  The vector of the ML estimates of unknown parameters. This vector has a form \deqn{y = (ln(10^3a), \beta , ln(\sigma ^2), 10b),}
#' where \eqn{a} and \eqn{b} are the slope and the shape parameters, \eqn{\beta } is the vector of the Cox-regression parameters,
#' and \eqn{\sigma ^2} is the variance of the frailty.
#' @details This function calculates the bias correction for the proportional hazards model with gamma frailty
#' in accordance with the Cox-Snell method.  The baseline cumulative hazard function is \deqn{H_{0}(t;a,b)=a(e^{bt}-1)}.
#'
#' @return The bias correction for vector \eqn{(ln(a), \beta , ln(\sigma ^2), b)}. If the calculation of the bias correction is
#' not possible returns NAs.
#'
#' @import MASS
#' @import Deriv
#'
#' @examples
#' \dontrun{
#' BiasCorrectionG(D,nf,ncl,para)
#'
#'}
#' @export
#'
BiasCorrectionG=function(D,nf,ncl,para){
  parnam=c("a",paste0("beta_",as.character(1:nf)),"g","b")
  nr=nrow(D)
  assign("delta",1)
  np=nf+3
  a.scale=rep(1,nr)
  g.scale=rep(1,nr)
  uR_beta=D[,1:nf]
  UR=cbind(a.scale,uR_beta,g.scale)
  Cens=D$event
  assign("a",para[1]+log(1e-3))
  assign("b",0.1*para[nf+3])
  assign("g",para[nf+2])
  fi="a"
  if (nf>0){
    par=c(a,para[2:(1+nf)],g,b)
    parnam=c('a',paste0("beta_",as.character(1:nf)),'g','b')
    for (i in 1:nf){
      assign(paste0('beta_',as.character(i)),par[1+i])
      fi=c(fi,paste0('beta_',as.character(i)))
    }
  } else {
    par=c(a,g,b)
    parnam=c('a','g','b')
  }
  assign("n",100)
  assign("d",100)
  par=c(par,c(n,d))
  parnam=c(parnam,c("n","d"))
  fi=c(fi,"g")
  list=sort(unique(D$cluster))
  RCC=array(0,c(np,np))
  RCCC=array(0,c(np,np,np))
  RCC_C=array(0,c(np,np,np))

  for (lj in 1:ncl){
    Parmen={}
    Parnam={}
    Parmen=c(Parmen,par)
    Parnam=c(Parnam,parnam)
    ind=which(D$cluster==list[lj])
    n<-length(ind)
    d<-sum(Cens[ind])
    Parmen[which(Parnam=="n")]=n
    Parmen[which(Parnam=="d")]=d
    li=length(ind)
    if (li==1)  {U=matrix(as.numeric(UR[ind,]),length(ind),ncol(UR))} else {
      U=UR[ind,]}

      WW1={}
      WW2={}
      WW3={}
      WW4={}
      WW5={}
      WW6={}
      for (as1 in 1:(nf+2)){
      work1=matrix(rep(matrix(rep(U[,as1],n),n,n),n),n^2,n)
      work2=t(matrix(rep(matrix(rep(U[,as1],n),n,n),n),n,n^2))
      work3=matrix(rep(t(matrix(rep(U[,as1],n),n,n)),n),n^2,n)
      work4=matrix(rep(U[,as1],n),n,n)
      work5=t(matrix(rep(U[,as1],n),n,n))
      work6=matrix(U[,as1],n,1)
      WW1=cbind(WW1,work1)
      WW2=cbind(WW2,work2)
      WW3=cbind(WW3,work3)
      WW4=cbind(WW4,work4)
      WW5=cbind(WW5,work5)
      WW6=cbind(WW6,work6)
      assign(paste0("U3_",as.character(as1)), work1)
      assign(paste0("tU3_",as.character(as1)),work2)
      assign(paste0("U3t_",as.character(as1)),work3)
      assign(paste0("U2_",as.character(as1)), work4)
      assign(paste0("tU2_",as.character(as1)),work5)
      assign(paste0("U1_",as.character(as1)), work6)
    }
    assign("k",100)
    assign("k1",100)
    assign("k2",100)
    assign("k3",100)
    assign("pl",100)
    assign("ND1",100)
    assign("ND2",100)
    assign("ND3",100)
    assign("Deg0",100)
    assign("c1",1)
    assign("c2",1)
    assign("c3",1)

    Parnam=c(Parnam,"k","k1","k2","k3","pl","ND1","ND2","ND3","Deg0","c1","c2","c3")
    Parmen=c(Parmen,k,k1,k2,k3,pl,ND1,ND2,ND3,Deg0,c1,c2,c3)
    nul=0
    u=c(rep(NA,nf+2))
    uu=matrix(NA,nf+2,nf+2)
    uuu=array(NA,c(nf+2,nf+2,nf+2))
    for (j1 in 1:(nf+2)){
      u[j1]=sum(U[,j1])
      for (j2 in 1:(nf+2)){
        uu[j1,j2]=sum(U[,j1]*U[,j2])
        for (j3 in 1:(nf+2)){
          uuu[j1,j2,j3]=sum(U[,j1]*U[,j2]*U[,j3])
        }
      }
    }

    lna="("
    for (i in 1:(nf+2)){
      if (i<(nf+2))  lna=paste0(lna,paste0("x_",as.character(i),""),"*",parnam[i],"+")
      if (i==(nf+2)) lna=paste0(lna,paste0("x_",as.character(i),""),"*",parnam[i],")")
    }
    ea1=paste0("(c1*exp",lna,")")
    lnaa=gsub("x_","x__",lna)
    eaa1=paste0("(c2*exp",lnaa,")")
    lnaaa=gsub("x_","x___",lna)
    eaaa1=paste0("(c3*exp",lnaaa,")")


    nu="exp(-g)"
    nu1="(exp(-g)+1)"
    nu2="(exp(-g)+2)"
    nu3="(exp(-g)+3)"
    nu4="(exp(-g)+4)"
    nu5="(exp(-g)+5)"
    nu6="(exp(-g)+6)"
    nud="(exp(-g)+d)"
    nun="(exp(-g)+n)"
    nun1="(exp(-g)+n+1)"
    nun2="(exp(-g)+n+2)"

    assign("mu",eval(parse(text=nu)))
    assign("mu1",eval(parse(text=nu1)))
    assign("mu2",eval(parse(text=nu2)))
    assign("mu3",eval(parse(text=nu3)))
    assign("mu4",eval(parse(text=nu4)))
    assign("mu5",eval(parse(text=nu5)))
    assign("mu6",eval(parse(text=nu6)))
    assign("mud",eval(parse(text=nud)))
    assign("mun",eval(parse(text=nun)))
    assign("mun1",eval(parse(text=nun1)))
    assign("mun2",eval(parse(text=nun2)))
    Parnam=c(Parnam,c("mu","mu1","mu2","mu3","mu4","mu5","mu6","mud","mun","mun1","mun2"))
    Parmen=c(Parmen,c(mu,mu1,mu2,mu3,mu4,mu5,mu6,mud,mun,mun1,mun2))
    for (im1 in 0:5){
      for (im2 in 0:3){
        for (im3 in 0:(nf+2)){
          if(im3!=0){
            assign(paste0("la1_",as.character(im3),"_",as.character(im1),"_",as.character(im2)),eval(parse(text=subsi1c(paste0("x_",as.character(im3),"*",lna,"^",as.character(im2),"*",ea1,"^",as.character(im1)),nf,Parnam,Parmen,WW6))))
            Parnam=c(Parnam,paste0("la1_",as.character(im3),"_",as.character(im1),"_",as.character(im2)))
            Parmen=c(Parmen,eval(parse(text=subsi1c(paste0("x_",as.character(im3),"*",lna,"^",as.character(im2),"*",ea1,"^",as.character(im1)),nf,Parnam,Parmen,WW6))))
            assign(paste0("la2_",as.character(im3),"_",as.character(im1),"_",as.character(im2)),eval(parse(text=subsi1c(paste0("x_",as.character(im3),"*",lna,"^",as.character(im2),"*",ea1,"^",as.character(im1)),nf,Parnam,Parmen,WW6))))
            Parnam=c(Parnam,paste0("la2_",as.character(im3),"_",as.character(im1),"_",as.character(im2)))
            Parmen=c(Parmen,eval(parse(text=subsi1c(paste0("x_",as.character(im3),"*",lna,"^",as.character(im2),"*",ea1,"^",as.character(im1)),nf,Parnam,Parmen,WW6))))
            assign(paste0("la3_",as.character(im3),"_",as.character(im1),"_",as.character(im2)),eval(parse(text=subsi1c(paste0("x_",as.character(im3),"*",lna,"^",as.character(im2),"*",ea1,"^",as.character(im1)),nf,Parnam,Parmen,WW6))))
            Parnam=c(Parnam,paste0("la3_",as.character(im3),"_",as.character(im1),"_",as.character(im2)))
            Parmen=c(Parmen,eval(parse(text=subsi1c(paste0("x_",as.character(im3),"*",lna,"^",as.character(im2),"*",ea1,"^",as.character(im1)),nf,Parnam,Parmen,WW6))))
          } else {
            assign(paste0("la1_0_",as.character(im1),"_",as.character(im2)),eval(parse(text=subsi1c(paste0(lna,"^",as.character(im2),"*",ea1,"^",as.character(im1)),nf,Parnam,Parmen,WW6))))
            Parnam=c(Parnam,paste0("la1_0_",as.character(im1),"_",as.character(im2)))
            Parmen=c(Parmen,eval(parse(text=subsi1c(paste0(lna,"^",as.character(im2),"*",ea1,"^",as.character(im1)),nf,Parnam,Parmen,WW6))))
            assign(paste0("la2_0_",as.character(im1),"_",as.character(im2)),eval(parse(text=subsi1c(paste0(lna,"^",as.character(im2),"*",ea1,"^",as.character(im1)),nf,Parnam,Parmen,WW6))))
            Parnam=c(Parnam,paste0("la2_0_",as.character(im1),"_",as.character(im2)))
            Parmen=c(Parmen,eval(parse(text=subsi1c(paste0(lna,"^",as.character(im2),"*",ea1,"^",as.character(im1)),nf,Parnam,Parmen,WW6))))
            assign(paste0("la3_0_",as.character(im1),"_",as.character(im2)),eval(parse(text=subsi1c(paste0(lna,"^",as.character(im2),"*",ea1,"^",as.character(im1)),nf,Parnam,Parmen,WW6))))
            Parnam=c(Parnam,paste0("la3_0_",as.character(im1),"_",as.character(im2)))
            Parmen=c(Parmen,eval(parse(text=subsi1c(paste0(lna,"^",as.character(im2),"*",ea1,"^",as.character(im1)),nf,Parnam,Parmen,WW6))))
          }
        }
      }
    }
    assign("la1_0_0_0",1)
    assign("la2_0_0_0",1)
    assign("la3_0_0_0",1)
    Parnam=c(Parnam,"la1_0_0_0","la2_0_0_0","la3_0_0_0")
    Parmen=c(Parmen,1,1,1)
    ea1nu=paste0("(1-",ea1,")^",nu)
    ea1ynu2=paste0("(1-",ea1,"*y_)^",nu2)
    ea1ynu3=paste0("(1-",ea1,"*y_)^",nu3)
    ea1ynu4=paste0("(1-",ea1,"*y_)^",nu4)
    ea2ynu3=paste0("(1-",ea1,"-",eaa1,"*y_)^(",nu3,"-k1)")
    ea2ynu3_=paste0("(1-",eaa1,"-",ea1,"*y_)^(",nu3,"-k2)")
    ea2ynu4=paste0("(1-",ea1,"-",eaa1,"*y_)^(",nu4,"-k1)")
    ea2ynu4_=paste0("(1-",eaa1,"-",ea1,"*y_)^(",nu4,"-k2)")
    ea2y_yynu4=paste0("(1-",ea1,"*y_-",eaa1,"*y__)^",nu4)
    ea2y_yynu5=paste0("(1-",ea1,"*y_-",eaa1,"*y__)^",nu5)
    ea3ynu4_1=paste0("(1-",ea1,"-",eaa1,"-",eaaa1,"*y_)^(",nu4,"-k1-k2)")
    ea3ynu4_2=paste0("(1-",ea1,"-",eaaa1,"-",eaa1,"*y_)^(",nu4,"-k1-k3)")
    ea3ynu4_3=paste0("(1-",eaa1,"-",eaaa1,"-",ea1,"*y_)^(",nu4,"-k2-k3)")
    ea3y_yynu5_1=paste0("(1-",ea1,"-",eaa1,"*y_-",eaaa1,"*y__)^(",nu5,"-k1)")
    ea3y_yynu5_2=paste0("(1-",eaa1,"-",ea1,"*y_-",eaaa1,"*y__)^(",nu5,"-k2)")
    ea3y_yynu5_3=paste0("(1-",eaaa1,"-",ea1,"*y_-",eaa1,"*y__)^(",nu5,"-k3)")
    ea1nuk=paste0("(1-",ea1,")^(",nu1,"-k)")
    ea2nuk=paste0("(1-",ea1,")^(",nu2,"-k)")
    ea3nuk=paste0("(1-",ea1,")^(",nu3,"-k)")
    ea4nuk=paste0("(1-",ea1,")^(",nu4,"-k)")
    ea2nuk1k2=paste0("(1-",ea1,"-",eaa1,")^(",nu2,"-k1-k2)")
    ea3nuk1k2=paste0("(1-",ea1,"-",eaa1,")^(",nu3,"-k1-k2)")
    ea3nuk1k2k3=paste0("(1-",ea1,"-",eaa1,"-",eaaa1,")^(",nu3,"-k1-k2-k3)")
    ea3nu6=paste0("(1-",ea1,"*y_-",eaa1,"*y__-",eaaa1,"*y___)^(",nu6,")")
    #Full begin
    I1=  paste0("(gamma(k+1)*gamma(",nu1,"-k)*(",ea1nuk,")^(-1)-gamma(",nu2,")*",ea1,"^(k+1)*(1/(k+1)+",nu2,"*",ea1,"/(k^2+3*k+2)+",nu2,"*",nu3,"*",ea1,"^2/((k+1)*(k+2)*(k+3))))*(b^pl*gamma(",nu,")*",nun,")^(-1)")
    I2_1=paste0("(gamma(k+1)*gamma(",nu2,"-k)*(",ea2nuk,"^(-1))-gamma(",nu3,")*",ea1,"^(k+1)*(1/(k+1)+",nu3,"*",ea1,"/(k^2+3*k+2)+",nu3,"*",nu4,"*",ea1,"^2/((k+1)*(k+2)*(k+3))))*(b^pl*gamma(",nu,")*",nun,"*",nun1,")^(-1)")
    I3_1=paste0("(gamma(k+1)*gamma(",nu3,"-k)*(",ea3nuk,"^(-1))-gamma(",nu4,")*",ea1,"^(k+1)*(1/(k+1)+",nu4,"*",ea1,"/(k^2+3*k+2)+",nu4,"*",nu5,"*",ea1,"^2/((k+1)*(k+2)*(k+3))))*(b^pl*gamma(",nu,")*",nun,"*",nun1,"*",nun2,")^(-1)")
    I2=  paste0("(gamma(k1+1)*gamma(k2+1)*gamma(",nu2,"-k1-k2)*(",ea2nuk1k2,")^(-1)-gamma(k1+1)*gamma(",nu3,"-k1)*",eaa1,"^(k2+1)*((1+(",nu3,"-k1)*",ea1,"+0.5*(",nu3,"-k1)*(",nu4,"-k1)*",ea1,"^2)/(k2+1)+((",nu3,"-k1)*",eaa1,"+(",nu3,"-k1)*(",nu4,"-k1)*",ea1,"*",eaa1,")/(k2^2+3*k2+2)+((",nu3,"-k1)*(",nu4,"-k1)*",eaa1,"^2)/((k2+1)*(k2+2)*(k2+3)))-gamma(k2+1)*gamma(",nu3,"-k2)*",ea1,"^(k1+1)*((1+(",nu3,"-k2)*",eaa1,"+0.5*(",nu3,"-k2)*(",nu4,"-k2)*",eaa1,"^2)/(k1+1)+((",nu3,"-k2)*",ea1,"+(",nu3,"-k2)*(",nu4,"-k2)*",eaa1,"*",ea1,")/(k1^2+3*k1+2)+((",nu3,"-k2)*(",nu4,"-k2)*",ea1,"^2)/((k1+1)*(k1+2)*(k1+3)))+gamma(",nu4,")*",ea1,"^(k1+1)*",eaa1,"^(k2+1)*(1/((k1+1)*(k2+1))+",nu4,"*",ea1,"/((k1^2+3*k1+2)*(k2+1))+",nu4,"*",eaa1,"/((k2^2+3*k2+2)*(k1+1))))*(b^pl*gamma(",nu,")*",nun,"*",nun1,")^(-1)")
    I3_2=paste0("(gamma(k1+1)*gamma(k2+1)*gamma(",nu3,"-k1-k2)*(",ea3nuk1k2,")^(-1)-gamma(k1+1)*gamma(",nu4,"-k1)*",eaa1,"^(k2+1)*((1+(",nu4,"-k1)*",ea1,"+0.5*(",nu4,"-k1)*(",nu5,"-k1)*",ea1,"^2)/(k2+1)+((",nu4,"-k1)*",eaa1,"+(",nu4,"-k1)*(",nu5,"-k1)*",ea1,"*",eaa1,")/(k2^2+3*k2+2)+((",nu4,"-k1)*(",nu5,"-k1)*",eaa1,"^2)/((k2+1)*(k2+2)*(k2+3)))-gamma(k2+1)*gamma(",nu4,"-k2)*",ea1,"^(k1+1)*((1+(",nu4,"-k2)*",eaa1,"+0.5*(",nu4,"-k2)*(",nu5,"-k2)*",eaa1,"^2)/(k1+1)+((",nu4,"-k2)*",ea1,"+(",nu4,"-k2)*(",nu5,"-k2)*",eaa1,"*",ea1,")/(k1^2+3*k1+2)+((",nu4,"-k2)*(",nu5,"-k2)*",ea1,"^2)/((k1+1)*(k1+2)*(k1+3)))+gamma(",nu5,")*",ea1,"^(k1+1)*",eaa1,"^(k2+1)*(1/((k1+1)*(k2+1))+",nu5,"*",ea1,"/((k1^2+3*k1+2)*(k2+1))+",nu5,"*",eaa1,"/((k2^2+3*k2+2)*(k1+1))))*(b^pl*gamma(",nu,")*",nun,"*",nun1,"*",nun2,")^(-1)")
    I3=paste0("(gamma(k1+1)*gamma(k2+1)*gamma(k3+1)*gamma(",nu3,"-k1-k2-k3)*(",ea3nuk1k2k3,")^(-1)-gamma(k1+1)*gamma(k2+1)*gamma(",nu4,"-k1-k2)*",eaaa1,"^(k3+1)*((1+(",ea1,"+",eaa1,")*(",nu4,"-k1-k2)+0.5*(",nu4,"-k1-k2)*(",nu5,"-k1-k2)*(",ea1,"+",eaa1,")^2)/(k3+1)+(",eaaa1,"*(",nu4,"-k1-k2)+(",ea1,"+",eaa1,")*",eaaa1,"*(",nu4,"-k1-k2)*(",nu5,"-k1-k2))/(k3^2+3*k3+2)+(",nu4,"-k1-k2)*(",nu5,"-k1-k2)*",eaaa1,"^2/((k3+1)*(k3+2)*(k3+3)))-gamma(k1+1)*gamma(k3+1)*gamma(",nu4,"-k1-k3)*",eaa1,"^(k2+1)*((1+(",ea1,"+",eaaa1,")*(",nu4,"-k1-k3)+0.5*(",nu4,"-k1-k3)*(",nu5,"-k1-k3)*(",ea1,"+",eaaa1,")^2)/(k2+1)+(",eaa1,"*(",nu4,"-k1-k3)+(",ea1,"+",eaaa1,")*",eaa1,"*(",nu4,"-k1-k3)*(",nu5,"-k1-k3))/(k2^2+3*k2+2)+(",nu4,"-k1-k3)*(",nu5,"-k1-k3)*",eaa1,"^2/((k2+1)*(k2+2)*(k2+3)))-gamma(k2+1)*gamma(k3+1)*gamma(",nu4,"-k2-k3)*",ea1,"^(k1+1)*((1+(",eaaa1,"+",eaa1,")*(",nu4,"-k3-k2)+0.5*(",nu4,"-k3-k2)*(",nu5,"-k3-k2)*(",eaaa1,"+",eaa1,")^2)/(k1+1)+(",ea1,"*(",nu4,"-k3-k2)+(",eaaa1,"+",eaa1,")*",ea1,"*(",nu4,"-k3-k2)*(",nu5,"-k3-k2))/(k1^2+3*k1+2)+(",nu4,"-k3-k2)*(",nu5,"-k3-k2)*",ea1,"^2/((k1+1)*(k1+2)*(k1+3)))+gamma(k1+1)*gamma(",nu5,"-k1)*",eaa1,"^(k2+1)*",eaaa1,"^(k3+1)*((1+(",nu5,"-k1)*",ea1,")/((k2+1)*(k3+1))+(",nu5,"-k1)*",eaa1,"/((k2^2+3*k2+2)*(k3+1))+(",nu5,"-k1)*",eaaa1,"/((k3^2+3*k3+2)*(k2+1)))+gamma(k2+1)*gamma(",nu5,"-k2)*",ea1,"^(k1+1)*",eaaa1,"^(k3+1)*((1+(",nu5,"-k2)*",eaa1,")/((k1+1)*(k3+1))+(",nu5,"-k2)*",ea1,"/((k1^2+3*k1+2)*(k3+1))+(",nu5,"-k2)*",eaaa1,"/((k3^2+3*k3+2)*(k1+1)))+gamma(k3+1)*gamma(",nu5,"-k3)*",ea1,"^(k1+1)*",eaa1,"^(k2+1)*((1+(",nu5,"-k3)*",eaaa1,")/((k2+1)*(k1+1))+(",nu5,"-k3)*",eaa1,"/((k2^2+3*k2+2)*(k1+1))+(",nu5,"-k3)*",ea1,"/((k1^2+3*k1+2)*(k2+1)))-gamma(",nu6,")*",ea1,"^(k1+1)*",eaa1,"^(k2+1)*",eaaa1,"^(k3+1))*(b^pl*gamma(",nu,")*",nun,"*",nun1,"*",nun2,")^(-1)")
    #Full end
    PartialI2_0=paste0(rep("gamma(k1+1)*gamma(k2+1)*gamma(mu2-k1-k2)*",10),rep("((b^pl*gamma(mu)*mun*mun1)^(-1))*",10),c(1,rep("((exp(-g)+2-k1-k2))",2),rep("0.5*((exp(-g)+2-k1-k2)*(exp(-g)+3-k1-k2))",2),rep("1*((exp(-g)+2-k1-k2)*(exp(-g)+3-k1-k2))",1),rep("(1/6)*((exp(-g)+2-k1-k2)*(exp(-g)+3-k1-k2)*(exp(-g)+4-k1-k2))",2),rep("0.5*((exp(-g)+2-k1-k2)*(exp(-g)+3-k1-k2)*(exp(-g)+4-k1-k2))",2)))
    PartialI2_1=c(paste0(rep("(-gamma(k1+1)*gamma(mu3-k1))*",3),rep("((b^pl*gamma(mu)*mun*mun1)^(-1))*",3),c("1/(k2+1)","(mu3-k1)/(k2+1)","(mu3-k1)/((k2+1)*(k2+2))")),paste0(rep("(-gamma(k2+1)*gamma(mu3-k2))*",3),rep("((b^pl*gamma(mu)*mun*mun1)^(-1))*",3),c("1/(k1+1)","(mu3-k2)/(k1+1)","(mu3-k2)/((k1+1)*(k1+2))")))
    PartialI3_0=paste0(rep("(gamma(k1+1)*gamma(k2+1)*gamma(k3+1)*gamma(exp(-g)+3-k1-k2-k3))*",20),rep("((b^pl*gamma(exp(-g))*(exp(-g)+n)*(exp(-g)+n+1)*(exp(-g)+n+2))^(-1))*",20),c("1",rep("((exp(-g)+3-k1-k2-k3))",3),rep("0.5*((exp(-g)+3-k1-k2-k3)*(exp(-g)+4-k1-k2-k3))",3),rep("1*((exp(-g)+3-k1-k2-k3)*(exp(-g)+4-k1-k2-k3))",3),rep("(1/6)*((exp(-g)+3-k1-k2-k3)*(exp(-g)+4-k1-k2-k3)*(exp(-g)+5-k1-k2-k3))",3),rep("0.5*((exp(-g)+3-k1-k2-k3)*(exp(-g)+4-k1-k2-k3)*(exp(-g)+5-k1-k2-k3))",6),"1*((exp(-g)+3-k1-k2-k3)*(exp(-g)+4-k1-k2-k3)*(exp(-g)+5-k1-k2-k3))"))
    PartialI3_1=c(paste0(rep("(-gamma(k1+1)*gamma(k2+1)*gamma(mu4-k1-k2))*",4),rep("((b^pl*gamma(exp(-g))*(exp(-g)+n)*(exp(-g)+n+1)*(exp(-g)+n+2))^(-1))*",4),c("1/(k3+1)",rep("((mu4-k1-k2)/(k3+1))",2),"((mu4-k1-k2)/((k3+1)*(k3+2)))")),paste0(rep("(-gamma(k1+1)*gamma(k3+1)*gamma(mu4-k1-k3))*",4),rep("((b^pl*gamma(exp(-g))*(exp(-g)+n)*(exp(-g)+n+1)*(exp(-g)+n+2))^(-1))*",4),c("1/(k2+1)",rep("((mu4-k1-k3)/(k2+1))",2),"((mu4-k1-k3)/((k2+1)*(k2+2)))")),paste0(rep("(-gamma(k3+1)*gamma(k2+1)*gamma(mu4-k3-k2))*",4),rep("((b^pl*gamma(exp(-g))*(exp(-g)+n)*(exp(-g)+n+1)*(exp(-g)+n+2))^(-1))*",4),c("1/(k1+1)",rep("((mu4-k3-k2)/(k1+1))",2),"((mu4-k3-k2)/((k1+1)*(k1+2)))")))
    AK1_0=rep(0,20)
    AK2_0=rep(0,20)
    AK3_0=rep(0,20)
    A1_0=c(0,1,0,0,2,0,0,1,1,0,3,0,0,2,1,2,1,0,0,1)#,4,0,0,3,1,3,1,0,0,2,2,0,2,1,1)
    A2_0=c(0,0,1,0,0,2,0,1,0,1,0,3,0,1,2,0,0,2,1,1)#,0,4,0,1,3,0,0,3,1,2,0,2,1,2,1)
    A3_0=c(0,0,0,1,0,0,2,0,1,1,0,0,3,0,0,1,2,1,2,1)#,0,0,4,0,0,1,3,1,3,0,2,2,1,1,2)
    LA1_0=rep(0,20)
    LA2_0=rep(0,20)
    LA3_0=rep(0,20)

    AK1_1=c(rep(0,4),rep(0,4),rep(1,4))
    AK2_1=c(rep(0,4),rep(1,4),rep(0,4))
    AK3_1=c(rep(1,4),rep(0,4),rep(0,4))
    A1_1=c(0,1,0,0,0,1,0,0,0,0,0,1)
    A2_1=c(0,0,1,0,0,0,0,1,0,0,1,0)
    A3_1=c(0,0,0,1,0,0,1,0,0,1,0,0)
    LA1_1=rep(0,12)
    LA2_1=rep(0,12)
    LA3_1=rep(0,12)

    III3<-sapply(data.frame(c(PartialI3_0,PartialI3_1),c(AK1_0,AK1_1),c(AK2_0,AK2_1),c(AK3_0,AK3_1),c(A1_0,A1_1),c(A2_0,A2_1),c(A3_0,A3_1),c(LA1_0,LA1_1),c(LA2_0,LA2_1),c(LA3_0,LA3_1)),as.character)

    AAK1_0=rep(0,10)
    AAK2_0=rep(0,10)
    AAK3_0=rep(0,10)
    AA1_0=c(0,1,0,2,0,1,3,0,2,1)
    AA2_0=c(0,0,1,0,2,1,0,3,1,2)
    AA3_0=c(0,0,0,0,0,0,0,0,0,0)
    LAA1_0=rep(0,10)
    LAA2_0=rep(0,10)
    LAA3_0=rep(0,10)

    AAK1_1=c(rep(0,3),rep(1,3))
    AAK2_1=c(rep(1,3),rep(0,3))
    AAK3_1=c(rep(0,3),rep(0,3))
    AA1_1=c(0,1,0,0,0,1)
    AA2_1=c(0,0,1,0,1,0)
    AA3_1=c(0,0,0,0,0,0)
    LAA1_1=rep(0,6)
    LAA2_1=rep(0,6)
    LAA3_1=rep(0,6)

    Part2_2=expression(-d/b^2)
    Part2_3=expression(2*d/b^3)
    Part3=expression(lgamma(exp(-g)+d)-lgamma(exp(-g)))
    Part3_2=as.expression(DD(Part3,"g",2))
    Part3_3=as.expression(DD(Part3,"g",3))
    ##########################################
    #Second order derivative
    ##########################################
    assign("pr",0)
    cc=array(expression(nul),c(np,np))
    cc[[np,np]]=eval(substitute(Part2_2,list(Part2_2=Part2_2[[1]])))
    cc_c=array(expression(nul),c(np,np,np))
    cc_c[[np,np,np]]=eval(Part2_3)

    Rcc=array(0,c(np,np))
    Rcc_c=array(0,c(np,np,np))
    ###########################F2.1
    sumug=u[nf+2]
    work=expression(exp(-g)*sumug/(exp(-g)+n))
    work=dpg(as.expression(eval(substitute(substitute(e, list(sumug=sumug)), list(e = work[[1]])))))
    cc[[(nf+2),(nf+2)]]=substitute(Part3_2+2*work,list(Part3_2=Part3_2[[1]],work=work[[1]]))

    for (i in 1:(np-2)){
      sumui=u[i]
      work=expression(exp(-g)*sumui/(exp(-g)+n))
      work=dpg(as.expression(eval(substitute(substitute(e, list(sumui=sumui)), list(e = work[[1]])))))
      cc[[(nf+2),i]]=substitute(work,list(work=work[[1]]))
      cc[[i,(nf+2)]]=cc[[(nf+2),i]]
    }
    E1=paste0(nu,"*",I1)
    E1k=paste0(deparse(DD(parse(text=E1),"k",1)),collapse="")
    E0=paste0(E1k,"-",lna,"*",E1)
    k<-1
    pl<-1
    c1<-1
    Parmen[which(Parnam=="k")]=k
    Parmen[which(Parnam=="pl")]=pl
    Parmen[which(Parnam=="c1")]=c1
    work=subsi1c(E0,nf,Parnam,Parmen,WW6)
    work1=as.expression(cc[[(nf+2),(nf+3)]])
    cc[[(nf+2),(nf+3)]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
    cc[[(nf+3),(nf+2)]]=cc[[(nf+2),(nf+3)]]

    for (i in 1:np){
      E0_i=paste0(deparse(DD(parse(text=E0),parnam[i],1)),collapse="")
      work1=as.expression(cc_c[[(nf+2),(nf+3),i]])
      work=subsi1c(E0_i,nf,Parnam,Parmen,WW6)
      cc_c[[(nf+2),(nf+3),i]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
      cc_c[[(nf+3),(nf+2),i]]=cc_c[[(nf+2),(nf+3),i]]
    }
    ###########################F2.2
    for (rj1 in 1:(np-1)){
      E1=paste0(nud,"*x_",as.character(rj1),"*",I1)
      E1k=paste0(deparse(DD(parse(text=E1),"k",1)),collapse="")
      E0=paste0(E1k,"-",lna,"*",E1)
      k<-1
      pl<-1
      c1<-1
      Parmen[which(Parnam=="k")]=k
      Parmen[which(Parnam=="pl")]=pl
      Parmen[which(Parnam=="c1")]=c1

      work=subsi1c(E0,nf,Parnam,Parmen,WW6)
      work1=as.expression(cc[[rj1,np]])
      cc[[rj1,np]]=substitute(-work+work1,list(work=work[[1]],work1=work1[[1]]))
      work1=as.expression(cc[[np,rj1]])
      cc[[np,rj1]]=substitute(-work+work1,list(work=work[[1]],work1=work1[[1]]))
      for (rj2 in 1:(np-1)){
        u12=uu[[rj1,rj2]]
        work=expression(-(exp(-g)+d)*u12/(exp(-g)+n))
        work=dpg(as.expression(eval(substitute(substitute(e, list(u12=u12)), list(e = work[[1]])))))
        work1=as.expression(cc[[rj1,rj2]])
        cc[[rj1,rj2]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
      }
      for (i in 1:np){
        E0_i=paste0(deparse(DD(parse(text=E0),parnam[i],1)),collapse="")
        work=subsi1c(E0_i,nf,Parnam,Parmen,WW6)
        work1=as.expression(cc_c[[rj1,np,i]])
        cc_c[[rj1,np,i]]=substitute(-work+work1,list(work=work[[1]],work1=work1[[1]]))
        cc_c[[np,rj1,i]]=cc_c[[rj1,np,i]]
      }
    }
    E1=paste0(nud,"*",I1)
    E1k=paste0(deparse(DD(parse(text=E1),"k",1)),collapse="")
    E2k=paste0(deparse(DD(parse(text=E1),"k",2)),collapse="")
    E0=paste0(E2k,"-2*",lna,"*",E1k,"+",lna,"^2*",E1)
    k<-1
    pl<-2
    c1<-1
    Parmen[which(Parnam=="k")]=k
    Parmen[which(Parnam=="pl")]=pl
    Parmen[which(Parnam=="c1")]=c1
    work=subsi1c(E0,nf,Parnam,Parmen,WW6)
    work1=eval(as.expression(cc[[np,np]]))
    cc[[np,np]]=substitute(-work+work1,list(work=work[[1]],work1=work1[[1]]))

    for (i in 1:np){
      E0_i=paste0(deparse(DD(parse(text=E0),parnam[i],1)),collapse="")
      work=subsi1c(E0_i,nf,Parnam,Parmen,WW6)
      work1=as.expression(cc_c[[np,np,i]])
      cc_c[[np,np,i]]=substitute(-work+work1,list(work=work[[1]],work1=work1[[1]]))
    }
    ###########################F2.3
    for (rj1 in 1:(np-1)){
      u1=u[rj1]
      E2_1=paste0(nud,"*x__",as.character(rj1),"*",I2)
      E2_1k1=paste0(deparse(DD(parse(text=E2_1),"k1",1)),collapse="")
      E2z=paste0(E2_1k1,"-",lna,"*",E2_1)
      k1<-1
      k2<-1
      pl<-1
      c1<-1
      c2<-0
      Parmen[which(Parnam=="k1")]=k1
      Parmen[which(Parnam=="k2")]=k2
      Parmen[which(Parnam=="pl")]=pl
      Parmen[which(Parnam=="c1")]=c1
      Parmen[which(Parnam=="c2")]=c2
      E2=paste0(nud,"*x_",as.character(rj1),"*",I2_1)
      E2k=paste0(deparse(DD(parse(text=E2),"k",1)),collapse="")
      E20=gsub("\n","",paste0(E2k,"-",lna,"*",E2))
      k<-2
      Parmen[which(Parnam=="k")]=k
      work=subsi2c(E2z,nf,Parnam,Parmen ,WW4,WW5)-subsi2c_i(E2z,nf,Parnam,Parmen,WW6)+subsi1c(E20,nf,Parnam,Parmen,WW6)
      E2_=paste0(nud,"*x_",as.character(rj1),"*",ea1,"*",I2_1)
      E2k_=paste0(deparse(DD(parse(text=E2_),"k",1)),collapse="")
      E20_=gsub("\n","",paste0(E2k_,"-",lna,"*",E2_))
      k<-1
      c1<-1
      Parmen[which(Parnam=="k")]=k
      Parmen[which(Parnam=="c1")]=c1
      work=work-subsi1c(E20_,nf,Parnam,Parmen,WW6)
      work1=as.expression(cc[[rj1,np]])
      cc[[rj1,np]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
      work1=as.expression(cc[[np,rj1]])
      cc[[np,rj1]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
      for (i in 1:np){
        E2z_i=paste0(deparse(DD(parse(text=E2z),parnam[i],1)),collapse="")
        E20_i=paste0(deparse(DD(parse(text=E20),parnam[i],1)),collapse="")
        k<-2
        Parmen[which(Parnam=="k")]=k
        work=subsi2c(E2z_i,nf,Parnam,Parmen,WW4,WW5)-subsi2c_i(E2z_i,nf,Parnam,Parmen,WW6)+subsi1c(E20_i,nf,Parnam,Parmen,WW6)
        E20_i_=paste0(deparse(DD(parse(text=E20_),parnam[i],1)),collapse="")
        k<-1
        c1<-1
        Parmen[which(Parnam=="k")]=k
        Parmen[which(Parnam=="c1")]=c1
        work=work-subsi1c(E20_i_,nf,Parnam,Parmen,WW6)
        work1=as.expression(cc_c[[rj1,np,i]])
        cc_c[[rj1,np,i]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
        work1=as.expression(cc_c[[np,rj1,i]])
        cc_c[[np,rj1,i]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
      }

      for (rj2 in 1:(np-1)){
        u2=u[rj2]
        u12=uu[rj1,rj2]
        work=expression(((exp(-g)+d)/((exp(-g)+n)*(exp(-g)+n+1)))*(u12+u1*u2))
        work=dpg(as.expression(eval(substitute(substitute(e, list(u1=u1,u2=u2,u12=u12)), list(e = work[[1]])))))
        work1=as.expression(cc[[rj1,rj2]])
        cc[[rj1,rj2]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
      }
    }
    E2=paste0(nud,"*",I2)
    E2k1k2=gsub("\n","",paste0(deparse(DD(DD(parse(text=E2),"k1",1),"k2",1)),collapse=""))
    k1<-1
    k2<-1
    pl<-2
    c1<-1
    c2<-1
    Parmen[which(Parnam=="k1")]=k1
    Parmen[which(Parnam=="c1")]=c1
    Parmen[which(Parnam=="k2")]=k2
    Parmen[which(Parnam=="c2")]=c2
    Parmen[which(Parnam=="pl")]=pl
    E2_2=paste0(nud,"*",lnaa,"*",I2)
    E2_2k1=paste0(deparse(DD(parse(text=E2_2),"k1",1)),collapse="")
    E2_3=paste0(nud,"*",lna,"*",lnaa,"*",I2)
    E2z=paste0(E2k1k2,"-2*",E2_2k1,"+",E2_3)
    work=subsi2c(E2z,nf,Parnam,Parmen,WW4,WW5)-subsi2c_i(E2z,nf,Parnam,Parmen,WW6)
    E2_0=paste0(nud,"*",lna,"^2*",I2_1)
    E2_1=paste0(nud,"*",lna,"*",I2_1)
    E2_2=paste0(nud,"*",I2_1)
    E2_1k1=paste0(deparse(DD(parse(text=E2_1),"k",1)),collapse="")
    E2k2=paste0(deparse(DD(parse(text=E2_2),"k",2)),collapse="")
    E2z_1=paste0(E2k2,"-2*",E2_1k1,"+",E2_0)
    k<-2
    Parmen[which(Parnam=="k")]=k
    work=work+subsi1c(E2z_1,nf,Parnam,Parmen,WW6)
    work1=eval(as.expression(cc[[np,np]]))
    cc[[np,np]]=work+work1
    for (i in 1:np){
      E2z_i=gsub("\n","",paste0(deparse(DD(parse(text=E2z),parnam[i],1)),collapse=""))
      E2z_1_i=gsub("\n","",paste0(deparse(DD(parse(text=E2z_1),parnam[i],1)),collapse=""))
      work=subsi2c(E2z_i,nf,Parnam,Parmen,WW4,WW5)-subsi2c_i(E2z_i,nf,Parnam,Parmen,WW6)
      work=work+subsi1c(E2z_1_i,nf,Parnam,Parmen,WW6)
      work1=as.expression(cc_c[[np,np,i]])
      cc_c[[np,np,i]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
    }
    Rcc=array(NA,c(np,np))
    Rcc_c=array(NA,c(np,np,np))
    for (i1 in 1:np){
      for (i2 in 1:np){
        Rcc[i1,i2]=eval(cc[[i1,i2]])
        if (i1!=np & i2!=np){
          for (i3 in 1:np){
            cc_c[[i1,i2,i3]]=DD(cc[[i1,i2]],parnam[i3],1)
            Rcc_c[[i1,i2,i3]]=eval(DD(cc[[i1,i2]],parnam[i3],1))}}
        else {
          for (i3 in 1:np){
            Rcc_c[[i1,i2,i3]]=eval(cc_c[[i1,i2,i3]])}
        }
      }
    }
    ###########################F2.0
    Rcc[nf+2,nf+2]=Rcc[nf+2,nf+2]-exp(-g)*sum(1/(exp(-g)+(0:(n-1))))
    Rcc_c[nf+2,nf+2,nf+2]=Rcc_c[nf+2,nf+2,nf+2]+exp(-g)*sum(1/(exp(-g)+(0:(n-1))))-exp(-2*g)*sum(1/((exp(-g)+(0:(n-1)))^2))
    #    Ti=print(c(as.character(round(as.numeric((proc.time() - ptm))[1],2))))
    ##########################################
      #Third order derivative
    ###########################F3.1
    ccc=array(expression(nul),c(np,np,np))
    ccc[[np,np,np]]=eval(substitute(Part2_3,list(Part2_3=Part2_3[[1]])))
    sumug=u[nf+2]
    work=expression(-exp(-g)*sumug/(exp(-g)+n))
    work=dpg(as.expression(eval(substitute(substitute(e, list(sumug=sumug)), list(e = work[[1]])))))
    ccc[[(np-1),(np-1),(np-1)]]=substitute(Part3_3+3*work,list(Part3_3=Part3_3[[1]],work=work[[1]]))
    for (i in 1:(np-2)){
      sumui=u[i]
      work=expression(-exp(-g)*sumui/(exp(-g)+n))
      work=dpg(as.expression(eval(substitute(substitute(e, list(sumui=sumui)), list(e = work[[1]])))))
      ccc[[(np-1),(np-1),i]]=substitute(work,list(work=work[[1]]))
      ccc[[i,(np-1),(np-1)]]=ccc[[(np-1),(np-1),i]]
      ccc[[(np-1),i,(np-1)]]=ccc[[(np-1),(np-1),i]]
    }

    E1=paste0(nu,"*",I1)
    E1k=paste0(deparse(DD(parse(text=E1),"k",1)),collapse="")
    E0=paste0(E1k,"-",lna,"*",E1)
    k<-1
    pl<-1
    c1<-1
    Parmen[which(Parnam=="k")]=k
    Parmen[which(Parnam=="pl")]=pl
    Parmen[which(Parnam=="c1")]=c1
    work=subsi1c(E0,nf,Parnam,Parmen,WW6)
    ccc[[(np-1),(np-1),np]]=substitute(-work,list(work=work[[1]]))
    ccc[[(np-1),np,(np-1)]]=substitute(-work,list(work=work[[1]]))
    ccc[[np,(np-1),(np-1)]]=substitute(-work,list(work=work[[1]]))
    ###########################F3.2
    for (rj1 in 1:(np-1)){
      u1=u[rj1]
      E1=paste0(nu,"*x_",as.character(rj1),"*",I1)
      E1k=paste0(deparse(DD(parse(text=E1),"k",1)),collapse="")
      E0=paste0(E1k,"-",lna,"*",E1)
      k<-1
      pl<-1
      c1<-1
      Parmen[which(Parnam=="k")]=k
      Parmen[which(Parnam=="pl")]=pl
      Parmen[which(Parnam=="c1")]=c1
      work=subsi1c(E0,nf,Parnam,Parmen,WW6)
      work1=as.expression(ccc[[rj1,np,np-1]])
      ccc[[rj1,np,np-1]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
      work1=as.expression(ccc[[np,rj1,np-1]])
      ccc[[np,rj1,np-1]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
      work1=as.expression(ccc[[rj1,np-1,np]])
      ccc[[rj1,np-1,np]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
      work1=as.expression(ccc[[np,np-1,rj1]])
      ccc[[np,np-1,rj1]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
      work1=as.expression(ccc[[np-1,rj1,np]])
      ccc[[np-1,rj1,np]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
      work1=as.expression(ccc[[np-1,np,rj1]])
      ccc[[np-1,np,rj1]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
      for (rj2 in 1:(np-1)){
        u12=uu[rj1,rj2]
        work=expression(exp(-g)*u12/(exp(-g)+n))
        work=dpg(as.expression(eval(substitute(substitute(e, list(u12=u12)), list(e = work[[1]])))))
        work1=as.expression(ccc[[np-1,rj1,rj2]])
        ccc[[np-1,rj1,rj2]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
        work1=as.expression(ccc[[rj1,np-1,rj2]])
        ccc[[rj1,np-1,rj2]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
        work1=as.expression(ccc[[rj1,rj2,np-1]])
        ccc[[rj1,rj2,np-1]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
      }
    }
    E1=paste0(nu,"*",I1)
    E1k=paste0(deparse(DD(parse(text=E1),"k",1)),collapse="")
    E2k=paste0(deparse(DD(parse(text=E1),"k",2)),collapse="")
    E0=paste0(E2k,"-2*",lna,"*",E1k,"+",E1,"*",lna,"^2")
    k<-1
    pl<-2
    c1<-1
    Parmen[which(Parnam=="k")]=k
    Parmen[which(Parnam=="pl")]=pl
    Parmen[which(Parnam=="c1")]=c1
    work=subsi1c(E0,nf,Parnam,Parmen,WW6)
    work1=as.expression(ccc[[np-1,np,np]])
    ccc[[np-1,np,np]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
    ccc[[np,np-1,np]]=ccc[[np-1,np,np]]
    ccc[[np,np,np-1]]=ccc[[np-1,np,np]]
    ###########################F3.3
    for (rj1 in 1:(np-1)){
      u1=u[rj1]
      E2_1=paste0(nu,"*x__",as.character(rj1),"*",I2)
      E2_1k1=paste0(deparse(DD(parse(text=E2_1),"k1",1)),collapse="")
      E2z=paste0(E2_1k1,"-",E2_1,"*",lna)
      k1<-1
      pl<-1
      c1<-1
      k2<-1
      c2<-0
      Parmen[which(Parnam=="k1")]=k1
      Parmen[which(Parnam=="pl")]=pl
      Parmen[which(Parnam=="c1")]=c1
      Parmen[which(Parnam=="k2")]=k2
      Parmen[which(Parnam=="c2")]=c2

      E2=paste0(nu,"*x_",as.character(rj1),"*",I2_1)
      E2k=paste0(deparse(DD(parse(text=E2),"k",1)),collapse="")
      E20=gsub("\n","",paste0(E2k,"-",lna,"*",E2))
      k<-2
      Parmen[which(Parnam=="k")]=k
      work=subsi2c(E2z,nf,Parnam,Parmen,WW4,WW5)-subsi2c_i(E2z,nf,Parnam,Parmen,WW6)+subsi1c(E20,nf,Parnam,Parmen,WW6)
      E2=paste0(nu,"*x_",as.character(rj1),"*",ea1,"*",I2_1)
      E2k=paste0(deparse(DD(parse(text=E2),"k",1)),collapse="")
      E20=gsub("\n","",paste0(E2k,"-",lna,"*",E2))
      k<-1
      c1<-1
      Parmen[which(Parnam=="k")]=k
      Parmen[which(Parnam=="c1")]=c1

      work=work-subsi1c(E20,nf,Parnam,Parmen,WW6)
      work1=as.expression(ccc[[rj1,np,np-1]])
      ccc[[rj1,np,np-1]]=substitute(-work+work1,list(work=work[[1]],work1=work1[[1]]))
      work1=as.expression(ccc[[np,rj1,np-1]])
      ccc[[np,rj1,np-1]]=substitute(-work+work1,list(work=work[[1]],work1=work1[[1]]))
      work1=as.expression(ccc[[rj1,np-1,np]])
      ccc[[rj1,np-1,np]]=substitute(-work+work1,list(work=work[[1]],work1=work1[[1]]))
      work1=as.expression(ccc[[np,np-1,rj1]])
      ccc[[np,np-1,rj1]]=substitute(-work+work1,list(work=work[[1]],work1=work1[[1]]))
      work1=as.expression(ccc[[np-1,rj1,np]])
      ccc[[np-1,rj1,np]]=substitute(-work+work1,list(work=work[[1]],work1=work1[[1]]))
      work1=as.expression(ccc[[np-1,np,rj1]])
      ccc[[np-1,np,rj1]]=substitute(-work+work1,list(work=work[[1]],work1=work1[[1]]))
      for (rj2 in 1:(np-1)){
        u2=u[rj2]
        u12=uu[rj1,rj2]
        work=expression((-exp(-g)/((exp(-g)+n)*(exp(-g)+n+1)))*(u12+u1*u2))
        work=dpg(as.expression(eval(substitute(substitute(e, list(u1=u1,u2=u2,u12=u12)), list(e = work[[1]])))))
        work1=as.expression(ccc[[rj1,rj2,np-1]])
        ccc[[rj1,rj2,np-1]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
        work1=as.expression(ccc[[rj1,np-1,rj2]])
        ccc[[rj1,np-1,rj2]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
        work1=as.expression(ccc[[np-1,rj1,rj2]])
        ccc[[np-1,rj1,rj2]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
      }
    }
    E2=paste0(nu,"*",I2)
    E2k1k2=gsub("\n","",paste0(deparse(DD(DD(parse(text=E2),"k1",1),"k2",1)),collapse=""))
    k1<-1
    pl<-2
    c1<-1
    k2<-1
    c2<-1
    Parmen[which(Parnam=="k1")]=k1
    Parmen[which(Parnam=="pl")]=pl
    Parmen[which(Parnam=="c1")]=c1
    Parmen[which(Parnam=="k2")]=k2
    Parmen[which(Parnam=="c2")]=c2

    E2_2=paste0(nu,"*",lnaa,"*",I2)
    E2_2k1=paste0(deparse(DD(parse(text=E2_2),"k1",1)),collapse="")
    E2_3=paste0(nu,"*",lna,"*",lnaa,"*",I2)
    E2z=paste0(E2k1k2,"-2*",E2_2k1,"+",E2_3)
    work=subsi2c(E2z,nf,Parnam,Parmen,WW4,WW5)-subsi2c_i(E2z,nf,Parnam,Parmen,WW6)

    E2_0=paste0(nu,"*",lna,"^2*",I2_1)
    E2_1=paste0(nu,"*",lna,"*",I2_1)
    E2_2=paste0(nu,"*",I2_1)
    E2_1k1=paste0(deparse(DD(parse(text=E2_1),"k",1)),collapse="")
    E2k2=paste0(deparse(DD(parse(text=E2_2),"k",2)),collapse="")
    E2z_1=paste0(E2k2,"-2*",E2_1k1,"+",E2_0)
    k<-2
    Parmen[which(Parnam=="k")]=k
    work=work+subsi1c(E2z_1,nf,Parnam,Parmen,WW6)
    work1=as.expression(ccc[[np,np,np-1]])
    ccc[[np,np,np-1]]=substitute(-work+work1,list(work=work[[1]],work1=work1[[1]]))
    work1=as.expression(ccc[[np,np-1,np]])
    ccc[[np,np-1,np]]=substitute(-work+work1,list(work=work[[1]],work1=work1[[1]]))
    work1=as.expression(ccc[[np-1,np,np]])
    ccc[[np-1,np,np]]=substitute(-work+work1,list(work=work[[1]],work1=work1[[1]]))
    ###########################F3.4
    for (rj1 in 1:(np-1)){
      u1=u[rj1]

      E0=paste0(nud,"*x_",as.character(rj1),"*",lna,"^2*",I1)
      E1=paste0(nud,"*x_",as.character(rj1),"*",lna,"*",I1)
      E2=paste0(nud,"*x_",as.character(rj1),"*",I1)
      E1_1k=paste0(deparse(DD(parse(text=E1),"k",1)),collapse="")
      E2_2k=paste0(deparse(DD(parse(text=E2),"k",2)),collapse="")
      Ez=paste0(E2_2k,"-2*",E1_1k,"+",E0)
      k<-1
      pl<-2
      c1<-1
      Parmen[which(Parnam=="k")]=k
      Parmen[which(Parnam=="pl")]=pl
      Parmen[which(Parnam=="c1")]=c1

      work=subsi1c(Ez,nf,Parnam,Parmen,WW6)
      work1=as.expression(ccc[[rj1,np,np]])
      ccc[[rj1,np,np]]=substitute(-work+work1,list(work=work[[1]],work1=work1[[1]]))
      work1=as.expression(ccc[[np,rj1,np]])
      ccc[[np,rj1,np]]=substitute(-work+work1,list(work=work[[1]],work1=work1[[1]]))
      work1=as.expression(ccc[[np,np,rj1]])
      ccc[[np,np,rj1]]=substitute(-work+work1,list(work=work[[1]],work1=work1[[1]]))
      for (rj2 in 1:(np-1)){
        uu12=uu[rj1,rj2]

        E1=paste0(nud,"*x_",as.character(rj1),"*x_",as.character(rj2),"*",I1)
        E0=paste0(E1,"*",lna)
        E1_1k=paste0(deparse(DD(parse(text=E1),"k",1)),collapse="")
        Ez=paste0(E1_1k,"-",E0)
        k<-1
        pl<-1
        c1<-1
        Parmen[which(Parnam=="k")]=k
        Parmen[which(Parnam=="pl")]=pl
        Parmen[which(Parnam=="c1")]=c1

        work=subsi1c(Ez,nf,Parnam,Parmen,WW6)
        work1=as.expression(ccc[[rj1,rj2,np]])
        ccc[[rj1,rj2,np]]=substitute(-work+work1,list(work=work[[1]],work1=work1[[1]]))
        work1=as.expression(ccc[[rj1,np,rj2]])
        ccc[[rj1,np,rj2]]=substitute(-work+work1,list(work=work[[1]],work1=work1[[1]]))
        work1=as.expression(ccc[[np,rj1,rj2]])
        ccc[[np,rj1,rj2]]=substitute(-work+work1,list(work=work[[1]],work1=work1[[1]]))
        for (rj3 in 1:(np-1)){
          uuu123=uuu[rj1,rj2,rj3]
          work=expression(-(exp(-g)+d)*uuu123/(exp(-g)+n))
          work=dpg(as.expression(eval(substitute(substitute(e, list(uuu123=uuu123)), list(e = work[[1]])))))
          work1=as.expression(ccc[[rj1,rj2,rj3]])
          ccc[[rj1,rj2,rj3]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
        }
      }
    }
    E0=paste0(nud,"*",lna,"^3*",I1)
    E1=paste0(nud,"*",lna,"^2*",I1)
    E2=paste0(nud,"*",lna,"*",I1)
    E3=paste0(nud,"*",I1)
    E1_1k=paste0(deparse(DD(parse(text=E1),"k",1)),collapse="")
    E2_2k=paste0(deparse(DD(parse(text=E2),"k",2)),collapse="")
    E3_3k=paste0(deparse(DD(parse(text=E3),"k",3)),collapse="")
    Ez=paste0(E3_3k,"-3*",E2_2k,"+3*",E1_1k,"-",E0)
    k<-1
    pl<-3
    c1<-1
    Parmen[which(Parnam=="k")]=k
    Parmen[which(Parnam=="pl")]=pl
    Parmen[which(Parnam=="c1")]=c1
    work=subsi1c(Ez,nf,Parnam,Parmen,WW6)
    work1=as.expression(ccc[[np,np,np]])
    ccc[[np,np,np]]=substitute(-work+work1,list(work=work[[1]],work1=work1[[1]]))
####################################F3.5
    for (rj1 in 1:(np-1)){
      u1=u[rj1]
      E2_1=paste0(nud,"*x_",as.character(rj1),"*",I2)
      E2=paste0(E2_1,"*",lna,"*",lnaa)
      E2_1k1k2=gsub("\n","",paste0(deparse(DD(DD(parse(text=E2_1),"k1",1),"k2",1)),collapse=""))
      E2_2=paste0(E2_1,"*",lna)
      E2_2k2=gsub("\n","",paste0(deparse(DD(parse(text=E2_2),"k2",1)),collapse=""))
      E2_3=paste0(E2_1,"*",lnaa)
      E2_3k1=gsub("\n","",paste0(deparse(DD(parse(text=E2_3),"k1",1)),collapse=""))
      k1<-1
      pl<-2
      c1<-1
      k2<-1
      c2<-1
      Parmen[which(Parnam=="k1")]=k1
      Parmen[which(Parnam=="pl")]=pl
      Parmen[which(Parnam=="c1")]=c1
      Parmen[which(Parnam=="k2")]=k2
      Parmen[which(Parnam=="c2")]=c2
      E2z=paste0(E2_1k1k2,"+",E2,"-",E2_2k2,"-",E2_3k1)
      work=2*subsi2c(E2z,nf,Parnam,Parmen,WW4,WW5)-2*subsi2c_i(E2z,nf,Parnam,Parmen,WW6)
      E2_1=paste0(nud,"*x_",as.character(rj1),"*",I2)
      E2=paste0(E2_1,"*",lnaa,"^2")
      E2_2=paste0(E2_1,"*",lnaa)
      E2_2k2=gsub("\n","",paste0(deparse(DD(parse(text=E2_2),"k2",1)),collapse=""))
      E2_1k22=gsub("\n","",paste0(deparse(DD(parse(text=E2_1),"k2",2)),collapse=""))
      E2z=paste0(E2_1k22,"+",E2,"-2*",E2_2k2)
      c1<-0
      c2<-1
      Parmen[which(Parnam=="c1")]=c1
      Parmen[which(Parnam=="c2")]=c2
      work=work+subsi2c(E2z,nf,Parnam,Parmen,WW4,WW5)-subsi2c_i(E2z,nf,Parnam,Parmen,WW6)
      E2_2=paste0(nud,"*x_",as.character(rj1),"*",I2_1)
      E2_0=paste0(E2_2,"*",lna,"^2")
      E2_1=paste0(E2_2,"*",lna)
      E2_1k=paste0(deparse(DD(parse(text=E2_1),"k",1)),collapse="")
      E2_2k=paste0(deparse(DD(parse(text=E2_2),"k",2)),collapse="")
      E2z_1=paste0(E2_2k,"+",E2_0,"-2*",E2_1k)
      k<-2
      c1<-1
      Parmen[which(Parnam=="k")]=k
      Parmen[which(Parnam=="c1")]=c1
      work=work+3*subsi1c(E2z_1,nf,Parnam,Parmen,WW6)
      E2_2_=paste0(nud,"*x_",as.character(rj1),"*",ea1,"*",I2_1)
      E2_0_=paste0(E2_2_,"*",lna,"^2")
      E2_1_=paste0(E2_2_,"*",lna)
      E2_1k_=paste0(deparse(DD(parse(text=E2_1_),"k",1)),collapse="")
      E2_2k_=paste0(deparse(DD(parse(text=E2_2_),"k",2)),collapse="")
      E2z_1_=paste0(E2_2k_,"+",E2_0_,"-2*",E2_1k_)
      k<-1
      c1<-1
      Parmen[which(Parnam=="k")]=k
      Parmen[which(Parnam=="c1")]=c1
      work=work-3*subsi1c(E2z_1_,nf,Parnam,Parmen,WW6)
      work1=as.expression(ccc[[rj1,np,np]])
      ccc[[rj1,np,np]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
      work1=as.expression(ccc[[np,rj1,np]])
      ccc[[np,rj1,np]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
      work1=as.expression(ccc[[np,np,rj1]])
      ccc[[np,np,rj1]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
      for (rj2 in 1:(np-1)){
        u2=u[rj2]
        u12=uu[rj1,rj2]

        E2_1=paste0(nud,"*x_",as.character(rj1),"*x__",as.character(rj2),"*",I2)
        E2=paste0(E2_1,"*",lna)
        E2k1=gsub("\n","",paste0(deparse(DD(parse(text=E2_1),"k1",1)),collapse=""))
        k1<-1
        pl<-1
        c1<-1
        k2<-1
        c2<-0
        Parmen[which(Parnam=="k1")]=k1
        Parmen[which(Parnam=="pl")]=pl
        Parmen[which(Parnam=="c1")]=c1
        Parmen[which(Parnam=="k2")]=k2
        Parmen[which(Parnam=="c2")]=c2
        E2z=paste0(E2k1,"-",E2)
        work=subsi2c(E2z,nf,Parnam,Parmen,WW4,WW5)-subsi2c_i(E2z,nf,Parnam,Parmen,WW6)
        E2_1=paste0(nud,"*x_",as.character(rj1),"*x_",as.character(rj2),"*",I2)
        E2=paste0(E2_1,"*",lnaa)
        E2k1=gsub("\n","",paste0(deparse(DD(parse(text=E2_1),"k2",1)),collapse=""))
        k1<-1
        c1<-0
        k2<-1
        c2<-1
        Parmen[which(Parnam=="k1")]=k1
        Parmen[which(Parnam=="c1")]=c1
        Parmen[which(Parnam=="k2")]=k2
        Parmen[which(Parnam=="c2")]=c2

        E2z=paste0(E2k1,"-",E2)
        work=work+subsi2c(E2z,nf,Parnam,Parmen,WW4,WW5)-subsi2c_i(E2z,nf,Parnam,Parmen,WW6)
        E2_1=paste0(nud,"*x_",as.character(rj1),"*x__",as.character(rj2),"*",I2)
        E2=paste0(E2_1,"*",lnaa)
        E2k1=gsub("\n","",paste0(deparse(DD(parse(text=E2_1),"k2",1)),collapse=""))
        k1<-1
        c1<-0
        k2<-1
        c2<-1
        Parmen[which(Parnam=="k1")]=k1
        Parmen[which(Parnam=="c1")]=c1
        Parmen[which(Parnam=="k2")]=k2
        Parmen[which(Parnam=="c2")]=c2

        E2z=paste0(E2k1,"-",E2)
        work=work+subsi2c(E2z,nf,Parnam,Parmen,WW4,WW5)-subsi2c_i(E2z,nf,Parnam,Parmen,WW6)
        E2_1=paste0(nud,"*x_",as.character(rj1),"*x_",as.character(rj2),"*",I2_1)
        E2=paste0(E2_1,"*",lna)
        E2_1k1=paste0(deparse(DD(parse(text=E2_1),"k",1)),collapse="")
        E2z_1=paste0(E2_1k1,"-",E2)
        k<-2
        c1<-1
        Parmen[which(Parnam=="k")]=k
        Parmen[which(Parnam=="c1")]=c1
        work=work+3*subsi1c(E2z_1,nf,Parnam,Parmen,WW6)

        E2_1=paste0(nud,"*x_",as.character(rj1),"*x_",as.character(rj2),"*",ea1,"*",I2_1)
        E2=paste0(E2_1,"*",lna)
        E2_1k1=paste0(deparse(DD(parse(text=E2_1),"k",1)),collapse="")
        E2z_1=paste0(E2_1k1,"-",E2)
        k<-1
        c1<-1
        Parmen[which(Parnam=="k")]=k
        Parmen[which(Parnam=="c1")]=c1
        work=work-3*subsi1c(E2z_1,nf,Parnam,Parmen,WW6)

        work1=as.expression(ccc[[rj1,rj2,np]])
        ccc[[rj1,rj2,np]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
        work1=as.expression(ccc[[rj1,np,rj2]])
        ccc[[rj1,np,rj2]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
        work1=as.expression(ccc[[np,rj1,rj2]])
        ccc[[np,rj1,rj2]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
        for (rj3 in 1:(np-1)){
          u123=uuu[rj1,rj2,rj3]
          u12=uu[rj1,rj2]
          u13=uu[rj1,rj3]
          u23=uu[rj2,rj3]
          u1=u[rj1]
          u2=u[rj2]
          u3=u[rj3]
          work=expression(((exp(-g)+d)/((exp(-g)+n)*(exp(-g)+n+1)))*(u12*u3+u13*u2+u23*u1+3*u123))
          work=dpg(as.expression(eval(substitute(substitute(e, list(u12=u12,u3=u3,u13=u13,u2=u2,u23=u23,u1=u1,u123=u123)), list(e = work[[1]])))))
          work1=as.expression(ccc[[rj1,rj2,rj3]])
          ccc[[rj1,rj2,rj3]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
        }
      }
    }
    E2=paste0(nud,"*",I2)
    E2k1k2=gsub("\n","",paste0(deparse(DD(DD(parse(text=E2),"k1",1),"k2",2)),collapse=""))
    k1<-1
    pl<-3
    c1<-1
    k2<-1
    c2<-1
    Parmen[which(Parnam=="k1")]=k1
    Parmen[which(Parnam=="pl")]=pl
    Parmen[which(Parnam=="c1")]=c1
    Parmen[which(Parnam=="k2")]=k2
    Parmen[which(Parnam=="c2")]=c2
    E2_3=paste0(nud,"*",lna,"*",I2)
    E2_3k2=gsub("\n","",paste0(deparse(DD(parse(text=E2_3),"k2",2)),collapse=""))

    E2_4=paste0(nud,"*",lnaa,"^2*",I2)
    E2_4k1=gsub("\n","",paste0(deparse(DD(parse(text=E2_4),"k1",1)),collapse=""))

    E2_5=paste0(nud,"*",lnaa,"*",I2)
    E2_5k1k2=gsub("\n","",paste0(deparse(DD(DD(parse(text=E2_5),"k1",1),"k2",1)),collapse=""))

    E2_6=paste0(nud,"*",lna,"*",lnaa,"*",I2)
    E2_6k2=gsub("\n","",paste0(deparse(DD(parse(text=E2_6),"k2",1)),collapse=""))

    E2_3=paste0(nud,"*",lna,"*",lnaa,"^2*",I2)
    E2z=paste0(E2k1k2,"-",E2_3,"-",E2_3k2,"+",E2_4k1,"-2*",E2_5k1k2,"+2*",E2_6k2)
    work=subsi2c(E2z,nf,Parnam,Parmen,WW4,WW5)-subsi2c_i(E2z,nf,Parnam,Parmen,WW6)
    E2_3=paste0(nud,"*",I2_1)
    E2_0=paste0(E2_3,"*",lna,"^3")
    E2_1=paste0(E2_3,"*",lna,"^2")
    E2_2=paste0(E2_3,"*",lna)
    E2_1k1=paste0(deparse(DD(parse(text=E2_1),"k",1)),collapse="")
    E2_2k2=paste0(deparse(DD(parse(text=E2_2),"k",2)),collapse="")
    E2_3k3=paste0(deparse(DD(parse(text=E2_3),"k",3)),collapse="")
    E2z_1=paste0(E2_3k3,"-3*",E2_2k2,"+3*",E2_1k1,"-",E2_0)
    k<-2
    Parmen[which(Parnam=="k")]=k
    work=work+subsi1c(E2z_1,nf,Parnam,Parmen,WW6)
    work1=as.expression(ccc[[np,np,np]])
    ccc[[np,np,np]]=substitute(3*work+work1,list(work=work[[1]],work1=work1[[1]]))
####################################F3.6
    for (rj1 in 1:(np-1)){
      u1=sum(u[rj1])
      E3_1=III3
      E3_1[,1]=paste0("2*mud*",E3_1[,1])
      E3_1_3=Diff(Diff(E3_1,"k2",1),"k3",1)
      E3_0_3=change(E3_1,0,1,1)
      E3_0_1=change(E3_1,0,0,1)
      E3_k3=Diff(E3_0_1,"k2",1)
      E3_0_2=change(E3_1,0,1,0)
      E3_k2=Diff(E3_0_2,"k3",1)
      k1<-1
      k2<-1
      k3<-1
      pl<-2
      c1<-0
      c2<-1
      c3<-1
      ND1<-k1+1
      ND2<-k2+1
      ND3<-k3+1
      Deg0<-3
      Parmen[which(Parnam=="k1")]=k1
      Parmen[which(Parnam=="k2")]=k2
      Parmen[which(Parnam=="k3")]=k3
      Parmen[which(Parnam=="pl")]=pl
      Parmen[which(Parnam=="c1")]=c1
      Parmen[which(Parnam=="c2")]=c2
      Parmen[which(Parnam=="c3")]=c3
      Parmen[which(Parnam=="ND1")]=ND1
      Parmen[which(Parnam=="ND2")]=ND2
      Parmen[which(Parnam=="ND3")]=ND3
      Parmen[which(Parnam=="Deg0")]=Deg0
E3con=convert(E3_1_3,rj1,0,Parnam,Parmen)-convert(E3_k2,rj1,0,Parnam,Parmen)-convert(E3_k3,rj1,0,Parnam,Parmen)+convert(E3_0_3,rj1,0,Parnam,Parmen)
      E3_1=paste0("2*x_",as.character(rj1),"*",nud,"*",I3)
      E3_1_3=gsub("\n","",paste0(deparse(DD(DD(parse(text=E3_1),"k2",1),"k3",1)),collapse=""))
      E3_0_3=paste0(E3_1,"*",lnaa,"*",lnaaa)
      E3_0_1=paste0(E3_1,"*",lnaaa)
      E3_k3=paste0(deparse(DD(parse(text=E3_0_1),"k2",1)),collapse="")
      E3_0_2=paste0(E3_1,"*",lnaa)
      E3_k2=paste0(deparse(DD(parse(text=E3_0_2),"k3",1)),collapse="")
      k1<-1
      k2<-1
      k3<-1
      pl<-2
      c1<-0
      c2<-1
      c3<-1
      Parmen[which(Parnam=="k1")]=k1
      Parmen[which(Parnam=="k2")]=k2
      Parmen[which(Parnam=="k3")]=k3
      Parmen[which(Parnam=="pl")]=pl
      Parmen[which(Parnam=="c1")]=c1
      Parmen[which(Parnam=="c2")]=c2
      Parmen[which(Parnam=="c3")]=c3
      E3=paste0(E3_1_3,"-",E3_k2,"-",E3_k3,"+",E3_0_3)
      work=E3con-subsi3c_i12(E3,nf,Parnam,Parmen,WW4,WW5)-subsi3c_i13(E3,nf,Parnam,Parmen,WW4,WW5)-subsi3c_i23(E3,nf,Parnam,Parmen,WW4,WW5)+2*subsi3c_i123(E3,nf,Parnam,Parmen,WW6)
      E3_1=paste0("2*x_",as.character(rj1),"*",nud,"*",I3_2)
      E3_1_3=gsub("\n","",paste0(deparse(DD(DD(parse(text=E3_1),"k1",1),"k2",1)),collapse=""))
      E3_0_3=paste0(E3_1,"*",lna,"*",lnaa)
      E3_0_1=paste0(E3_1,"*",lnaa)
      E3_k3=paste0(deparse(DD(parse(text=E3_0_1),"k1",1)),collapse="")
      E3_0_2=paste0(E3_1,"*",lna)
      E3_k2=paste0(deparse(DD(parse(text=E3_0_2),"k2",1)),collapse="")
      k1<-2
      k2<-1
      c1<-1
      c2<-1
      Parmen[which(Parnam=="k1")]=k1
      Parmen[which(Parnam=="k2")]=k2
      Parmen[which(Parnam=="c1")]=c1
      Parmen[which(Parnam=="c2")]=c2
      E3=paste0(E3_1_3,"-",E3_k2,"-",E3_k3,"+",E3_0_3)
      work=work+2*subsi2c(E3,nf,Parnam,Parmen,WW4,WW5)-2*subsi2c_i(E3,nf,Parnam,Parmen,WW6)
      E3_1=paste0("2*x_",as.character(rj1),"*",nud,"*",ea1,"*",I3_2)
      E3_1_3=gsub("\n","",paste0(deparse(DD(DD(parse(text=E3_1),"k1",1),"k2",1)),collapse=""))
      E3_0_3=paste0(E3_1,"*",lna,"*",lnaa)
      E3_0_1=paste0(E3_1,"*",lnaa)
      E3_k3=paste0(deparse(DD(parse(text=E3_0_1),"k1",1)),collapse="")
      E3_0_2=paste0(E3_1,"*",lna)
      E3_k2=paste0(deparse(DD(parse(text=E3_0_2),"k2",1)),collapse="")
      k1<-1
      k2<-1
      c1<-1
      c2<-1
      Parmen[which(Parnam=="k1")]=k1
      Parmen[which(Parnam=="k2")]=k2
      Parmen[which(Parnam=="c1")]=c1
      Parmen[which(Parnam=="c2")]=c2
      E3=paste0(E3_1_3,"-",E3_k2,"-",E3_k3,"+",E3_0_3)
      work=work-2*subsi2c(E3,nf,Parnam,Parmen,WW4,WW5)+2*subsi2c_i(E3,nf,Parnam,Parmen,WW6)
      k1<-1
      k2<-2
      c1<-0
      c2<-1
      Parmen[which(Parnam=="k1")]=k1
      Parmen[which(Parnam=="k2")]=k2
      Parmen[which(Parnam=="c1")]=c1
      Parmen[which(Parnam=="c2")]=c2
      E2_12_1=paste0("2*x_",as.character(rj1),"*",nud,"*",I3_2)
      E2_12_0=paste0(E2_12_1,"*",lnaa)
      E2_12_2=paste0(E2_12_1,"*",lnaa,"^2")
      E2_12k2=paste0(deparse(DD(parse(text=E2_12_1),"k2",2)),collapse="")
      E2_13k2=paste0(deparse(DD(parse(text=E2_12_0),"k2",1)),collapse="")
      E2_1=paste0(E2_12k2,"-2*",E2_13k2,"+",E2_12_2)
      work=work+subsi2c(E2_1,nf,Parnam,Parmen,WW4,WW5)-subsi2c_i(E2_1,nf,Parnam,Parmen,WW6)
      k<-3
      c1<-1
      Parmen[which(Parnam=="k")]=k
      Parmen[which(Parnam=="c1")]=c1

      E1_3_1=paste0("2*x_",as.character(rj1),"*",nud,"*",I3_1)
      E1_3_0=paste0(E1_3_1,"*",lna)
      E1_3_2=paste0(E1_3_1,"*",lna,"^2")
      E1_3k=paste0(deparse(DD(parse(text=E1_3_1),"k",2)),collapse="")
      E1_3k1=paste0(deparse(DD(parse(text=E1_3_0),"k",1)),collapse="")
      E1_3=paste0(E1_3k,"-2*",E1_3k1,"+",E1_3_2)
      work=work+subsi1c(E1_3,nf,Parnam,Parmen,WW6)
      k<-2
      c1<-1
      Parmen[which(Parnam=="k")]=k
      Parmen[which(Parnam=="c1")]=c1
      E1_3_1=paste0("2*x_",as.character(rj1),"*",nud,"*",ea1,"*",I3_1)
      E1_3_0=paste0(E1_3_1,"*",lna)
      E1_3_2=paste0(E1_3_1,"*",lna,"^2")
      E1_3k=paste0(deparse(DD(parse(text=E1_3_1),"k",2)),collapse="")
      E1_3k1=paste0(deparse(DD(parse(text=E1_3_0),"k",1)),collapse="")
      E1_3=paste0(E1_3k,"-2*",E1_3k1,"+",E1_3_2)
      work=work-subsi1c(E1_3,nf,Parnam,Parmen,WW6)
      work1=as.expression(ccc[[rj1,np,np]])
      ccc[[rj1,np,np]]=substitute(-work+work1,list(work=work[[1]],work1=work1[[1]]))
      work1=as.expression(ccc[[np,rj1,np]])
      ccc[[np,rj1,np]]=substitute(-work+work1,list(work=work[[1]],work1=work1[[1]]))
      work1=as.expression(ccc[[np,np,rj1]])
      ccc[[np,np,rj1]]=substitute(-work+work1,list(work=work[[1]],work1=work1[[1]]))
      for (rj2 in 1:(np-1)){
        u2=u[rj2]
        u12=uu[rj1,rj2]
        ptm <- proc.time()
        E3_1=III3
        E3_1[,1]=paste0("2*mud*",E3_1[,1])
        E3_0=change(E3_1,0,0,1)
        E3_k3=Diff(E3_1,"k3",1)
        k1<-1
        k2<-1
        k3<-1
        pl<-1
        c1<-0
        c2<-0
        c3<-1
        ND1<-k1+1
        ND2<-k2+1
        ND3<-k3+1
        Deg0<-3
        Parmen[which(Parnam=="k1")]=k1
        Parmen[which(Parnam=="k2")]=k2
        Parmen[which(Parnam=="k3")]=k3
        Parmen[which(Parnam=="pl")]=pl
        Parmen[which(Parnam=="c1")]=c1
        Parmen[which(Parnam=="c2")]=c2
        Parmen[which(Parnam=="c3")]=c3
        Parmen[which(Parnam=="ND1")]=ND1
        Parmen[which(Parnam=="ND2")]=ND2
        Parmen[which(Parnam=="ND3")]=ND3
        Parmen[which(Parnam=="Deg0")]=Deg0
        E3con=convert(E3_k3,rj1,rj2,Parnam,Parmen)-convert(E3_0,rj1,rj2,Parnam,Parmen)
        E3_1=paste0("2*x_",as.character(rj1),"*x__",as.character(rj2),"*",nud,"*",I3)
        E3_0=paste0(E3_1,"*",lnaaa)
        E3_k3=paste0(deparse(DD(parse(text=E3_1),"k3",1)),collapse="")
        k1<-1
        k2<-1
        k3<-1
        pl<-1
        c1<-0
        c2<-0
        c3<-1
        Parmen[which(Parnam=="k1")]=k1
        Parmen[which(Parnam=="k2")]=k2
        Parmen[which(Parnam=="k3")]=k3
        Parmen[which(Parnam=="pl")]=pl
        Parmen[which(Parnam=="c1")]=c1
        Parmen[which(Parnam=="c2")]=c2
        Parmen[which(Parnam=="c3")]=c3

        E3=paste0(E3_k3,"-",E3_0)
        work=E3con-subsi3c_i12(E3,nf,Parnam,Parmen,WW4,WW5)-subsi3c_i13(E3,nf,Parnam,Parmen,WW4,WW5)-subsi3c_i23(E3,nf,Parnam,Parmen,WW4,WW5)+2*subsi3c_i123(E3,nf,Parnam,Parmen,WW6)
        k1<-2
        k2<-1
        c1<-0
        c2<-1
        Parmen[which(Parnam=="k1")]=k1
        Parmen[which(Parnam=="k2")]=k2
        Parmen[which(Parnam=="c1")]=c1
        Parmen[which(Parnam=="c2")]=c2
        E2_12_1=paste0("2*x_",as.character(rj1),"*","x_",as.character(rj2),"*",nud,"*",I3_2)
        E2_12_0=paste0(E2_12_1,"*",lnaa)
        E2_12k2=paste0(deparse(DD(parse(text=E2_12_1),"k2",1)),collapse="")
        E2_1=paste0(E2_12k2,"-",E2_12_0)
        work=work+subsi2c(E2_1,nf,Parnam,Parmen,WW4,WW5)-subsi2c_i(E2_1,nf,Parnam,Parmen,WW6)
        k1<-2
        k2<-1
        c1<-1
        c2<-0
        Parmen[which(Parnam=="k1")]=k1
        Parmen[which(Parnam=="k2")]=k2
        Parmen[which(Parnam=="c1")]=c1
        Parmen[which(Parnam=="c2")]=c2
        E2_13_1=paste0("2*x_",as.character(rj1),"*x__",as.character(rj2),"*",nud,"*",I3_2)
        E2_13_0=paste0(E2_13_1,"*",lna)
        E2_13k1=paste0(deparse(DD(parse(text=E2_13_1),"k1",1)),collapse="")
        E2_2=paste0(E2_13k1,"-",E2_13_0)
        work=work+subsi2c(E2_2,nf,Parnam,Parmen,WW4,WW5)-subsi2c_i(E2_2,nf,Parnam,Parmen,WW6)
        k1<-1
        k2<-1
        c1<-1
        c2<-0
        Parmen[which(Parnam=="k1")]=k1
        Parmen[which(Parnam=="k2")]=k2
        Parmen[which(Parnam=="c1")]=c1
        Parmen[which(Parnam=="c2")]=c2
        E2_13_1=paste0("2*x_",as.character(rj1),"*x__",as.character(rj2),"*",nud,"*",ea1,"*",I3_2)
        E2_13_0=paste0(E2_13_1,"*",lna)
        E2_13k1=paste0(deparse(DD(parse(text=E2_13_1),"k1",1)),collapse="")
        E2_2=paste0(E2_13k1,"-",E2_13_0)
        work=work-subsi2c(E2_2,nf,Parnam,Parmen,WW4,WW5)+subsi2c_i(E2_2,nf,Parnam,Parmen,WW6)

        k1<-1
        k2<-2
        c1<-0
        c2<-1
        Parmen[which(Parnam=="k1")]=k1
        Parmen[which(Parnam=="k2")]=k2
        Parmen[which(Parnam=="c1")]=c1
        Parmen[which(Parnam=="c2")]=c2

        E2_23_1=paste0("2*x_",as.character(rj1),"*x__",as.character(rj2),"*",nud,"*",I3_2)
        E2_23_0=paste0(E2_23_1,"*",lnaa)
        E2_23k2=paste0(deparse(DD(parse(text=E2_23_1),"k2",1)),collapse="")
        E2_3=paste0(E2_23k2,"-",E2_23_0)
        work=work+subsi2c(E2_3,nf,Parnam,Parmen,WW4,WW5)-subsi2c_i(E2_3,nf,Parnam,Parmen,WW6)
        k1<-1
        k2<-1
        c1<-0
        c2<-1
        Parmen[which(Parnam=="k1")]=k1
        Parmen[which(Parnam=="k2")]=k2
        Parmen[which(Parnam=="c1")]=c1
        Parmen[which(Parnam=="c2")]=c2

        E2_23_1=paste0("2*x_",as.character(rj1),"*x__",as.character(rj2),"*",nud,"*",eaa1,"*",I3_2)
        E2_23_0=paste0(E2_23_1,"*",lnaa)
        E2_23k2=paste0(deparse(DD(parse(text=E2_23_1),"k2",1)),collapse="")
        E2_3=paste0(E2_23k2,"-",E2_23_0)
        work=work-subsi2c(E2_3,nf,Parnam,Parmen,WW4,WW5)+subsi2c_i(E2_3,nf,Parnam,Parmen,WW6)


        k<-3
        c1<-1
        Parmen[which(Parnam=="k")]=k
        Parmen[which(Parnam=="c1")]=c1

        E1_3_1=paste0("2*x_",as.character(rj1),"*x_",as.character(rj2),"*",nud,"*",I3_1)
        E1_3_0=paste0(E1_3_1,"*",lna)
        E1_3k=paste0(deparse(DD(parse(text=E1_3_1),"k",1)),collapse="")
        E1_3=paste0(E1_3k,"-",E1_3_0)
        work=work+subsi1c(E1_3,nf,Parnam,Parmen,WW6)
        k<-2
        c1<-1
        Parmen[which(Parnam=="k")]=k
        Parmen[which(Parnam=="c1")]=c1
        E1_3_1=paste0("2*x_",as.character(rj1),"*x_",as.character(rj2),"*",nud,"*",ea1,"*",I3_1)
        E1_3_0=paste0(E1_3_1,"*",lna)
        E1_3k=paste0(deparse(DD(parse(text=E1_3_1),"k",1)),collapse="")
        E1_3=paste0(E1_3k,"-",E1_3_0)
        work=work-subsi1c(E1_3,nf,Parnam,Parmen,WW6)
        work1=as.expression(ccc[[rj1,rj2,np]])
        ccc[[rj1,rj2,np]]=substitute(work1-work,list(work1=work1[[1]],work=work[[1]]))
        work1=as.expression(ccc[[rj1,np,rj2]])
        ccc[[rj1,np,rj2]]=substitute(work1-work,list(work1=work1[[1]],work=work[[1]]))
        work1=as.expression(ccc[[np,rj1,rj2]])
        ccc[[np,rj1,rj2]]=substitute(work1-work,list(work1=work1[[1]],work=work[[1]]))
        for (rj3 in 1:(np-1)){
          u123=uuu[rj1,rj2,rj3]
          u12=uu[rj1,rj2]
          u13=uu[rj1,rj3]
          u23=uu[rj2,rj3]
          u1=u[rj1]
          u2=u[rj2]
          u3=u[rj3]
          work=expression((-2*(exp(-g)+d)/((exp(-g)+n)*(exp(-g)+n+1)*(exp(-g)+n+2)))*(u1*u2*u3+u12*u3+u13*u2+u23*u1+2*u123))
          work=dpg(as.expression(eval(substitute(substitute(e, list(u12=u12,u13=u13,u23=u23,u1=u1,u2=u2,u3=u3,u123=u123)), list(e = work[[1]])))))

          work1=as.expression(ccc[[rj1,rj2,rj3]])
          ccc[[rj1,rj2,rj3]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
        }
      }
    }
    E3_1=III3
    E3_1[,1]=paste0("2*mud*",E3_1[,1])
    E3_1_3=Diff(Diff(Diff(E3_1,"k1",1),"k2",1),"k3",1)
    E3_0_0=change(E3_1,1,1,1)
    E3_0_1=change(E3_1,1,1,0)
    E3_0_2=change(E3_1,1,0,1)
    E3_0_3=change(E3_1,0,1,1)
    E3_0_1k3=Diff(E3_0_1,"k3",1)

    E3_0_2k2=Diff(E3_0_2,"k2",1)
    E3_0_3k1=Diff(E3_0_3,"k1",1)
    E3_2_1=change(E3_1,1,0,0)
    E3_2_2=change(E3_1,0,1,0)
    E3_2_3=change(E3_1,0,0,1)
    E3_2_1k23=Diff(Diff(E3_2_1,"k2",1),"k3",1)
    E3_2_2k13=Diff(Diff(E3_2_2,"k1",1),"k3",1)
    E3_2_3k12=Diff(Diff(E3_2_3,"k1",1),"k2",1)
    k1<-1
    k2<-1
    k3<-1
    pl<-3
    ND1<-k1+1
    ND2<-k2+1
    ND3<-k3+1
    Deg0<-3
    Parmen[which(Parnam=="k1")]=k1
    Parmen[which(Parnam=="ND1")]=ND1
    Parmen[which(Parnam=="k2")]=k2
    Parmen[which(Parnam=="ND2")]=ND2
    Parmen[which(Parnam=="k3")]=k3
    Parmen[which(Parnam=="ND3")]=ND3
    Parmen[which(Parnam=="pl")]=pl
    Parmen[which(Parnam=="Deg0")]=Deg0
    E3con=convert(E3_1_3,0,0,Parnam,Parmen)-convert(E3_2_1k23,0,0,Parnam,Parmen)-convert(E3_2_2k13,0,0,Parnam,Parmen)-convert(E3_2_3k12,0,0,Parnam,Parmen)+convert(E3_0_1k3,0,0,Parnam,Parmen)+convert(E3_0_2k2,0,0,Parnam,Parmen)+convert(E3_0_3k1,0,0,Parnam,Parmen)-convert(E3_0_0,0,0,Parnam,Parmen)
    E3_1=paste0("2*",nud,"*",I3)
    E3_1_3=gsub("\n","",paste0(deparse(DD(DD(DD(parse(text=E3_1),"k1",1),"k2",1),"k3",1)),collapse=""))
    E3_0_0=paste0(E3_1,"*",lna,"*",lnaa,"*",lnaaa)
    E3_0_1=paste0(E3_1,"*",lna,"*",lnaa)
    E3_0_2=paste0(E3_1,"*",lna,"*",lnaaa)
    E3_0_3=paste0(E3_1,"*",lnaa,"*",lnaaa)
    E3_0_1k3=paste0(deparse(DD(parse(text=E3_0_1),"k3",1)),collapse="")
    E3_0_2k2=paste0(deparse(DD(parse(text=E3_0_2),"k2",1)),collapse="")
    E3_0_3k1=paste0(deparse(DD(parse(text=E3_0_3),"k1",1)),collapse="")
    E3_2_1=paste0(E3_1,"*",lna)
    E3_2_2=paste0(E3_1,"*",lnaa)
    E3_2_3=paste0(E3_1,"*",lnaaa)
    E3_2_1k23=paste0(deparse(DD(DD(parse(text=E3_2_1),"k2",1),"k3",1)),collapse="")
    E3_2_2k13=paste0(deparse(DD(DD(parse(text=E3_2_2),"k1",1),"k3",1)),collapse="")
    E3_2_3k12=paste0(deparse(DD(DD(parse(text=E3_2_3),"k1",1),"k2",1)),collapse="")
    k1<-1
    k2<-1
    k3<-1
    pl<-3
    c1<-1
    c2<-1
    c3<-1
    Parmen[which(Parnam=="k1")]=k1
    Parmen[which(Parnam=="c1")]=c1
    Parmen[which(Parnam=="k2")]=k2
    Parmen[which(Parnam=="c2")]=c2
    Parmen[which(Parnam=="k3")]=k3
    Parmen[which(Parnam=="c3")]=c3
    Parmen[which(Parnam=="pl")]=pl
    E3=paste0(E3_1_3,"-",E3_2_1k23,"-",E3_2_2k13,"-",E3_2_3k12,"+",E3_0_1k3,"+",E3_0_2k2,"+",E3_0_3k1,"-",E3_0_0)
    work=E3con-subsi3c_i12(E3,nf,Parnam,Parmen,WW4,WW5)-subsi3c_i13(E3,nf,Parnam,Parmen,WW4,WW5)-subsi3c_i23(E3,nf,Parnam,Parmen,WW4,WW5)+2*subsi3c_i123(E3,nf,Parnam,Parmen,WW6)
    E3_1=paste0("2","*",nud,"*",I3_2)
    E3_1_3k12=gsub("\n","",paste0(deparse(DD(DD(parse(text=E3_1),"k1",2),"k2",1)),collapse=""))
    E3_0_3=paste0("2*",E3_1,"*",lna,"*",lnaa)
    E3_k=paste0(deparse(DD(parse(text=E3_0_3),"k1",1)),collapse="")
    E3_0_1=paste0(E3_1,"*",lnaa)
    E3_k3=paste0(deparse(DD(parse(text=E3_0_1),"k1",2)),collapse="")
    E3_0_2=paste0(E3_1,"*",lna,"^2*",lnaa)
    E3_0_5=paste0(E3_1,"*",lna,"^2")
    E3_5k2=paste0(deparse(DD(parse(text=E3_0_5),"k2",1)),collapse="")
    E3_0_4=paste0("2*",E3_1,"*",lna)
    E3_1_3kk12=gsub("\n","",paste0(deparse(DD(DD(parse(text=E3_0_4),"k1",1),"k2",1)),collapse=""))
    k1<-2
    k2<-1
    Parmen[which(Parnam=="k1")]=k1
    Parmen[which(Parnam=="k2")]=k2
    E3=paste0(E3_1_3k12,"-",E3_k3,"-",E3_1_3kk12,"+",E3_k,"+",E3_5k2,"-",E3_0_2)
    work=work+3*subsi2c(E3,nf,Parnam,Parmen,WW4,WW5)-3*subsi2c_i(E3,nf,Parnam,Parmen,WW6)
    k<-3
    Parmen[which(Parnam=="k")]=k
    E1_3_1=paste0("2","*",nud,"*",I3_1)
    E1_3_0=paste0(E1_3_1,"*",lna)
    E1_3_2=paste0(E1_3_1,"*",lna,"^2")
    E1_3_3=paste0(E1_3_1,"*",lna,"^3")
    E1_3k=paste0(deparse(DD(parse(text=E1_3_1),"k",3)),collapse="")
    E1_3k1=paste0(deparse(DD(parse(text=E1_3_0),"k",2)),collapse="")
    E1_3k2=paste0(deparse(DD(parse(text=E1_3_2),"k",1)),collapse="")
    E1_3=paste0(E1_3k,"-3*",E1_3k1,"+3*",E1_3k2,"-",E1_3_3)
    work=work+subsi1c(E1_3,nf,Parnam,Parmen,WW6)
    work0=as.expression(ccc[[np,np,np]])
    ccc[[np,np,np]]=substitute(work0-work,list(work0=work0[[1]],work=work[[1]]))

    Rccc=array(NA,c(np,np,np))
    for (i1 in 1:np){
      for (i2 in 1:np){
        for (i3 in 1:np){
          Rccc[[i1,i2,i3]]=eval(ccc[[i1,i2,i3]])
        }
      }
    }
    ###########################FFinal
    Rccc[nf+2,nf+2,nf+2]=Rccc[nf+2,nf+2,nf+2]+exp(-g)*sum(1/(exp(-g)+(0:(n-1))))
    RCCC=RCCC+Rccc
    RCC=RCC+Rcc
    RCC_C=RCC_C+Rcc_c
  }
  if (!(any(is.na(RCC)) | any(is.na(RCCC)))){
    RCC_=ginv(RCC)
    kra=RCC_C-0.5*RCCC
    RCC_=ginv(RCC)
    Bias={}
    np=nf+3
    for (ir in 1:np){
      work=0
      for (ir1 in 1:np){
        for (ir2 in 1:np){
          for (ir3 in 1:np){
            work=work+kra[ir1,ir2,ir3]*RCC_[ir,ir1]*RCC_[ir2,ir3]
          }
        }
      }
      Bias=c(Bias,work)
    }
  } else {
    Bias=rep(NA,np)
  }
  return(Bias)
}
