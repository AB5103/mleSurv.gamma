#' @title BiasCorrectionW
#'
#' @description   Returns the bias correction for ML estimates.
#'
#' @param D The data frame. The first 'nf' columns are the covariates under study. Other columns must include fields
#' 'time' for time-to-failure, 'Cens' for censoring and 'cluster' for cluster identification. The number of rows is equal to the number of subjects.
#' @param nf The number of covariates under study.
#' @param ncl The number of clusters.
#' @param para  The vector of the ML estimates of unknown parameters. This vector has a form \deqn{y = (ln(a), \beta , ln(\sigma ^2), b),}
#' where \eqn{a} and \eqn{b} are the slope and the shape parameters, \eqn{\beta } is the vector of the Cox-regression parameters,
#' and \eqn{\sigma ^2} is the variance of the frailty.
#'
#' @details This function calculates the bias correction for the proportianal hazards model with gamma frailty
#' in accordance with the Cox-Snell method. The baseline cumulative hazard function is \deqn{H_{0}(t;a,b)=at^b}.
#'
#' @return The bias correction for vector \eqn{y}. If the calculation of the bias correction is
#' not possible returns NAs.
#'
#' @import MASS
#'
#' @examples
 #' \dontrun{
#' BiasCorrectionW(D,nf,ncl,para)
#'
#'}
#' @export
#'
BiasCorrectionW=function(D,nf,ncl,para){
  parnam=c("a",paste0("beta_",as.character(1:nf)),"g","b")
  nr=nrow(D)
  np=nf+3
  a.scale=rep(1,nr)
  g.scale=rep(1,nr)
  uR_beta=D[,1:nf]
  UR=cbind(a.scale,uR_beta,g.scale)
  Cens=D$event
  a=para[1]
  b=para[nf+3]
  g=para[nf+2]
  par=c(a,para[2:(1+nf)],g,b)
  parnam=c("a",paste0("beta_",as.character(1:nf)),"g","b")
  fi="a"
  if (nf>0){
    for (i in 1:nf){
      assign(paste0("beta_",as.character(i)),par[1+i])
      fi=c(fi,paste0("beta_",as.character(i)))
    }
  }
  fi=c(fi,"g")
  list=sort(unique(D$cluster))
  RCC=array(0,c(np,np))
  RCCC=array(0,c(np,np,np))
  RCC_C=array(0,c(np,np,np))

  for (lj in 1:ncl){
    ind=which(D$cluster==list[lj])
    n=length(ind)
    d=sum(Cens[ind])
    li=length(ind)
    if (li==1)  {U=matrix(as.numeric(UR[ind,]),length(ind),ncol(UR))} else {
      U=UR[ind,]}
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

    ufiu=c(rep(NA,nf+2))
    uufiu=array(expression(nul),c((nf+2),(nf+2)))
    ufiufiu=c(rep(NA,nf+2))
    fiufiu=NA
    work2={}
    work1={}
    work4={}
    j=0
    for (j1 in 1:(nf+2)){
      work={}
      work2={}
      for (j2 in 1:(nf+2)){
        j=j+1
        if (j2<(nf+2)) work=paste0(work,as.character(uu[j1,j2]),"*",fi[j2],"+")
        if (j2==(nf+2)) work=paste0(work,as.character(uu[j1,j2]),"*",fi[j2])
        if (!(j1==(nf+2) & (j2==(nf+2)))) work1=paste0(work1,as.character(uu[j1,j2]),"*",fi[j1],"*",fi[j2],"+")
        if (j1==(nf+2) & j2==(nf+2)) work1=paste0(work1,as.character(uu[j1,j2]),"*",fi[j1],"*",fi[j2])
        work0={}
        for (j3 in 1:(nf+2)){
          if (j3<(nf+2)) work0=paste0(work0,as.character(uuu[j1,j2,j3]),"*",fi[j3],"+")
          if (j3==(nf+2)) work0=paste0(work0,as.character(uuu[j1,j2,j3]),"*",fi[j3])
          if (j3==(nf+2) & j2==(nf+2)) work2=paste0(work2,as.character(uuu[j1,j2,j3]),"*",fi[j2],"*",fi[j3])
          if (!(j3==(nf+2) & j2==(nf+2))) work2=paste0(work2,as.character(uuu[j1,j2,j3]),"*",fi[j2],"*",fi[j3],"+")
          if (!(j3==(nf+2) & j2==(nf+2) & j1==(nf+2))) work4=paste0(work4,as.character(uuu[j1,j2,j3]),"*",fi[j1],"*",fi[j2],"*",fi[j3],"+")
          if (j3==(nf+2) & j2==(nf+2) & j1==(nf+2)) work4=paste0(work4,as.character(uuu[j1,j2,j3]),"*",fi[j1],"*",fi[j2],"*",fi[j3])
        }
        uufiu[[j1,j2]]=parse(text=work0)
      }
      ufiu[j1]=parse(text=work)
      ufiufiu[j1]=parse(text=work2)
    }
    fiufiu=parse(text=work1)
    fiufiufiu=parse(text=work4)

    fiufi={}
    for (i in 1:(nf+2)){
      if (i<(nf+2)) fiufi=paste0(fiufi,as.character(u[i]),"*",fi[i],"+")
      if (i==(nf+2)) fiufi=paste0(fiufi,as.character(u[i]),"*",fi[i])
    }
    fiufi=parse(text=fiufi)
    CoxF=expression(a+Cox+g)

    Part2_2=expression(-d/b^2)
    Part2_3=expression(2*d/b^3)

    Part3=expression(lgamma(exp(-g)+d)-lgamma(exp(-g)))
    Part3_2=as.expression(DD(Part3,"g",2))
    Part3_3=as.expression(DD(Part3,"g",3))
    ##########################################
    #Second order derivative
    ##########################################
    cc=array(expression(nul),c(np,np))
    cc[[np,np]]=substitute(Part2_2,list(Part2_2=Part2_2[[1]]))
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
    work=expression((exp(-g)/((exp(-g)+n)*b))*(n*(digamma(2)-digamma(exp(-g)))-fiufi))
    work=dpg(as.expression(eval(substitute(substitute(e, list(fiufi=fiufi)), list(e = work[[1]])))))
    cc[[(nf+2),(nf+3)]]=substitute(work,list(work=work[[1]]))
    cc[[(nf+3),(nf+2)]]=cc[[(nf+2),(nf+3)]]
    ###########################F2.2
    for (j1 in 1:(np-1)){
      u1=u[j1]
      ufiu1=ufiu[j1]
      work=expression(-(exp(-g)+d)*((digamma(2)-digamma(exp(-g)))*u1-ufiu1)/((exp(-g)+n)*b))
      work=dpg(as.expression(eval(substitute(substitute(e, list(u1=u1,ufiu1=ufiu1)), list(e = work[[1]])))))
      work1=as.expression(cc[[j1,np]])
      cc[[j1,np]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
      work1=as.expression(cc[[np,j1]])
      cc[[np,j1]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
      for (j2 in 1:(np-1)){
        u12=uu[[j1,j2]]
        work=expression(-(exp(-g)+d)*u12/(exp(-g)+n))
        work=dpg(as.expression(eval(substitute(substitute(e, list(u12=u12)), list(e = work[[1]])))))
        work1=as.expression(cc[[j1,j2]])
        cc[[j1,j2]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
      }
    }
    work0=expression(-(exp(-g)+d)*n*((digamma(2)-digamma(exp(-g)))^2+trigamma(2)+trigamma(exp(-g)))/((exp(-g)+n)*b^2))
    work1=as.expression(cc[[np,np]])
    work=expression(2*(exp(-g)+d)*((digamma(2)-digamma(exp(-g)))*fiufi)/((exp(-g)+n)*b^2))
    work=dpg(as.expression(eval(substitute(substitute(e, list(fiufi=fiufi)), list(e = work[[1]])))))
    work2=expression(-(exp(-g)+d)*fiufiu/((exp(-g)+n)*b^2))
    work2=dpg(as.expression(eval(substitute(substitute(e, list(fiufiu=fiufiu)), list(e = work2[[1]])))))


    cc[[np,np]]=substitute(work0+work+work1+work2,list(work0=work0[[1]],work=work[[1]],work1=work1[[1]],work2=work2[[1]]))
    ###########################F2.3
    for (j1 in 1:(np-1)){
      u1=u[j1]
      ufiu1=ufiu[j1]
      work=expression(((exp(-g)+d)/(b*(exp(-g)+n)*(exp(-g)+n+1)))*(((2*digamma(3)-digamma(2)-digamma(exp(-g)))*u1-ufiu1)+n*(digamma(2)-digamma(exp(-g)))*u1-fiufi*u1))
      work=dpg(as.expression(eval(substitute(substitute(e, list(u1=u1,ufiu1=ufiu1,fiufi=fiufi)), list(e = work[[1]])))))
      work1=as.expression(cc[[j1,np]])
      cc[[j1,np]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
      work1=as.expression(cc[[np,j1]])
      cc[[np,j1]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
      for (j2 in 1:(np-1)){
        u2=u[j2]
        u12=uu[j1,j2]
        work=expression(((exp(-g)+d)/((exp(-g)+n)*(exp(-g)+n+1)))*(u12+u1*u2))
        work=dpg(as.expression(eval(substitute(substitute(e, list(u1=u1,u2=u2,u12=u12)), list(e = work[[1]])))))
        work1=as.expression(cc[[j1,j2]])
        cc[[j1,j2]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
      }
    }
    work0=expression((exp(-g)+d)*((n^2-n)*((digamma(2)-digamma(exp(-g)))^2+trigamma(exp(-g)))+2*n*((digamma(3)-digamma(exp(-g)))^2+trigamma(3)+trigamma(exp(-g))))/((exp(-g)+n)*(exp(-g)+n+1)*b^2))
    work1=as.expression(cc[[np,np]])
    work=expression(-2*(exp(-g)+d)*(((n-1)*(digamma(2)-digamma(exp(-g)))+2*(digamma(3)-digamma(exp(-g))))*fiufi)/((exp(-g)+n)*(exp(-g)+n+1)*b^2))
    work=dpg(as.expression(eval(substitute(substitute(e, list(fiufi=fiufi)), list(e = work[[1]])))))
    work2=expression((exp(-g)+d)*(fiufiu+fiufi^2)/((exp(-g)+n)*(exp(-g)+n+1)*b^2))
    work2=dpg(as.expression(eval(substitute(substitute(e, list(fiufiu=fiufiu,fiufi=fiufi)), list(e = work2[[1]])))))
    cc[[np,np]]=substitute(work0+work+work1+work2,list(work0=work0[[1]],work=work[[1]],work1=work1[[1]],work2=work2[[1]]))
    ###########################################
    Rcc=array(NA,c(np,np))
    cc_c=array(expression(1),c(np,np,np))
    Rcc_c=array(NA,c(np,np,np))
    for (i1 in 1:np){
      for (i2 in 1:np){
        Rcc[i1,i2]=eval(cc[[i1,i2]])
        for (i3 in 1:np){
          Rcc_c[[i1,i2,i3]]=eval(DD(cc[[i1,i2]],parnam[i3],1))
          cc_c[[i1,i2,i3]]=DD(cc[[i1,i2]],parnam[i3],1)
        }
      }
    }
    Rcc[nf+2,nf+2]=Rcc[nf+2,nf+2]-exp(-g)*sum(1/(exp(-g)+(0:(n-1))))
    Rcc_c[nf+2,nf+2,nf+2]=Rcc_c[nf+2,nf+2,nf+2]+exp(-g)*sum(1/(exp(-g)+(0:(n-1))))-exp(-2*g)*sum(1/((exp(-g)+(0:(n-1)))^2))
    ##########################################
    #Third order derivative
    ##########################################
    ###########################F3.1
    ccc=array(expression(nul),c(np,np,np))
    ccc[[np,np,np]]=substitute(Part2_3,list(Part2_3=Part2_3[[1]]))
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
    work=expression((-exp(-g)/((exp(-g)+n)*b))*(n*(digamma(2)-digamma(exp(-g)))-fiufi))
    work=dpg(as.expression(eval(substitute(substitute(e, list(fiufi=fiufi)), list(e = work[[1]])))))
    ccc[[(np-1),(np-1),np]]=substitute(work,list(work=work[[1]]))
    ccc[[(np-1),np,(np-1)]]=substitute(work,list(work=work[[1]]))
    ccc[[np,(np-1),(np-1)]]=substitute(work,list(work=work[[1]]))
    ###########################F3.2
    for (j1 in 1:(np-1)){
      u1=u[j1]
      ufiu1=ufiu[j1]
      work=expression(exp(-g)*((digamma(2)-digamma(exp(-g)))*u1-ufiu1)/((exp(-g)+n)*b))
      work=dpg(as.expression(eval(substitute(substitute(e, list(u1=u1,ufiu1=ufiu1)), list(e = work[[1]])))))
      work1=as.expression(ccc[[j1,np,np-1]])
      ccc[[j1,np,np-1]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
      work1=as.expression(ccc[[np,j1,np-1]])
      ccc[[np,j1,np-1]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
      work1=as.expression(ccc[[j1,np-1,np]])
      ccc[[j1,np-1,np]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
      work1=as.expression(ccc[[np,np-1,j1]])
      ccc[[np,np-1,j1]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
      work1=as.expression(ccc[[np-1,j1,np]])
      ccc[[np-1,j1,np]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
      work1=as.expression(ccc[[np-1,np,j1]])
      ccc[[np-1,np,j1]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
      for (j2 in 1:(np-1)){
        u12=uu[j1,j2]
        work=expression(exp(-g)*u12/(exp(-g)+n))
        work=dpg(as.expression(eval(substitute(substitute(e, list(u12=u12)), list(e = work[[1]])))))
        work1=as.expression(ccc[[np-1,j1,j2]])
        ccc[[np-1,j1,j2]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
        work1=as.expression(ccc[[j1,np-1,j2]])
        ccc[[j1,np-1,j2]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
        work1=as.expression(ccc[[j1,j2,np-1]])
        ccc[[j1,j2,np-1]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
      }
    }
    work0=expression(exp(-g)*n*((digamma(2)-digamma(exp(-g)))^2+trigamma(2)+trigamma(exp(-g)))/((exp(-g)+n)*b^2))
    work=expression(-2*exp(-g)*((digamma(2)-digamma(exp(-g)))*fiufi)/((exp(-g)+n)*b^2))
    work=dpg(as.expression(eval(substitute(substitute(e, list(fiufi=fiufi)), list(e = work[[1]])))))
    work2=expression(exp(-g)*fiufiu/((exp(-g)+n)*b^2))
    work2=dpg(as.expression(eval(substitute(substitute(e, list(fiufiu=fiufiu)), list(e = work2[[1]])))))
    work1=as.expression(ccc[[np-1,np,np]])
    ccc[[np-1,np,np]]=substitute(work0+work+work1+work2,list(work0=work0[[1]],work=work[[1]],work1=work1[[1]],work2=work2[[1]]))
    ccc[[np,np-1,np]]=ccc[[np-1,np,np]]
    ccc[[np,np,np-1]]=ccc[[np-1,np,np]]
    ###########################F3.3
    for (j1 in 1:(np-1)){
      u1=u[j1]
      ufiu1=ufiu[j1]
      work=expression((-exp(-g)/((exp(-g)+n)*(exp(-g)+n+1)*b))*(((2*digamma(3)-digamma(2)-digamma(exp(-g)))*u1-ufiu1)+n*(digamma(2)-digamma(exp(-g)))*u1-fiufi*u1))
      work=dpg(as.expression(eval(substitute(substitute(e, list(u1=u1,ufiu1=ufiu1,fiufi=fiufi)), list(e = work[[1]])))))
      work1=as.expression(ccc[[j1,np,np-1]])
      ccc[[j1,np,np-1]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
      work1=as.expression(ccc[[np,j1,np-1]])
      ccc[[np,j1,np-1]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
      work1=as.expression(ccc[[j1,np-1,np]])
      ccc[[j1,np-1,np]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
      work1=as.expression(ccc[[np,np-1,j1]])
      ccc[[np,np-1,j1]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
      work1=as.expression(ccc[[np-1,j1,np]])
      ccc[[np-1,j1,np]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
      work1=as.expression(ccc[[np-1,np,j1]])
      ccc[[np-1,np,j1]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
      for (j2 in 1:(np-1)){
        u2=u[j2]
        u12=uu[j1,j2]
        work=expression((-exp(-g)/((exp(-g)+n)*(exp(-g)+n+1)))*(u12+u1*u2))
        work=dpg(as.expression(eval(substitute(substitute(e, list(u1=u1,u2=u2,u12=u12)), list(e = work[[1]])))))
        work1=as.expression(ccc[[j1,j2,np-1]])
        ccc[[j1,j2,np-1]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
        work1=as.expression(ccc[[j1,np-1,j2]])
        ccc[[j1,np-1,j2]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
        work1=as.expression(ccc[[np-1,j1,j2]])
        ccc[[np-1,j1,j2]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
      }
    }

    work0=expression(-exp(-g)*((n^2-n)*((digamma(2)-digamma(exp(-g)))^2+trigamma(exp(-g)))+2*n*((digamma(3)-digamma(exp(-g)))^2+trigamma(3)+trigamma(exp(-g))))/((exp(-g)+n)*(exp(-g)+n+1)*b^2))
    work=expression(2*exp(-g)*(((n-1)*(digamma(2)-digamma(exp(-g)))+2*(digamma(3)-digamma(exp(-g))))*fiufi)/((exp(-g)+n)*(exp(-g)+n+1)*b^2))
    work=dpg(as.expression(eval(substitute(substitute(e, list(fiufi=fiufi)), list(e = work[[1]])))))
    work2=expression(-exp(-g)*(fiufiu+fiufi^2)/((exp(-g)+n)*(exp(-g)+n+1)*b^2))
    work2=dpg(as.expression(eval(substitute(substitute(e, list(fiufiu=fiufiu,fiufi=fiufi)), list(e = work2[[1]])))))
    work1=as.expression(ccc[[np,np,np-1]])
    ccc[[np,np,np-1]]=substitute(work0+work+work1+work2,list(work0=work0[[1]],work=work[[1]],work1=work1[[1]],work2=work2[[1]]))
    work1=as.expression(ccc[[np,np-1,np]])
    ccc[[np,np-1,np]]=substitute(work0+work+work1+work2,list(work0=work0[[1]],work=work[[1]],work1=work1[[1]],work2=work2[[1]]))
    work1=as.expression(ccc[[np-1,np,np]])
    ccc[[np-1,np,np]]=substitute(work0+work+work1+work2,list(work0=work0[[1]],work=work[[1]],work1=work1[[1]],work2=work2[[1]]))
    ###########################F3.4
    for (j1 in 1:(np-1)){
      u1=u[j1]
      ufiu1=ufiu[j1]
      ufiufiu1=ufiufiu[j1]
      work0=expression(-(exp(-g)+d)*u1*((digamma(2)-digamma(exp(-g)))^2+trigamma(2)+trigamma(exp(-g)))/((exp(-g)+n)*b^2))
      work0=dpg(as.expression(eval(substitute(substitute(e, list(u1=u1)), list(e = work0[[1]])))))
      work=expression(2*(exp(-g)+d)*((digamma(2)-digamma(exp(-g)))*ufiu1)/((exp(-g)+n)*b^2))
      work=dpg(as.expression(eval(substitute(substitute(e, list(ufiu1=ufiu1)), list(e = work[[1]])))))
      work2=expression(-(exp(-g)+d)*ufiufiu1/((exp(-g)+n)*b^2))
      work2=dpg(as.expression(eval(substitute(substitute(e, list(ufiufiu1=ufiufiu1)), list(e = work2[[1]])))))
      work1=as.expression(ccc[[j1,np,np]])
      ccc[[j1,np,np]]=substitute(work0+work+work1+work2,list(work0=work0[[1]],work=work[[1]],work1=work1[[1]],work2=work2[[1]]))
      work1=as.expression(ccc[[np,j1,np]])
      ccc[[np,j1,np]]=substitute(work0+work+work1+work2,list(work0=work0[[1]],work=work[[1]],work1=work1[[1]],work2=work2[[1]]))
      work1=as.expression(ccc[[np,np,j1]])
      ccc[[np,np,j1]]=substitute(work0+work+work1+work2,list(work0=work0[[1]],work=work[[1]],work1=work1[[1]],work2=work2[[1]]))
      for (j2 in 1:(np-1)){
        uu12=uu[j1,j2]
        uufiu12=uufiu[[j1,j2]]
        work=expression(-(exp(-g)+d)*((digamma(2)-digamma(exp(-g)))*uu12-uufiu12)/((exp(-g)+n)*b))
        work=dpg(as.expression(eval(substitute(substitute(e, list(uu12=uu12,uufiu12=uufiu12)), list(e = work[[1]])))))
        work1=as.expression(ccc[[j1,j2,np]])
        ccc[[j1,j2,np]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
        work1=as.expression(ccc[[j1,np,j2]])
        ccc[[j1,np,j2]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
        work1=as.expression(ccc[[np,j1,j2]])
        ccc[[np,j1,j2]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
        for (j3 in 1:(np-1)){
          uuu123=uuu[j1,j2,j3]
          work=expression(-(exp(-g)+d)*uuu123/(exp(-g)+n))
          work=dpg(as.expression(eval(substitute(substitute(e, list(uuu123=uuu123)), list(e = work[[1]])))))
          work1=as.expression(ccc[[j1,j2,j3]])
          ccc[[j1,j2,j3]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
        }
      }
    }
    work0=expression(-(exp(-g)+d)*n*((digamma(2)-digamma(exp(-g)))^3+3*(trigamma(exp(-g))+trigamma(2))*(digamma(2)-digamma(exp(-g)))+psigamma(2, 2L)-psigamma(exp(-g), 2L))/((exp(-g)+n)*b^3))
    work=expression(3*(exp(-g)+d)*(((digamma(2)-digamma(exp(-g)))^2+(trigamma(2)+trigamma(exp(-g))))*fiufi/((exp(-g)+n)*b^3)))
    work=dpg(as.expression(eval(substitute(substitute(e, list(fiufi=fiufi)), list(e = work[[1]])))))
    work2=expression(-3*(exp(-g)+d)*fiufiu*(digamma(2)-digamma(exp(-g)))/((exp(-g)+n)*b^3))
    work2=dpg(as.expression(eval(substitute(substitute(e, list(fiufiu=fiufiu)), list(e = work2[[1]])))))
    work3=expression((exp(-g)+d)*fiufiufiu/((exp(-g)+n)*b^3))
    work3=dpg(as.expression(eval(substitute(substitute(e, list(fiufiufiu=fiufiufiu)), list(e = work3[[1]])))))
    work1=as.expression(ccc[[np,np,np]])
    ccc[[np,np,np]]=substitute(work0+work+work1+work2+work3,list(work0=work0[[1]],work=work[[1]],work1=work1[[1]],work2=work2[[1]],work3=work3[[1]]))
    ####################################F3.5
    for (j1 in 1:(np-1)){
      ufiufiu1=ufiufiu[j1]
      ufiu1=ufiu[j1]
      u1=u[j1]
      work11=expression(((exp(-g)+d)/((exp(-g)+n)*(exp(-g)+n+1)*b^2))*((n*(digamma(2)-digamma(exp(-g)))-fiufi)*((digamma(2)-digamma(exp(-g)))*u1-ufiu1)+n*trigamma(exp(-g))*u1))
      work11=dpg(as.expression(eval(substitute(substitute(e, list(u1=u1,fiufi=fiufi,ufiu1=ufiu1)), list(e = work11[[1]])))))
      work12=expression(((exp(-g)+d)/((exp(-g)+n)*(exp(-g)+n+1)*b^2))*((n*((digamma(2)-digamma(exp(-g)))^2+trigamma(2)+trigamma(exp(-g)))+fiufiu-2*(digamma(2)-digamma(exp(-g)))*fiufi)*u1))
      work12=dpg(as.expression(eval(substitute(substitute(e, list(u1=u1,fiufi=fiufi,ufiu1=ufiu1,fiufiu=fiufiu)), list(e = work12[[1]])))))
      work21=expression(((exp(-g)+d)/((exp(-g)+n)*(exp(-g)+n+1)*b^2))*(u1*((digamma(2)-digamma(exp(-g)))^2+trigamma(exp(-g)))+ufiufiu1-2*ufiu1*(digamma(2)-digamma(exp(-g)))))
      work21=dpg(as.expression(eval(substitute(substitute(e, list(u1=u1,ufiufiu1=ufiufiu1,ufiu1=ufiu1)), list(e = work21[[1]])))))
      work22=expression(((exp(-g)+d)/((exp(-g)+n)*(exp(-g)+n+1)*b^2))*(u1*((digamma(2)-digamma(exp(-g)))^2+trigamma(2)+trigamma(exp(-g)))+ufiufiu1-2*ufiu1*(digamma(2)-digamma(exp(-g)))))
      work22=dpg(as.expression(eval(substitute(substitute(e, list(u1=u1,ufiufiu1=ufiufiu1,ufiu1=ufiu1)), list(e = work22[[1]])))))
      work3=expression(((exp(-g)+d)/((exp(-g)+n)*(exp(-g)+n+1)*b^2))*2*(u1*((digamma(3)-digamma(exp(-g)))^2+trigamma(3)+trigamma(exp(-g)))+ufiufiu1-2*ufiu1*(digamma(3)-digamma(exp(-g)))))
      work3=dpg(as.expression(eval(substitute(substitute(e, list(ufiufiu1=ufiufiu1,ufiu1=ufiu1,u1=u1)), list(e = work3[[1]])))))
      work=as.expression(ccc[[j1,np,np]])
      ccc[[j1,np,np]]=substitute(work+2*work11-2*work21+work12-work22+3*work3,list(work=work[[1]],work11=work11[[1]],work12=work12[[1]],work21=work21[[1]],work22=work22[[1]],work3=work3[[1]]))
      work=as.expression(ccc[[np,j1,np]])
      ccc[[np,j1,np]]=substitute(work+2*work11-2*work21+work12-work22+3*work3,list(work=work[[1]],work11=work11[[1]],work12=work12[[1]],work21=work21[[1]],work22=work22[[1]],work3=work3[[1]]))
      work=as.expression(ccc[[np,np,j1]])
      ccc[[np,np,j1]]=substitute(work+2*work11-2*work21+work12-work22+3*work3,list(work=work[[1]],work11=work11[[1]],work12=work12[[1]],work21=work21[[1]],work22=work22[[1]],work3=work3[[1]]))
      for (j2 in 1:(np-1)){
        u2=u[j2]
        ufiu1=ufiu[j1]
        ufiu2=ufiu[j2]
        uufiu12=uufiu[[j1,j2]]
        u12=uu[j1,j2]
        work1=expression(((exp(-g)+d)/((exp(-g)+n)*(exp(-g)+n+1)*b))*((n*(digamma(2)-digamma(exp(-g)))-fiufi)*u12))
        work1=dpg(as.expression(eval(substitute(substitute(e, list(u12=u12,fiufi=fiufi)), list(e = work1[[1]])))))
        work2=expression(((exp(-g)+d)/((exp(-g)+n)*(exp(-g)+n+1)*b))*(((digamma(2)-digamma(exp(-g)))*u1-ufiu1)*u2))
        work2=dpg(as.expression(eval(substitute(substitute(e, list(u1=u1,u2=u2,ufiu1=ufiu1)), list(e = work2[[1]])))))
        work3=expression(((exp(-g)+d)/((exp(-g)+n)*(exp(-g)+n+1)*b))*(((digamma(2)-digamma(exp(-g)))*u2-ufiu2)*u1))
        work3=dpg(as.expression(eval(substitute(substitute(e, list(u1=u1,u2=u2,ufiu2=ufiu2)), list(e = work3[[1]])))))
        work=expression(3*((2*digamma(3)-digamma(2)-digamma(exp(-g)))*u12-uufiu12)*(exp(-g)+d)/((exp(-g)+n)*(exp(-g)+n+1)*b))
        work=dpg(as.expression(eval(substitute(substitute(e, list(u12=u12,uufiu12=uufiu12)), list(e = work[[1]])))))
        work4=as.expression(ccc[[j1,j2,np]])
        ccc[[j1,j2,np]]=substitute(work+work1+work2+work3+work4,list(work=work[[1]],work1=work1[[1]],work2=work2[[1]],work3=work3[[1]],work4=work4[[1]]))
        work4=as.expression(ccc[[j1,np,j2]])
        ccc[[j1,np,j2]]=substitute(work+work1+work2+work3+work4,list(work=work[[1]],work1=work1[[1]],work2=work2[[1]],work3=work3[[1]],work4=work4[[1]]))
        work4=as.expression(ccc[[np,j1,j2]])
        ccc[[np,j1,j2]]=substitute(work+work1+work2+work3+work4,list(work=work[[1]],work1=work1[[1]],work2=work2[[1]],work3=work3[[1]],work4=work4[[1]]))
        for (j3 in 1:(np-1)){
          u123=uuu[j1,j2,j3]
          u12=uu[j1,j2]
          u13=uu[j1,j3]
          u23=uu[j2,j3]
          u1=u[j1]
          u2=u[j2]
          u3=u[j3]
          work=expression(((exp(-g)+d)/((exp(-g)+n)*(exp(-g)+n+1)))*(u12*u3+u13*u2+u23*u1+3*u123))
          work=dpg(as.expression(eval(substitute(substitute(e, list(u12=u12,u3=u3,u13=u13,u2=u2,u23=u23,u1=u1,u123=u123)), list(e = work[[1]])))))
          work1=as.expression(ccc[[j1,j2,j3]])
          ccc[[j1,j2,j3]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
        }
      }
    }

    work0=expression(((exp(-g)+d)/((exp(-g)+n)*(exp(-g)+n+1)*b^3))*((n*((digamma(2)-digamma(exp(-g)))^2+trigamma(2)+trigamma(exp(-g)))-2*fiufi*(digamma(2)-digamma(exp(-g)))+fiufiu)*(n*(digamma(2)-digamma(exp(-g)))-fiufi)-n^2*psigamma(exp(-g),2L)+2*n^2*trigamma(exp(-g))*(digamma(2)-digamma(exp(-g)))-2*n*trigamma(exp(-g))*fiufi))
    work0=dpg(as.expression(eval(substitute(substitute(e, list(fiufi=fiufi,fiufiu=fiufiu)), list(e = work0[[1]])))))
    work2=expression(-((exp(-g)+d)/((exp(-g)+n)*(exp(-g)+n+1)*b^3))*(n*(digamma(2)-digamma(exp(-g)))^3-3*(digamma(2)-digamma(exp(-g)))^2*fiufi+3*(digamma(2)-digamma(exp(-g)))*fiufiu-fiufiufiu+n*(3*trigamma(exp(-g))+trigamma(2))*(digamma(2)-digamma(exp(-g)))-fiufi*(trigamma(2)+3*trigamma(exp(-g)))-n*psigamma(exp(-g),2L)))
    work2=dpg(as.expression(eval(substitute(substitute(e, list(fiufi=fiufi,fiufiu=fiufiu,fiufiufiu=fiufiufiu)), list(e = work2[[1]])))))
    work3=expression(2*((exp(-g)+d)/((exp(-g)+n)*(exp(-g)+n+1)*b^3))*(n*(digamma(3)-digamma(exp(-g)))^3-3*(digamma(3)-digamma(exp(-g)))^2*fiufi+3*(digamma(3)-digamma(exp(-g)))*fiufiu-fiufiufiu+3*n*(trigamma(exp(-g))+trigamma(3))*(digamma(3)-digamma(exp(-g)))-3*fiufi*(trigamma(exp(-g))+trigamma(3))+n*psigamma(3,2L)-n*psigamma(exp(-g),2L)))
    work3=dpg(as.expression(eval(substitute(substitute(e, list(fiufi=fiufi,fiufiu=fiufiu,fiufiufiu=fiufiufiu)), list(e = work3[[1]])))))
    work=as.expression(ccc[[np,np,np]])
    ccc[[np,np,np]]=substitute(3*work0+work+3*work2+3*work3,list(work0=work0[[1]],work=work[[1]],work2=work2[[1]],work3=work3[[1]]))
    ######################################################################################F3.6
    for (j1 in 1:(np-1)){
      ufiufiu1=ufiufiu[j1]
      ufiu1=ufiu[j1]
      u1=sum(u[j1])

      work1=expression((-2*(exp(-g)+d)/((exp(-g)+n)*(exp(-g)+n+1)*(exp(-g)+n+2)*b^2))*(u1*(n*(digamma(2)-digamma(exp(-g)))-fiufi)^2+n^2*trigamma(exp(-g))*u1))
      work1=dpg(as.expression(eval(substitute(substitute(e, list(u1=u1,fiufi=fiufi)), list(e = work1[[1]])))))
      work2_1=expression((-2*(exp(-g)+d)/((exp(-g)+n)*(exp(-g)+n+1)*(exp(-g)+n+2)*b^2))*((n*(digamma(2)-digamma(exp(-g)))^2-2*fiufi*(digamma(2)-digamma(exp(-g)))+fiufiu)*u1+n*trigamma(exp(-g))*u1))
      work2_1=dpg(as.expression(eval(substitute(substitute(e, list(u1=u1,fiufi=fiufi,fiufiu=fiufiu)), list(e = work2_1[[1]])))))
      work2_2=expression((-2*(exp(-g)+d)/((exp(-g)+n)*(exp(-g)+n+1)*(exp(-g)+n+2)*b^2))*2*(((digamma(2)-digamma(exp(-g)))*u1-ufiu1)*(n*(digamma(2)-digamma(exp(-g)))-fiufi)+n*trigamma(exp(-g))*u1))
      work2_2=dpg(as.expression(eval(substitute(substitute(e, list(u1=u1,ufiu1=ufiu1,fiufi=fiufi)), list(e = work2_2[[1]])))))
      work2_3=expression((-2*(exp(-g)+d)/((exp(-g)+n)*(exp(-g)+n+1)*(exp(-g)+n+2)*b^2))*2*(u1*(digamma(2)-digamma(exp(-g)))^2-2*ufiu1*(digamma(2)-digamma(exp(-g)))+ufiufiu1+trigamma(exp(-g))*u1))
      work2_3=dpg(as.expression(eval(substitute(substitute(e, list(u1=u1,ufiu1=ufiu1,ufiufiu1=ufiufiu1)), list(e = work2_3[[1]])))))

      work3_1=expression((-2*(exp(-g)+d)/((exp(-g)+n)*(exp(-g)+n+1)*(exp(-g)+n+2)*b^2))*((n*((digamma(3)-digamma(exp(-g)))^2+trigamma(3)+trigamma(exp(-g)))-2*fiufi*(digamma(3)-digamma(exp(-g)))+fiufiu)*u1))
      work3_1=dpg(as.expression(eval(substitute(substitute(e, list(u1=u1,fiufi=fiufi,fiufiu=fiufiu)), list(e = work3_1[[1]])))))
      work3_2=expression((-2*(exp(-g)+d)/((exp(-g)+n)*(exp(-g)+n+1)*(exp(-g)+n+2)*b^2))*((digamma(3)-digamma(exp(-g)))*u1-ufiu1)*(n*(digamma(2)-digamma(exp(-g)))-fiufi))
      work3_2=dpg(as.expression(eval(substitute(substitute(e, list(u1=u1,ufiu1=ufiu1,fiufi=fiufi)), list(e = work3_2[[1]])))))
      work3_3=expression((-2*(exp(-g)+d)/((exp(-g)+n)*(exp(-g)+n+1)*(exp(-g)+n+2)*b^2))*((digamma(3)-digamma(exp(-g)))*u1-ufiu1)*(n*(digamma(2)-digamma(exp(-g)))-fiufi))
      work3_3=dpg(as.expression(eval(substitute(substitute(e, list(u1=u1,ufiu1=ufiu1,fiufi=fiufi)), list(e = work3_3[[1]])))))
      work3_4=expression((-2*(exp(-g)+d)/((exp(-g)+n)*(exp(-g)+n+1)*(exp(-g)+n+2)*b^2))*2*n*u1*trigamma(exp(-g)))
      work3_4=dpg(as.expression(eval(substitute(substitute(e, list(u1=u1)), list(e = work3_4[[1]])))))

      work4_1=expression((-2*(exp(-g)+d)/((exp(-g)+n)*(exp(-g)+n+1)*(exp(-g)+n+2)*b^2))*(((digamma(3)-digamma(exp(-g)))^2+trigamma(3)+trigamma(exp(-g)))*u1-2*ufiu1*(digamma(3)-digamma(exp(-g)))+ufiufiu1))
      work4_1=dpg(as.expression(eval(substitute(substitute(e, list(u1=u1,ufiu1=ufiu1,ufiufiu1=ufiufiu1)), list(e = work4_1[[1]])))))
      work4_2=expression((-2*(exp(-g)+d)/((exp(-g)+n)*(exp(-g)+n+1)*(exp(-g)+n+2)*b^2))*2*(((digamma(3)-digamma(exp(-g)))*(digamma(2)-digamma(exp(-g)))+trigamma(exp(-g)))*u1-ufiu1*(digamma(3)+digamma(2)-2*digamma(exp(-g)))+ufiufiu1))
      work4_2=dpg(as.expression(eval(substitute(substitute(e, list(u1=u1,ufiu1=ufiu1,ufiufiu1=ufiufiu1)), list(e = work4_2[[1]])))))

      work5=expression((-2*(exp(-g)+d)/((exp(-g)+n)*(exp(-g)+n+1)*(exp(-g)+n+2)*b^2))*6*(((digamma(4)-digamma(exp(-g)))^2+trigamma(4)+trigamma(exp(-g)))*u1-2*ufiu1*(digamma(4)-digamma(exp(-g)))+ufiufiu1))
      work5=dpg(as.expression(eval(substitute(substitute(e, list(u1=u1,ufiu1=ufiu1,ufiufiu1=ufiufiu1)), list(e = work5[[1]])))))

      work=as.expression(ccc[[j1,np,np]])
      ccc[[j1,np,np]]=substitute(work+work1-work2_1-work2_2+work2_3+2*work3_1+2*work3_2+2*work3_3+2*work3_4-2*work4_1-2*work4_2+work5,list(work=work[[1]],work1=work1[[1]],work2_1=work2_1[[1]],work2_2=work2_2[[1]],work2_3=work2_3[[1]],work3_1=work3_1[[1]],work3_2=work3_2[[1]],work3_3=work3_3[[1]],work3_4=work3_4[[1]],work4_1=work4_1[[1]],work4_2=work4_2[[1]],work5=work5[[1]]))
      work=as.expression(ccc[[np,j1,np]])
      ccc[[np,j1,np]]=substitute(work+work1-work2_1-work2_2+work2_3+2*work3_1+2*work3_2+2*work3_3+2*work3_4-2*work4_1-2*work4_2+work5,list(work=work[[1]],work1=work1[[1]],work2_1=work2_1[[1]],work2_2=work2_2[[1]],work2_3=work2_3[[1]],work3_1=work3_1[[1]],work3_2=work3_2[[1]],work3_3=work3_3[[1]],work3_4=work3_4[[1]],work4_1=work4_1[[1]],work4_2=work4_2[[1]],work5=work5[[1]]))
      work=as.expression(ccc[[np,np,j1]])
      ccc[[np,np,j1]]=substitute(work+work1-work2_1-work2_2+work2_3+2*work3_1+2*work3_2+2*work3_3+2*work3_4-2*work4_1-2*work4_2+work5,list(work=work[[1]],work1=work1[[1]],work2_1=work2_1[[1]],work2_2=work2_2[[1]],work2_3=work2_3[[1]],work3_1=work3_1[[1]],work3_2=work3_2[[1]],work3_3=work3_3[[1]],work3_4=work3_4[[1]],work4_1=work4_1[[1]],work4_2=work4_2[[1]],work5=work5[[1]]))
      for (j2 in 1:(np-1)){
        u2=u[j2]
        ufiu2=ufiu[j2]
        ufiu12=uufiu[[j1,j2]]
        u12=uu[j1,j2]

        work1=expression((-2*(exp(-g)+d)/((exp(-g)+n)*(exp(-g)+n+1)*(exp(-g)+n+2)*b))*((n*(digamma(2)-digamma(exp(-g)))-fiufi)*u2*u1))
        work1=dpg(as.expression(eval(substitute(substitute(e, list(u1=u1,u2=u2,fiufi=fiufi)), list(e = work1[[1]])))))

        work2_1=expression((-2*(exp(-g)+d)/((exp(-g)+n)*(exp(-g)+n+1)*(exp(-g)+n+2)*b))*((u2*(digamma(2)-digamma(exp(-g)))-ufiu2)*u1))
        work2_1=dpg(as.expression(eval(substitute(substitute(e, list(u1=u1,u2=u2,ufiu2=ufiu2)), list(e = work2_1[[1]])))))

        work2_2=expression((-2*(exp(-g)+d)/((exp(-g)+n)*(exp(-g)+n+1)*(exp(-g)+n+2)*b))*((u1*(digamma(2)-digamma(exp(-g)))-ufiu1)*u2))
        work2_2=dpg(as.expression(eval(substitute(substitute(e, list(u2=u2,ufiu1=ufiu1,u1=u1)), list(e = work2_2[[1]])))))

        work2_3=expression((-2*(exp(-g)+d)/((exp(-g)+n)*(exp(-g)+n+1)*(exp(-g)+n+2)*b))*((n*(digamma(2)-digamma(exp(-g)))-fiufi)*u12))
        work2_3=dpg(as.expression(eval(substitute(substitute(e, list(u12=u12,fiufi=fiufi)), list(e = work2_3[[1]])))))

        work3_1=expression((-2*(exp(-g)+d)/((exp(-g)+n)*(exp(-g)+n+1)*(exp(-g)+n+2)*b))*2*((u2*(digamma(3)-digamma(exp(-g)))-ufiu2)*u1))
        work3_1=dpg(as.expression(eval(substitute(substitute(e, list(u1=u1,ufiu2=ufiu2,u2=u2)), list(e = work3_1[[1]])))))

        work3_2=expression((-2*(exp(-g)+d)/((exp(-g)+n)*(exp(-g)+n+1)*(exp(-g)+n+2)*b))*2*((u1*(digamma(3)-digamma(exp(-g)))-ufiu1)*u2))
        work3_2=dpg(as.expression(eval(substitute(substitute(e, list(u2=u2,ufiu1=ufiu1,u1=u1)), list(e = work3_2[[1]])))))

        work3_3=expression((-2*(exp(-g)+d)/((exp(-g)+n)*(exp(-g)+n+1)*(exp(-g)+n+2)*b))*2*((n*(digamma(2)-digamma(exp(-g)))-fiufi)*u12))
        work3_3=dpg(as.expression(eval(substitute(substitute(e, list(u12=u12,fiufi=fiufi)), list(e = work3_3[[1]])))))

        work3_4=expression((-2*(exp(-g)+d)/((exp(-g)+n)*(exp(-g)+n+1)*(exp(-g)+n+2)*b))*2*((digamma(2)-digamma(exp(-g)))*u12-ufiu12))
        work3_4=dpg(as.expression(eval(substitute(substitute(e, list(u12=u12,ufiu12=ufiu12)), list(e = work3_4[[1]])))))

        work4_1=expression((-2*(exp(-g)+d)/((exp(-g)+n)*(exp(-g)+n+1)*(exp(-g)+n+2)*b))*4*((digamma(3)-digamma(exp(-g)))*u12-ufiu12))
        work4_1=dpg(as.expression(eval(substitute(substitute(e, list(u12=u12,ufiu12=ufiu12)), list(e = work4_1[[1]])))))

        work4_2=expression((-2*(exp(-g)+d)/((exp(-g)+n)*(exp(-g)+n+1)*(exp(-g)+n+2)*b))*2*((digamma(2)-digamma(exp(-g)))*u12-ufiu12))
        work4_2=dpg(as.expression(eval(substitute(substitute(e, list(u12=u12,ufiu12=ufiu12)), list(e = work4_2[[1]])))))

        work5=expression((-2*(exp(-g)+d)/((exp(-g)+n)*(exp(-g)+n+1)*(exp(-g)+n+2)*b))*6*((digamma(4)-digamma(exp(-g)))*u12-ufiu12))
        work5=dpg(as.expression(eval(substitute(substitute(e, list(u12=u12,ufiu12=ufiu12)), list(e = work5[[1]])))))

        work4=as.expression(ccc[[j1,j2,np]])
        ccc[[j1,j2,np]]=substitute(work1-work2_1-work2_2-work2_3+work3_1+work3_2+work3_3+work3_4-work4_1-work4_2+work5+work4,list(work1=work1[[1]],work2_1=work2_1[[1]],work2_2=work2_2[[1]],work2_3=work2_3[[1]],work3_1=work3_1[[1]],work3_2=work3_2[[1]],work3_3=work3_3[[1]],work3_4=work3_4[[1]],work4_1=work4_1[[1]],work4_2=work4_2[[1]],work5=work5[[1]],work4=work4[[1]]))
        work4=as.expression(ccc[[j1,np,j2]])
        ccc[[j1,np,j2]]=substitute(work1-work2_1-work2_2-work2_3+work3_1+work3_2+work3_3+work3_4-work4_1-work4_2+work5+work4,list(work1=work1[[1]],work2_1=work2_1[[1]],work2_2=work2_2[[1]],work2_3=work2_3[[1]],work3_1=work3_1[[1]],work3_2=work3_2[[1]],work3_3=work3_3[[1]],work3_4=work3_4[[1]],work4_1=work4_1[[1]],work4_2=work4_2[[1]],work5=work5[[1]],work4=work4[[1]]))
        work4=as.expression(ccc[[np,j1,j2]])
        ccc[[np,j1,j2]]=substitute(work1-work2_1-work2_2-work2_3+work3_1+work3_2+work3_3+work3_4-work4_1-work4_2+work5+work4,list(work1=work1[[1]],work2_1=work2_1[[1]],work2_2=work2_2[[1]],work2_3=work2_3[[1]],work3_1=work3_1[[1]],work3_2=work3_2[[1]],work3_3=work3_3[[1]],work3_4=work3_4[[1]],work4_1=work4_1[[1]],work4_2=work4_2[[1]],work5=work5[[1]],work4=work4[[1]]))
        for (j3 in 1:(np-1)){
          u123=uuu[j1,j2,j3]
          u12=uu[j1,j2]
          u13=uu[j1,j3]
          u23=uu[j2,j3]
          u1=u[j1]
          u2=u[j2]
          u3=u[j3]
          work=expression((-2*(exp(-g)+d)/((exp(-g)+n)*(exp(-g)+n+1)*(exp(-g)+n+2)))*(u1*u2*u3+u12*u3+u13*u2+u23*u1+2*u123))
          work=dpg(as.expression(eval(substitute(substitute(e, list(u12=u12,u13=u13,u23=u23,u1=u1,u2=u2,u3=u3,u123=u123)), list(e = work[[1]])))))
          work1=as.expression(ccc[[j1,j2,j3]])
          ccc[[j1,j2,j3]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
        }
      }
    }
    work=expression((-2*(exp(-g)+d)/((exp(-g)+n)*(exp(-g)+n+1)*(exp(-g)+n+2)*b^3))*((n*(digamma(2)-digamma(exp(-g)))-fiufi)^3+3*n^2*trigamma(exp(-g))*(n*(digamma(2)-digamma(exp(-g)))-fiufi)-n^3*psigamma(exp(-g),2L)))
    work=dpg(as.expression(eval(substitute(substitute(e, list(fiufi=fiufi)), list(e = work[[1]])))))
    work1=expression((-2*(exp(-g)+d)/((exp(-g)+n)*(exp(-g)+n+1)*(exp(-g)+n+2)*b^3))*((n*(digamma(2)-digamma(exp(-g)))-fiufi)*(n*(digamma(2)-digamma(exp(-g)))^2-2*fiufi*(digamma(2)-digamma(exp(-g)))+fiufiu)+3*n*trigamma(exp(-g))*(n*(digamma(2)-digamma(exp(-g)))-fiufi)-n^2*psigamma(exp(-g),2L)))
    work1=dpg(as.expression(eval(substitute(substitute(e, list(fiufi=fiufi,fiufiu=fiufiu)), list(e = work1[[1]])))))
    work2=expression((-2*(exp(-g)+d)/((exp(-g)+n)*(exp(-g)+n+1)*(exp(-g)+n+2)*b^3))*(n*(digamma(2)-digamma(exp(-g)))^3-3*(digamma(2)-digamma(exp(-g)))^2*fiufi+3*(digamma(2)-digamma(exp(-g)))*fiufiu-fiufiufiu+3*trigamma(exp(-g))*(n*(digamma(2)-digamma(exp(-g)))-fiufi)-n*psigamma(exp(-g),2L)))
    work2=dpg(as.expression(eval(substitute(substitute(e, list(fiufi=fiufi,fiufiu=fiufiu,fiufiufiu=fiufiufiu)), list(e = work2[[1]])))))
    work2_1=expression((-2*(exp(-g)+d)/((exp(-g)+n)*(exp(-g)+n+1)*(exp(-g)+n+2)*b^3))*(n*(digamma(2)-digamma(exp(-g)))-fiufi)*(n*(digamma(3)-digamma(exp(-g)))^2-2*fiufi*(digamma(3)-digamma(exp(-g)))+fiufiu))
    work2_1=dpg(as.expression(eval(substitute(substitute(e, list(fiufi=fiufi,fiufiu=fiufiu)), list(e = work2_1[[1]])))))
    work2_2=expression((-2*(exp(-g)+d)/((exp(-g)+n)*(exp(-g)+n+1)*(exp(-g)+n+2)*b^3))*(2*n*trigamma(exp(-g))*(n*(digamma(3)-digamma(exp(-g)))-fiufi)-n^2*psigamma(exp(-g),2L)))
    work2_2=dpg(as.expression(eval(substitute(substitute(e, list(fiufi=fiufi)), list(e = work2_2[[1]])))))
    work2_3=expression((-2*(exp(-g)+d)/((exp(-g)+n)*(exp(-g)+n+1)*(exp(-g)+n+2)*b^3))*n*(trigamma(exp(-g))+trigamma(3))*(n*(digamma(2)-digamma(exp(-g)))-fiufi))
    work2_3=dpg(as.expression(eval(substitute(substitute(e, list(fiufi=fiufi)), list(e = work2_3[[1]])))))
    work3_1=expression((-2*(exp(-g)+d)/((exp(-g)+n)*(exp(-g)+n+1)*(exp(-g)+n+2)*b^3))*(n*((digamma(3)-digamma(exp(-g)))^2*(digamma(2)-digamma(exp(-g)))+2*trigamma(exp(-g))*(digamma(3)-digamma(exp(-g)))-psigamma(exp(-g),2L)+(trigamma(3)+trigamma(exp(-g)))*(digamma(2)-digamma(exp(-g))))))
    work3_2=expression((-2*(exp(-g)+d)/((exp(-g)+n)*(exp(-g)+n+1)*(exp(-g)+n+2)*b^3))*(-fiufi*((digamma(3)-digamma(exp(-g)))^2+2*(digamma(3)-digamma(exp(-g)))*(digamma(2)-digamma(exp(-g)))+(3*trigamma(exp(-g))+trigamma(3)))+fiufiu*(digamma(2)-digamma(exp(-g)))+2*fiufiu*(digamma(3)-digamma(exp(-g)))-fiufiufiu))
    work3_2=dpg(as.expression(eval(substitute(substitute(e, list(fiufi=fiufi,fiufiufiu=fiufiufiu,fiufiu=fiufiu,fiufi=fiufi)), list(e = work3_2[[1]])))))
    work4_1=expression((-2*(exp(-g)+d)/((exp(-g)+n)*(exp(-g)+n+1)*(exp(-g)+n+2)*b^3))*(n*(digamma(4)-digamma(exp(-g)))^3-3*(digamma(4)-digamma(exp(-g)))^2*fiufi+3*(digamma(4)-digamma(exp(-g)))*fiufiu-fiufiufiu))
    work4_1=dpg(as.expression(eval(substitute(substitute(e, list(fiufi=fiufi,fiufiufiu=fiufiufiu,fiufiu=fiufiu,fiufi=fiufi)), list(e = work4_1[[1]])))))
    work4_2=expression((-2*(exp(-g)+d)/((exp(-g)+n)*(exp(-g)+n+1)*(exp(-g)+n+2)*b^3))*(3*(n*(digamma(4)-digamma(exp(-g)))-fiufi)*(trigamma(4)+trigamma(exp(-g)))+n*(psigamma(4,2L)-psigamma(exp(-g),2L))))
    work4_2=dpg(as.expression(eval(substitute(substitute(e, list(fiufi=fiufi)), list(e = work4_2[[1]])))))

    work0=as.expression(ccc[[np,np,np]])
    ccc[[np,np,np]]=substitute(work0+work-3*work1+2*work2+6*work2_1+6*work2_2+6*work2_3-6*work3_1-6*work3_2+6*work4_1+6*work4_2,list(work0=work0[[1]],work=work[[1]],work1=work1[[1]],work2=work2[[1]],work2_1=work2_1[[1]],work2_2=work2_2[[1]],work2_3=work2_3[[1]],work3_1=work3_1[[1]],work3_2=work3_2[[1]],work4_1=work4_1[[1]],work4_2=work4_2[[1]]))
    Rccc=array(NA,c(np,np,np))
    for (i1 in 1:np){
      for (i2 in 1:np){
        for (i3 in 1:np){
          #        print(c(i1,i2,i3))
          Rccc[[i1,i2,i3]]=eval(ccc[[i1,i2,i3]])
        }
      }
    }
    ###########################Final
    Rccc[nf+2,nf+2,nf+2]=Rccc[nf+2,nf+2,nf+2]+exp(-g)*sum(1/(exp(-g)+(0:(n-1))))
    RCCC=RCCC+Rccc
    RCC=RCC+Rcc
    RCC_C=RCC_C+Rcc_c
  }
    RCC_=ginv(RCC)
    kra=RCC_C-0.5*RCCC
    if (!(any(is.na(RCC)) | any(is.na(RCCC)))){
    RCC_=ginv(RCC)
    Bias={}
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
