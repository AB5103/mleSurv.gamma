#' @title BiasCorrectionE
#'
#' @description   Returns the bias correction for ML estimates.
#'
#' @param D The data frame. The first 'nf' columns are the covariates under study. Other columns must include fields
#' 'time' for time-to-failure, 'Cens' for censoring and 'cluster' for cluster identification. The number of rows is equal to the number of subjects.
#' @param nf The number of covariates under study.
#' @param ncl The number of clusters.
#' @param para  The vector of the ML estimates of unknown parameters. This vector has a form
#' \deqn{y = (ln(a), \beta , ln(\sigma ^2)),}
#' where \eqn{a} is the slope, \eqn{\beta } is the vector of the Cox-regression parameters,
#' and \eqn{\sigma ^2} is the variance of frailty.
#' @details This function calculates the bias correction for the proportional hazards model with gamma frailty
#' in accordance with the Cox-Snell method.  The baseline cumulative hazard function is \deqn{H_{0}(t;a,b)=at}.
#'
#' @return The bias correction for vector \eqn{y}. If the calculation of the bias correction is
#' not possible returns NAs.
#'
#' @import MASS
#'
#' @examples
#' \dontrun{
#' BiasCorrectionE(D,nf,ncl,para)
#'
#'}
#' @export
#'
BiasCorrectionE=function(D,nf,ncl,para){
  parnam=c("a",paste0("beta_",as.character(1:nf)),"g")
  nr=nrow(D)
  np=nf+2
  a.scale=rep(1,nr)
  g.scale=rep(1,nr)
  uR_beta=D[,1:nf]
  UR=cbind(a.scale,uR_beta,g.scale)
  Cens=D$event
  a=para[1]
  g=para[nf+2]
  par=c(a,para[2:(1+nf)],g)
  parnam=c("a",paste0("beta_",as.character(1:nf)),"g")
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

    Part3=expression(lgamma(exp(-g)+d)-lgamma(exp(-g)))
    Part3_2=as.expression(DD(Part3,"g",2))
    Part3_3=as.expression(DD(Part3,"g",3))
    ##########################################
    #Second order derivative
    ##########################################
    cc=array(expression(nul),c(np,np))
    ###########################F2.1
    sumug=u[np]
    work=expression(exp(-g)*sumug/(exp(-g)+n))
    work=dpg(as.expression(eval(substitute(substitute(e, list(sumug=sumug)), list(e = work[[1]])))))
    cc[[np,np]]=substitute(Part3_2+2*work,list(Part3_2=Part3_2[[1]],work=work[[1]]))
    for (i in 1:(np-1)){
      sumui=u[i]
      work=expression(exp(-g)*sumui/(exp(-g)+n))
      work=dpg(as.expression(eval(substitute(substitute(e, list(sumui=sumui)), list(e = work[[1]])))))
      cc[[np,i]]=substitute(work,list(work=work[[1]]))
      cc[[i,np]]=cc[[np,i]]
    }
    ###########################F2.2
    for (j1 in 1:np){
      u1=u[j1]
      for (j2 in 1:np){
        u12=uu[[j1,j2]]
        work=expression(-(exp(-g)+d)*u12/(exp(-g)+n))
        work=dpg(as.expression(eval(substitute(substitute(e, list(u12=u12)), list(e = work[[1]])))))
        work1=as.expression(cc[[j1,j2]])
        cc[[j1,j2]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
      }
    }
    ###########################F2.3
    for (j1 in 1:np){
      u1=u[j1]
      for (j2 in 1:np){
        u2=u[j2]
        u12=uu[j1,j2]
        work=expression(((exp(-g)+d)/((exp(-g)+n)*(exp(-g)+n+1)))*(u12+u1*u2))
        work=dpg(as.expression(eval(substitute(substitute(e, list(u1=u1,u2=u2,u12=u12)), list(e = work[[1]])))))
        work1=as.expression(cc[[j1,j2]])
        cc[[j1,j2]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
      }
    }
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
    Rcc[np,np]=Rcc[np,np]-exp(-g)*sum(1/(exp(-g)+(0:(n-1))))
    Rcc_c[np,np,np]=Rcc_c[np,np,np]+exp(-g)*sum(1/(exp(-g)+(0:(n-1))))-exp(-2*g)*sum(1/((exp(-g)+(0:(n-1)))^2))
    ##########################################
    #Third order derivative
    ##########################################
    ###########################F3.1
    ccc=array(expression(nul),c(np,np,np))
    sumug=u[np]
    work=expression(-exp(-g)*sumug/(exp(-g)+n))
    work=dpg(as.expression(eval(substitute(substitute(e, list(sumug=sumug)), list(e = work[[1]])))))
    ccc[[np,np,np]]=substitute(Part3_3+3*work,list(Part3_3=Part3_3[[1]],work=work[[1]]))
    for (i in 1:(np-1)){
      sumui=u[i]
      work=expression(-exp(-g)*sumui/(exp(-g)+n))
      work=dpg(as.expression(eval(substitute(substitute(e, list(sumui=sumui)), list(e = work[[1]])))))
      ccc[[np,np,i]]=substitute(work,list(work=work[[1]]))
      ccc[[i,np,np]]=ccc[[np,np,i]]
      ccc[[np,i,np]]=ccc[[np,np,i]]
    }
    ###########################F3.2
    for (j1 in 1:np){
      u1=u[j1]
      for (j2 in 1:np){
        u12=uu[j1,j2]
        work=expression(exp(-g)*u12/(exp(-g)+n))
        work=dpg(as.expression(eval(substitute(substitute(e, list(u12=u12)), list(e = work[[1]])))))
        work1=as.expression(ccc[[np,j1,j2]])
        ccc[[np,j1,j2]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
        work1=as.expression(ccc[[j1,np,j2]])
        ccc[[j1,np,j2]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
        work1=as.expression(ccc[[j1,j2,np]])
        ccc[[j1,j2,np]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
      }
    }
    ###########################F3.3
    for (j1 in 1:np){
      u1=u[j1]
      for (j2 in 1:np){
        u2=u[j2]
        u12=uu[j1,j2]
        work=expression((-exp(-g)/((exp(-g)+n)*(exp(-g)+n+1)))*(u12+u1*u2))
        work=dpg(as.expression(eval(substitute(substitute(e, list(u1=u1,u2=u2,u12=u12)), list(e = work[[1]])))))
        work1=as.expression(ccc[[j1,j2,np]])
        ccc[[j1,j2,np]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
        work1=as.expression(ccc[[j1,np,j2]])
        ccc[[j1,np,j2]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
        work1=as.expression(ccc[[np,j1,j2]])
        ccc[[np,j1,j2]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
      }
    }
    ###########################F3.4
    for (j1 in 1:np){
      for (j2 in 1:np){
        for (j3 in 1:np){
          uuu123=uuu[j1,j2,j3]
          work=expression(-(exp(-g)+d)*uuu123/(exp(-g)+n))
          work=dpg(as.expression(eval(substitute(substitute(e, list(uuu123=uuu123)), list(e = work[[1]])))))
          work1=as.expression(ccc[[j1,j2,j3]])
          ccc[[j1,j2,j3]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
        }
      }
    }
    ####################################F3.5
    for (j1 in 1:np){
      u1=u[j1]
      for (j2 in 1:np){
        u2=u[j2]
        for (j3 in 1:np){
          u123=uuu[j1,j2,j3]
          u12=uu[j1,j2]
          u13=uu[j1,j3]
          u23=uu[j2,j3]
          u3=u[j3]
          work=expression(((exp(-g)+d)/((exp(-g)+n)*(exp(-g)+n+1)))*(u12*u3+u13*u2+u23*u1+3*u123))
          work=dpg(as.expression(eval(substitute(substitute(e, list(u12=u12,u3=u3,u13=u13,u2=u2,u23=u23,u1=u1,u123=u123)), list(e = work[[1]])))))
          work1=as.expression(ccc[[j1,j2,j3]])
          ccc[[j1,j2,j3]]=substitute(work+work1,list(work=work[[1]],work1=work1[[1]]))
        }
      }
    }
    ####################################F3.6
    for (j1 in 1:np){
      u1=sum(u[j1])
      for (j2 in 1:np){
        u2=u[j2]
        u12=uu[j1,j2]

        for (j3 in 1:np){
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
    Rccc=array(NA,c(np,np,np))
    for (i1 in 1:np){
      for (i2 in 1:np){
        for (i3 in 1:np){
          Rccc[[i1,i2,i3]]=eval(ccc[[i1,i2,i3]])
        }
      }
    }
    ###########################Final
    Rccc[np,np,np]=Rccc[np,np,np]+exp(-g)*sum(1/(exp(-g)+(0:(n-1))))
    RCCC=RCCC+Rccc
    RCC=RCC+Rcc
    RCC_C=RCC_C+Rcc_c
  }
  if (!(any(is.na(RCC)) | any(is.na(RCCC)))){
  RCC_=ginv(RCC)
  kra=RCC_C-0.5*RCCC
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
