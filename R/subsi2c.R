#' @title subsi2c
#'
#' @description   Calculates the sum of elements for the matrix of expressions.
#'
#' @param Zag A character object. Includes variables \eqn{x_i} and \eqn{x__i} corresponding to covariates.
#' @param nf The number of covariates.
#' @param Parnam The vector of the names of parameters needed to calculate expression.
#' @param Parmen The vector of the values of parameters needed to calculate expression.
#' @param WW4 The vector of the matrices of covariates. Each such matrix is the \eqn{n}-times repeated vector-column
#'            standing for a covariate, where \eqn{n} is the number of subjects.
#' @param WW5 The vector of the transposed matrices of covariates.
#' @details Text 'Zag' stands for an expression. Variables \eqn{x_i} and \eqn{x__i} are replaced
#' with matrix of covariates and its transposed, respectively. Then the full sum of expression over all subjects is calculated.
#' This function is used in function 'BiasCorrectionG' only.
#'
#' @return The full sum of expression over all subjects
#'
#' @examples
#' \dontrun{
#' subsi2c(Zag,nf,Parnam,Parmen,WW4,WW5)
#' }
#'
#' @export
#'
subsi2c=function(Zag,nf,Parnam,Parmen,WW4,WW5){

  ni=length(Parnam)
  for(ik in 1:ni){
  assign(Parnam[ik],Parmen[ik])
  }
  assign("n",Parmen[which(Parnam=="n")])
  assign("Deg0",Parmen[which(Parnam=="Deg0")])

  for (ik in 0:(nf+1)){
    assign(paste0("U2_",as.character(ik+1)), WW4[,(n*ik+1):(n*(ik+1))])
    assign(paste0("tU2_",as.character(ik+1)),WW5[,(n*ik+1):(n*(ik+1))])
  }

  Zag=paste0("sum(",Zag,")")
  for (i in 1:(nf+2)){
    Zag=gsub(paste0("x_",as.character(i)),paste0("U2_",as.character(i)),Zag)
    Zag=gsub(paste0("x__",as.character(i)),paste0("tU2_",as.character(i)),Zag)
  }
  Zagy=eval(parse(text=Zag))

  return(Zagy)
}
