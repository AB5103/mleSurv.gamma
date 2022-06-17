#' @title subsi2c
#'
#' @description   Calculates the sum of elements for the matrix of expressions.
#'
#' @param Zag A character object. Includes variables \eqn{x_i} and \eqn{x__i} corresponding to covariates.
#' @param nf The number of covariates.
#' @param Parnam The vector of the names of parameters needed to calculate expression.
#' @param Parmen The vector of the values of parameters needed to calculate expression.
#' @param WW6 The matrix of covariates.
#' @details Text 'Zag' is a character representation of an expression. Variables \eqn{x_i} and \eqn{x__i} are replaced
#' with respective matrices of covariates and the full sum of expression over all subjects is calculated.
#' This function is used in function 'BiasCorrectionG' only.
#'
#' @return The full sum of expression over all subjects
#'
#' @examples
#' \dontrun{
#' subsi2c_i(Zag,nf,Parnam,Parmen,WW6)
#' }
#'
#' @export
#'
subsi2c_i=function(Zag,nf,Parnam,Parmen,WW6){

  ni=length(Parnam)
  for(ik in 1:ni){
    assign(Parnam[ik],Parmen[ik])
  }

  for (ik in 1:(nf+2)){
    assign(paste0("U1_",as.character(ik)), WW6[,ik])
  }

  Zag=paste0("sum(",Zag,")")
  for (i in 1:(nf+2)){
    Zag=gsub(paste0("x_",as.character(i)),paste0("U1_",as.character(i)),Zag)
    Zag=gsub(paste0("x__",as.character(i)),paste0("U1_",as.character(i)),Zag)
  }
  Zagy=eval(parse(text=Zag))
  return(Zagy)
}
