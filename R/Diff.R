#' @title Diff
#'
#' @description   Differentiates and modifies (if necessary) an expression given in variable 'Text'.
#'
#' @param Text A character matrix, where the first column contains expressions in text format that must be differentiated.
#' Other 9 columns contain non-negative integer numbers.
#' @param x character vector, giving the variable name including numerical index (e.g. 'x1') with respect to which derivatives will be computed.
#' @param order Order of differentiation. A positive integer number.
#' @details This function is written specially for this package and takes into account the special structure of matrix 'Text'.
#' Modification of the matrix 'Text' depends on the content of the columns 2-10 of matrix 'Text' and on the index part of variable 'x'.
#'
#' @return Changed character matrix 'Text'.
#'
#' @import Deriv
#'
#' @examples
#' \dontrun{
#' Text=matrix(c("x1^2","x1^3",rep(c("2","1"),9)),2,10)
#' Diff(Text,"x1",1)
#' }
#'
#' @export
#'
Diff=function(Text,x,order=1){
  if(order < 1) stop("'order' must be >= 1")
  if(order == 1) {
    Text0=Text
    Text=Diff0(Text,x)
    iw=as.numeric(substring(x,2,2))+1
    ns=nrow(Text0)
    count=0
    for (i in 1:ns){
      if (as.numeric(Text0[i,iw])>0){
        count=count+1
        P2=Text0[i,]
        P2[as.numeric(substring(x,2,2))+7]=as.character(as.numeric(P2[as.numeric(substring(x,2,2))+7])+1)
        Text=rbind(Text,P2)
      }
    }
    return(Text)
  } else {
    return(Diff(Diff0(Text,x),x,order-1))
  }
}
