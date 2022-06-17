#' @title Diff0
#'
#' @description   Provides derivatives of first order.
#'
#' @param Text A character matrix, where the first column contains expressions in text format that must be differentiated.
#' Other 9 columns contain non-negative integer numbers.
#' @param x A character vector, giving the variable name  with respect to which derivatives will be computed.
#'
#'
#' @return Diff0 return a call.
#'
#' @import Deriv
#'
#' @examples
#' \dontrun{
#' Text=matrix(c("x^2","x^3",rep(c("2","1"),9)),2,10)
#' Diff0(Text,"x")
#' }
#'
#' @export
#'
Diff0=function(Text,x) {
  ns=nrow(Text)
  for (i in 1:ns) {
    Text[i,1]=paste0(deparse(DD(parse(text=Text[i,1]),x,1)),collapse="")
  }
  return(Text)
}
