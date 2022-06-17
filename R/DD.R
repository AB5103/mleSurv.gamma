#' @title DD
#'
#' @description Symbolic differentiation of any order of simple expressions.
#' @param expr  An expression.
#' @param name  Character vector, giving the variable names (only one for DD()) with respect to which derivatives will be computed.
#' @param order Order of differentiation. A positive integer number.
#'
#' @return The result of differentiation (an expession object).
#'
#' @import Deriv
#'
#' @examples
#' \dontrun{
#' DD(expression(x^2+y^3),y,2)
#' }
#'
#' @export
#'
DD <- function(expr, name, order = 1) {
  if(order < 1) stop("'order' must be >= 1")
  if(order == 1) D(expr, name)
  else DD(D(expr, name), name, order - 1)
}
