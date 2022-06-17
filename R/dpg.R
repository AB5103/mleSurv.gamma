#' @title dpg
#'
#' @description   Returns the parsed but unevaluted expression.
#'
#' @param t An expression.
#' @details Convert expression 't' into parsed object.    
#'
#' @return Parsed object.
#'
#' @examples
#' \dontrun{
#' t=expression(a+b*x)
#' dpg(t) 
#' }
#'
#' @export
#'
dpg=function(t) parse(text=gsub("expression","",deparse(t)))
