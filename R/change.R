#' @title change
#'
#' @description   Adds values of i1, i2 and i3 to 8th, 9th and 10th columns of character matrix 'Text'.
#'
#' @param Text A character matrix, where 2-10 columns are numbers in character format.
#' @param i1 A number.
#' @param i2 A number.
#' @param i3 A number.
#'
#' @return Changed character matrix 'Text'.
#'
#' @examples
#' \dontrun{
#' change(Text,i1,i2,i3)
#' }
#'
#' @export
#'
change=function(Text,i1,i2,i3){
  Text[,8]=as.character(as.numeric(Text[,8])+i1)
  Text[,9]=as.character(as.numeric(Text[,9])+i2)
  Text[,10]=as.character(as.numeric(Text[,10])+i3)
  return(Text)
}
