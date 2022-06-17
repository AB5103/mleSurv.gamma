#' @title convert
#'
#' @description   Calculates the sum of values over expressions.
#'
#' @param Text A character matrix, where the first column contains expressions in text format.
#' Other 9 columns contain non-negative integer numbers.
#' @param i1 first index.
#' @param i2 second index.
#' @param Parnam The vector of the names of parameters needed to calculate expressions.
#' @param Parmen The vector of the values of parameters needed to calculate expressions.
#' @details For each row of the matrix 'Text' this function calculates the value of the following expression:
#' \deqn{Expr=n^Deg*T_1*c_1^(T_2*ND1+T_5)*la1_i1_(T_2*ND1+T_5)_T8*c_2^(T_3*ND2+T_6)*la2_i2_(T_3*ND1+T_6)_T9*c_3^(T_4*ND1+T_7)*la3_0_(T_4*ND1+T_7)_T10}
#' and then takes the sum over all rows. Here \eqn{T_1} (first column in parameter 'Text') is a character string standing for an expression,
#' \eqn{T_1}, \eqn{i=2,...,10} (last 9 columns in parameter 'Table') are the numbers in character format,
#' other variables (with exception of \eqn{i1} and \eqn{i2}) are the global parameters.
#' @return The sum of values over all expressions in the first column 'Text'.
#'
#' @examples
#' \dontrun{
#' convert(Text,i1,i2)
#' }
#'
#' @export
#'
convert=function(Text,i1,i2,Parnam,Parmen){
  assign("n",Parmen[which(Parnam=="n")])
  ni=length(Parnam)
  assign("Deg0",Parmen[which(Parnam=="Deg0")])
  assign("ND1",Parmen[which(Parnam=="ND1")])
  assign("ND2",Parmen[which(Parnam=="ND2")])
  assign("ND3",Parmen[which(Parnam=="ND3")])
  for(ik in 1:ni){
    assign(Parnam[ik],Parmen[ik])
  }
  N1=((as.numeric(Text[,2])+as.numeric(Text[,5])+as.numeric(Text[,8]))>0)
  N2=((as.numeric(Text[,3])+as.numeric(Text[,6])+as.numeric(Text[,9]))>0)
  N3=((as.numeric(Text[,4])+as.numeric(Text[,7])+as.numeric(Text[,10]))>0)
  if (i1!=0) N1=((N1+1)>0)
  if (i2!=0) N2=((N2+1)>0)
  Deg=as.character(Deg0-N1-N2-N3)
  work=paste0("(n^",Deg,"*(",Text[,1],"))*((c1^(",as.character(as.numeric(Text[,2])*ND1+as.numeric(Text[,5])),"))*la1_",as.character(i1),"_",as.character(as.numeric(Text[,2])*ND1+as.numeric(Text[,5])),"_",Text[,8],"*(c2^(",as.character(as.numeric(Text[,3])*ND2+as.numeric(Text[,6])),"))*la2_",as.character(i2),"_",as.character(as.numeric(Text[,3])*ND2+as.numeric(Text[,6])),"_",Text[,9],"*(c3^(",as.character(as.numeric(Text[,4])*ND3+as.numeric(Text[,7])),"))*la3_0_",as.character(as.numeric(Text[,4])*ND3+as.numeric(Text[,7])),"_",Text[,10],")")
  R=0
  for (i in 1:nrow(Text)){
    R=R+eval(parse(text=work[i]))
  }
  return(R)
}
