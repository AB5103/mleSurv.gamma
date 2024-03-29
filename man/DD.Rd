% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DD.R
\name{DD}
\alias{DD}
\title{DD}
\usage{
DD(expr, name, order = 1)
}
\arguments{
\item{expr}{An expression.}

\item{name}{Character vector, giving the variable names (only one for DD()) with respect to which derivatives will be computed.}

\item{order}{Order of differentiation. A positive integer number.}
}
\value{
The result of differentiation (an expession object).
}
\description{
Symbolic differentiation of any order of simple expressions.
}
\examples{
\dontrun{
DD(expression(x^2+y^3),y,2)
}

}
