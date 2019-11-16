% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/toNum.R
\name{toNum}
\alias{toNum}
\title{To Numeric.}
\usage{
toNum(x)
}
\arguments{
\item{x}{a vector of character, factor, or numeric}
}
\value{
\code{toNum} Returns a numeric vector. \cr
}
\description{
This function checks if your vector is not a character and will then converts it to a numeric.
}
\examples{

 
 # vector of factors
 fact <- C("1", "2", "3", "4"))
 # convert to character
 num <- toNum(fact)
 
}