#' To Numeric.
#'
#' This function checks if your vector is not a character and will then converts it to a numeric.
#'
#'
#' @rdname toNum
#' @name toNum
#'
#' @param x a vector of character, factor, or numeric
#'
#' @return \code{toNum} Returns a numeric vector. \cr
#'
#'
#' @examples
#' 
#'  
#'  # vector of factors
#'  fact <- C("1", "2", "3", "4"))
#'  # convert to character
#'  num <- toNum(fact)
#'  
#' @export


toNum <- function(x) {
  # this function checks if your vector is not a numeric and will then convert it to a numeric
  # args
  # x = vector that is character, factor, or numeric
  # Returns
  # Vector as a numeric 
  if(class(x) == "character") return(as.numeric(x))
  
  if(class(x) == "factor") return(as.numeric(levels(x))[x])
  if(class(x) == "numeric") return(x)
} 

