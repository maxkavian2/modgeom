#'@title Finding the interval index of a knot support
#'@description \code{bspline_find_interval_index} finds the interval
#'index of the provided parameter value.
#'@details The evaluation of B-splines at the parameter value \code{x} is bound to the
#'interval where \code{x} lies within the B-spline support (i.e. the sequence of
#'ruputure points or knots). This function computes this interval for the single
#'value of \code{x}.
#'@param x The parameter value (single value).
#'@param U The support, as a vector.
#'@param native If \code{TRUE} native code will be used instead (default).
#'@return the interval index, as an integer.
#'@author Máximo Sánchez-Aragón
bspline_find_interval_index <- function(x , U, native=TRUE){

  if(!native){
    
    r = 0; # 0 index would generate an error (see below)
    for(i in 1:length(U)-1)
      if(x >= U[i] && x < U[i+1]){ # this condition cannot be satisfied for empty intervals
        r=i
        break
      }
  
    # adjusts r to the first not null interval
    i=0
    if(r == 0 && x >= U[length(U)])
      while(U[length(U)-i] == U[length(U)-i-1]){
        i = i +1
        r = r - 1
      }
    return(r)
    
  }else{
    # here goes the native code (check bspline_core.cpp file)
    errcode <- 0
    r <- 0
    result <- .C("bspline_find_interval_index", as.numeric(x),
                                      as.numeric(U),
                                      as.integer(length(U)),
                                      as.integer(r),
                                      as.integer(errcode))
    return(result[[4]])
  }

}
