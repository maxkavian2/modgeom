#'@title Building a clamped knot support
#'@description \code{bspline_support} makes a clamped knot support for splines
#'@details This function makes a clamped knot support out of an increasing sequence of
#'knots. The support creates the necessary multiplicities (but not all)
#'for the provided degree (i.e. the degree of the polynomials that build the spline).
#'These multiplicities determine at last the maximum curve class for the spline.
#'@param degree The degree of the polynomial for the current spline.
#'@param knots The non-redundant sequence of knots (rupture points) in the parameter space.
#'@param multiplicities A vector of integers, with same length as \code{knots}
#'that specifies the multiplicities at each position of \code{knots}
#'@param native If \code{TRUE} native code is used.
#'@return A support that is adequate for clamped spline evaluation
#'@author Máximo Sánchez-Aragón
#'@export
bspline_support <- function( degree,
                             knots=c(0:5),
                             multiplicities=c(degree+1,
                                              rep(1,length(knots)-2),
                                              degree+1), 
                             native=FALSE ){

  
  if(length(knots) != length(multiplicities))
    stop("the number of knots must equal the number of multiplicities")
  
  if(!native){

    for(i in 1:(length(knots)-1) )
      if(knots[i] >= knots[i+1])
        stop("knots must be unique and increasing")
  
    r <- c()
    for( i in 1:length(knots) )
      r <- c(r, rep(knots[i], times=multiplicities[i]) )
  
    return(r)
  
  }else{
  
    support <- rep(0, times=sum(multiplicities))
    errcode <- 0
    result <- .C("bspline_support", as.integer(degree), as.numeric(knots), 
                as.integer(length(knots)), 
                as.integer(multiplicities),
                as.numeric(support),
                as.numeric(errcode)
            )
    
    if(result[[6]] > 0)
      stop ("knots must be unique and increasing")
    else
      return(result[[5]])
    
  }
  
}


