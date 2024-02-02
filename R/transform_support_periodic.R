#'@title Turning knot support into periodic.
#'@description \code{transform_support_periodic} turns a clamped support into 
#' a periodic support.
#'@param degree The degree of spline polynomials that use this support.
#'@param U A clamped support, created by the function \code{bspline_support}.
#'@param native If TRUE native code is used (default).
#'@return A periodic support vector.
#'@author Máximo Sánchez-Aragón
#'@export
transform_support_periodic <- function(degree = 3, 
                                       U = bspline_support(degree, knots=seq(from=0,to=1,
                                                                             length.out=10) ), 
                                       native = FALSE) {
  sup_per <- U
  if(degree <= 0)
    return(sup_per)
  
  if( !native ){
    
    n <- length(U)
    lint <- U[(n-degree):(n-2*degree+1)] - U[(n-degree-1):(n-2*degree)]
    
    lint <- lint[length(lint):1]
    lint2 <- lint
    
    for(i in 1:length(lint)) lint2[i] <- sum(lint[i:length(lint)])
    lint <- lint2
    
    sup_per[1:degree] <- rep(U[degree], times = degree) - lint
    
    lint <- U[(2+degree):(2*degree+1)] - U[(1+degree):(2*degree)]
    lint2 <- lint
    
    for(i in 1:length(lint)) lint2[i] <- sum(lint[1:i])
    lint <- lint2
    
    sup_per[(n-degree+1):n] <- rep(U[n-degree], times=degree) + lint
    
    #print(sup_per) # TEST LINE
    
    return(sup_per)
    
  }else{
    
    result <- .C("transform_support_periodic", 
                 as.integer(degree),
                 as.numeric(sup_per),
                 as.integer(length(sup_per)) )
    
    return(result[[2]])
    
  }
}