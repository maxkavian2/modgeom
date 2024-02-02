#'@title Spline evaluation
#'@description \code{bspline_eval} evaluates at once the parameter values of vector \code{x}
#'@details This function evaluates all B-splines at the values of the vector \code{x}.
#'It is important to remark that no warning will be thrown when \code{native} equals \code{TRUE}
#'for parametric values (\code{x}) that have reached the support highest value.
#'@param x The parameter values, as a numeric vector.
#'@param degree The degree of the polynomials that define the spline
#'@param U The support, as a numeric vector.
#'@param native If \code{TRUE} a native implementation is used instead
#'@return A matrix containing all B-splines evaluations. Evaluations for the same
#'parameter value are contained in individual rows. Columns represent different evaluations
#'of the same B-spline at different parameter values.
#'@author Máximo Sánchez-Aragón
#'@example examples/bspline_eval_example2.R
#'@export
bspline_eval <- function( x, 
                          degree=3,
                          U=bspline_support(degree),
                          native=TRUE ){

  B <- c()
  if(!native){
    
    for(xi in x)
      B <- rbind( B, bspline(xi, k=-1, degree=degree, U=U, native=native ) )
    return(B)
    
  }else{
    
    xn <- length(U)-degree-1
    B <- rep(0, times=xn*length(x))
    
    result  <- .C("bspline_eval", as.numeric(x),
                                  as.integer(degree),
                                  as.numeric(U),
                                  as.integer(length(U)),
                                  as.numeric(B),
                                  as.integer(xn),
                                  as.integer(length(x))
                               )
    
    B <- matrix(result[[5]], ncol=xn, nrow=length(x))
    return(B)
    
  }
}
