#'@title Recalculating control points for the spline derivative.
#'@description \code{recompute_control_points} computes the control points
#'for the derivative of a spline.
#'@details The derivative of a spline built on k-degree polynomials
#' can be expressed as the same spline built on (k-1) polynomials and
#' a different set of control points. The new control points (calculated by this function)
#' depend on the knot sequence and the current control points.
#'@param v Control points of the current spline, as a row-vector matrix.
#'@param U The support, as a vector; it must be equal in length than the number of control points.
#'@param native If TRUE native code is used (default).
#'@param degree The degree of the polynomial for the current spline.
#'@return The control points for the derivative, as a matrix of row-vector coordinates.
#'@author Máximo Sánchez-Aragón
#'@export
recompute_control_points <- function(v, U, degree=1, native=TRUE){

  if(!native){
    
    a <- matrix()
    for(i in 1:nrow(v)-1){
      if(i <= 1)
        a <- ( degree / ( U[i+degree+1] - U[i+1] ) ) * (v[i+1,] - v[i,])
      else
        a <- rbind(a, ( degree / ( U[i+degree+1] - U[i+1] ) ) * (v[i+1,] - v[i,]))
    }
    row.names(a) <- NULL
    return(a)
    
  }else{
    
   result <- .C("recompute_control_points", as.numeric( v ),
                                            as.integer(nrow(v)),
                                      as.integer(ncol(v)),
                                      as.numeric(U),
                                      as.integer(degree)  )
   R <- t(matrix(result[[1]], nrow=ncol(v), ncol=nrow(v)))[1:(nrow(v)-1), ]
   return(R)
  }
}




