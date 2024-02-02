#' @title Projection angles computation.
#' @description \code{projection_angles} computes the angles of data point projections onto the assigned footpoints in a
#' spline curve. These footpoints are not necessarilly optimized.
#' @param ap_u Set of parameter values defining the points the curve
#' @param d_point A row-vector matrix with the data points to be projected.
#' @param c_point Control points of the spline curve.
#' @param degree The degree of the spline curve.
#' @param unit The angle unit, i.e. "radian" or "degree"
#' @param ... Other arguments passed to \code{interpolate_bspline}
#' @return a vector with the computed angles
#' @author Máximo Sánchez-Aragón
#' @export
 
projection_angles <- function(ap_u, d_point, c_point, degree = 3,
                              unit="degree", ...){
  Rd <- interpolate_bspline(ap_u, c_point, degree = degree, ...)
  Rd_der <- interpolate_bspline(ap_u, c_point, degree = degree, derivate = 1, ...)
  
  r <- c()
  
  for(i in 1:nrow(d_point)){
    a <- (d_point[i,]- Rd[i,]) * Rd_der[i,]
    m1 <- sqrt( sum( (d_point[i,]- Rd[i,]) * (d_point[i,]- Rd[i,]) ) )
    m2 <- sqrt( sum( Rd_der[i,] * Rd_der[i,] ) )
    
    if(unit == "degree")
    r[i] <- acos(sum(a)/(m1 * m2))*360/(2*pi)
    else
    r[i] <- acos(sum(a)/(m1 * m2))
  }
  
  r
    
}
