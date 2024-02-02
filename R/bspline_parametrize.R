#'@title Finding parameter values for a set of interpolation points
#'@description \code{bspline_parametrize} finds the parameter values for the spatial points
#'provided, i.e. points belonging to the spline curve, by one of several methods.
#'@param interp.points The interpolation points, as a matrix of row-vector positions
#'@param method The interpolation method. Currently there are only
#'three options: 'uniform', 'chordal' or 'centripetal' (default)
#'@param limit The upper bound for the knot sequence, which is set to 1 by default.
#'Lower bound is set to zero in all cases.
#'@return paramater values for the provided curve points
#'@author Máximo Sánchez-Aragón
#'@export
bspline_parametrize <- function ( interp.points ,
                                  method="centripetal", limit = 1){


  r <- seq( from=0, to=limit, length.out=nrow(interp.points) )

  if( method == "centripetal"){

    s <- 0
    for( i in 2:nrow(interp.points) ){
      delta <- sqrt(sqrt(sum((interp.points[i,] - interp.points[i-1,])^2)))
      r[i] <- r[i-1] + delta * limit
      s <- s + delta
    }
    r <- r / s

  }else if( method == "chordal" ){

    s <- 0
    for( i in 2:nrow(interp.points) ){
      delta <- sqrt(sum((interp.points[i,] - interp.points[i-1,])^2))
      r[i] <- r[i-1] + delta * limit
      s <- s + delta
    }
    r <- r / s

  }else if( method == "uniform" ){

  }else
    stop("[Max] unknown method selected.")

  r
}
