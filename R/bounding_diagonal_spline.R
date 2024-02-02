#'@section Auxiliary functions:
#'@title Makes a straight diagonal spline curve
#'@description \code{bounding_diagonal_spline} makes a diagonal straight spline from the data point bounding box.
#'This function returns the control points of the spline, as a matrix of row-vectors. It has been used 
#'as onset for spline fitting routines in some examples, but it should not be standarized - each specific problem
#'requires careful consideration.
#'@param x The data points coordinates, as a matrix of row-vectors
#'@param degree The degree of the spline
#'@param f A coefficient that determines how interior are the points of the computed diagonal. Higher values
#'correspond with more interior points
#'@return A matrix of row-vectors representing the control points of a spline curve.
#'@author Máximo Sánchez-Aragón
#'@export
bounding_diagonal_spline <- function(x, degree, f=1.2){

  # gets the bounding diagonal as two points, the minimum and the maximum
  # in row vector notation
  bx <- c()
  for(j in 1:ncol(x))    
    bx <- cbind(bx, c( min(x[,j]), max(x[,j]) ))
  
  # now there are <degree> - 1 points missing, for which the full length
  # has to be divided <degree> times  
  uvec <- (bx[2,] - bx[1,]) / sqrt(sum((bx[2,] - bx[1,])^2))
  bx[1,] <- bx[1,] + uvec*f
  bx[2,] <- bx[2,] - uvec*f
  step <- sqrt(sum((bx[2,] - bx[1,])^2)) / degree
  
  cont_p <- c()
  for(j in 0:degree){
    nc <- bx[1,] + uvec*step*j
    cont_p <- rbind(cont_p, nc)
  }
  
  row.names(cont_p) <- NULL
  
  cont_p
  
}
