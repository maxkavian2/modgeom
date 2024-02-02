#' @title Finding the length of a spline curve by a parameter subinterval
#' @description It computes the length of a spline curve by the extended sum of the distances of a sequence
#' of interpolated points in the interval of parameter values specified by the user. The sequence of
#' interpolated points is chosen from an evenly distributed set of parameter values in the subinterval
#' that which is intended to measure.
#' @param c_point the row-vector matrix of control points
#' @param degree the spline polynomial degree
#' @param resolution The number of chunks used for length computation. This number is applied to a 
#' parameter domain [0,1] and it is rescaled accordingly with the arguments of \code{begin} and \code{end}
#' @param begin the lower bound of the parameter value interval from which the length of the spline curve is computed
#' @param end the upper bound of the parameter value interval from which the length of the spline curve is computed
#' @author Máximo Sánchez-Aragón
#' @export
bspline_length    <- function(c_point, degree=3, resolution=100, 
                              begin=0, end=1, ...)
         {
                      
              "
                computes the bspline length between the values of the parameters 
                      specified.
                    
                c_point : control points, as row vector points
                degree  : the degree of the bspline
                resolution : the number of chunks used for the length computation
                begin : the staring parameter value
                end   : the ending parameter value
                ...   : other parameters pased to interpolate_bspline
                  
              "
              knots = seq(from = 0, to = 1, length.out=nrow(c_point) - degree + 1)
              u <- seq(from = begin, to = end, 
                       length.out=(resolution * (end-begin) +1)  )
              
              points <- interpolate_bspline(u, c_point, degree=degree, knots = knots, ...)
              rm(u)
              curve.length <- 0
              for(i in 1:nrow(points)-1) {
                a <- points[i+1,] - points[i,]                           
                distance <- sqrt(sum(a * a))
                curve.length <- curve.length + distance
              }
              rm(a, distance,i)
              curve.length
           }
