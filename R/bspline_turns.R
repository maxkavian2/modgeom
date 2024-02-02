#'@title Turns of a spline curve relative to a reference vector
#'@description This function returns a set of parameter values corresponding with the 
#'points in a spline curve where a change of pi/2 degrees in the tangent relative to a 
#'reference vector occurs. 
#'
#'This function is deprecated after the functions \link{bspline_uz_sign}, 
#'\link{bspline_parity_sign} and \link{bspline_parity_turns}. 
#'@param c_point The current control points, expressed as a matrix of row vectors.
#'@param degree The maximum polynomial degree of the spline curve.
#'@param initial.scan.points Number of initial points that are scanned.
#'@param knots The sequence of knots of the spline curve support.
#'@param verbose If \code{TRUE} it shows information while computing.
#'@param max.iterations  the maximum number of allowed iterations allowed.
#' If the value is exceeded, the algorithm interprets this as a non-orthogonal
#' minimum, and discards it with a warning.
#'@param ref_vector The reference vector used. If \code{NULL} the tangent vector at 
#'u = 0 will be used.
#'@param pick.sides Picks maximum and minimum values only from minimum candidates 
#'for u.
#'@param ... Other parameters passed to interpolate_bspline.
#'@return A \code{data.frame} that includes the parameter of a turn, along 
#'other data such as the actual dot product, the second derivative and the 
#'number of NRA iterations.
#'@author Máximo Sánchez-Aragón
#'@export
bspline_turns     <- function(c_point, 
                              degree = 3,
                              tolerance=1e-9,
                              initial.scan.points=nrow(c_point)*8,
                              knots = seq(from = 0, to = 1, length.out=nrow(c_point) - degree + 1),
                              verbose=FALSE,
                              max.iterations = 100,
                              ref_vector = NULL,
                              pick.sides = FALSE,
                              ...){
  
  "
  
  finds the turns in the curve relative to the initial vector
  i.e. the number of rows of the x matrix). One is provided by default

  c_point : the current control points, expressed as a matrix of row vectors
  initial.scan.points : number of initial points that are scanned when
  u is not null.
  degree : the degree of the bspline curve
  tolerance : the threshold, as absolute dot product value, below which
              the routine ends.
  max.iterations : the maximum number of iterations allowed. If the value
              is exceeded, the algorithm interprets this as a non-orthogonal
              minimum, and discards it with a warning,
  verbose : shows information while computing
  ref_vector : the reference vector used. If null the vector at u = 0 
               will be used
  eliminate.redundant : eliminates redundant solutions
  pick.sides : picks maximum and minimum values only from minimum candidates for u
  ... : other parameters passed to interpolate_bspline

  
  returned value: a data.frame that includes the parameter of a turn, along 
  other data such as the actual dot product, the second derivative and the 
  number of Newton-Raphson algorithm iterations.
  
  "
  
  #defines the user knots and the support (clamped)
  #knots <- seq(from = 0, to = 1, length.out=nrow(c_point) - degree + 1)
    
  if(verbose){
    print("finding turns ...")    
  }
  
  # INITIALIZATION ...
  #u <- 0 # the tested point
  u.found <- c()
  der1 <- c()
  der2 <- c()
  iterations <- c()

  mink <- min(knots)
  maxk <- max(knots)
      
  # 1.-) finds the minimum reference for the next iteration  
  if(is.null(ref_vector))
  xd_ref <- interpolate_bspline(0,c_point, degree=degree, derivate=1, 
                                unitary=F, ...)    
  else
  xd_ref <- ref_vector
  
  v <- seq(from=0, to=1, length.out=initial.scan.points) 
  pnt_derivates_unit <- interpolate_bspline(v,c_point, degree=degree, derivate=1, 
                                            unitary=F, ...) 
  CP <- as.vector(pnt_derivates_unit %*% t(xd_ref) )  
  #u <- v[ which( abs(CP) == min(abs(CP)) ) ]  # NOTE finds one minimum only
  
  # finds all local minima
  u.min <- c()
  # CP <- CP^2
  for(i in 1:(length(CP)-1) )
    if( (0<CP[i+1] & 0>=CP[i]) | (0>CP[i+1] & 0<=CP[i]) )
      u.min <- c(u.min, v[i])    
  
  
  if(length(u.min) <= 0){
    warning("[Max] no suggested minima found. returning an empty value")    
    return(NULL)
  }  
  
  if(pick.sides){
    u.min2 <- c(min(u.min),max(u.min))
    u.min <- u.min2
    rm(u.min2)
  }
    
  
  if(verbose)
    print(paste(length(u.min)," minima suggested",sep=""))
  
  # iterates over all minima
  for(u in u.min){
  
  # RECURRENCE ...  
  # 2.-) performs Newton-Gauss minimization on each point
  tol <- Inf 
  count <- 0 
  b <- TRUE
  while(tol > tolerance){ # this loop must test the iteration quality
    
    # evaluates the point     
    xd   <- interpolate_bspline(u, c_point, degree = degree, derivate = 1, 
                              knots = knots, ...)
    xdd  <- interpolate_bspline(u, c_point, degree = degree, derivate = 2, 
                               knots = knots, ...)
    #xddd <- interpolate_bspline(u, c_point, degree = degree, derivate = 3, 
    #                            knots = knots, ...)       
    
    # multiplies by the reference
    xd   <- xd_ref %*% t(xd)    
    xdd  <- xd_ref %*% t(xdd)
    #xddd <- xd_ref %*% t(xddd)
    
    # newtonian expression
    fd  <- as.numeric( xd )  # as.numeric( xd  * xdd )
    fdd <- as.numeric( xdd ) # as.numeric( xd  * xddd + xdd * xdd ) 
    
    inc_u   <- -fd / fdd
        
    count <- count+1
    tol <- xd * xd #abs(xd)
    print(as.numeric(tol))
    
    u <- u + inc_u  
                 
    if(count > max.iterations || u < mink || u > maxk){
      warning("[Max] non-orthogonal minimum found. discarding ... ")
      b <- FALSE
      break
    }
              
  }
  
  # adds the found minimum
  if(b){
    
    u.found <- c(u.found, u)
    iterations <- c(iterations, count)
    der1 <- c(der1, xd)
    der2 <- c(der2, xdd)
      
    if(verbose)
      print(paste("minima found ...","it.:",count," sec.der.:",xdd))      
    
  }  
  # ENDS of RECURRENCE ...
  
  }
  
  
  result <- data.frame(u.found, u0 = u.min, iterations, der1, der2)
    
  # returns the found minima
  result
    
}


