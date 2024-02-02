#'@title B-spline evaluation
#'@description \code{bspline} evaluates all B-splines for the single value \code{x}
#'@details This function evaluates all B-splines at a single value of the parameter \code{x}
#'@param x The parameter value (single value).
#'@param k The index of the containing interval for x. If it is
#'equal to -1 the interval is computed.
#'@param U The knot support, as a vector.
#'@param native If \code{TRUE} Native code will be used instead.
#'@return A vector of B-splines evaluations
#'@author Máximo Sánchez-Aragón
#'@export
bspline <- function(
  x, k=-1,
  degree=3,
  U=bspline_support(degree),
  native=TRUE
){

  if( x < U[1] || x > U[length(U)] )
    stop("interpolating point out of reach!");

  # makes some consistency tests
  if(degree < 0) {
    degree = 0;
    warning("degree below zero has been adjusted to zero")
  }

  if(!native){
   
     # makes a default spline if the last point has been reached
    # i.e. the user tries to interpolate the last point
    if(x == U[length(U)]){
      warning("last point reached! setting default b-spline")
      B<-c(rep(0, times=length(U)-degree-1))
      B[length(U)-degree-1]=1
      return(B)
    }
  
    # finds the interval index in case it is not specified
    if(k==-1)
      k = bspline_find_interval_index(x, U)
  
    #1.) builds the 0-degree B-spline
    B<-c( rep(0, times=length(U)-degree-1) )
    B[k] = 1
  
    #2.) evaluates using the recursion theorem for B-splines (de Boor)
    p = 0               # the current degree
    lenB <- length(B)   # the length of the spline
  
    while(p < degree ){
  
      p <- p+1
  
      # it evaluates the middle b-splines. Both terms may be not-null
      for(i in (k-p):k) {
        t1 <- 0
        t2 <- 0
  
        # test the first term of the recursion
        if( U[i+p] > U[i] )
          t1 <- ( (x - U[i]) / ( U[i+p] - U[i] ) ) * B[i]
        
        # test the second term of the recursion
        if( U[i+p+1] > U[i+1] ) # idem
          t2 <- ( ( U[i+p+1] - x ) / ( U[i+p+1] - U[i+1] ) ) * B[ (i %% lenB)+1 ]
  
        B[i] <- t1 + t2
      }
    }
    
    return(B)
  
  }else {
    # here goes the native code (check bspline_core.cpp file)
    errcode <- 0;
    B <- c( rep(0, times=length(U)-degree-1) ) # makes the null spline evaluation
    result<- .C("bspline", as.numeric(x),
                    as.integer(k),
                    as.integer(degree),
                    as.numeric(U),
                    as.integer(length(U)),
                    as.numeric(B),
                    as.integer(length(B)),
                    as.integer(errcode) )
    
    if(result[[8]] == 2)
      warning("last point reached! setting default b-spline")
    
    
    return(result[[6]])
    
  }
  
  
}
