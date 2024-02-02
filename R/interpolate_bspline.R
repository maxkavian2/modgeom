#'@title Spline interpolation
#'@description \code{interpolate_bspline} computes the curve coordinates of
#' a numeric vector of parameter values or its derivatives.
#'@details This function evaluates a spline and finds the coordinates of the
#'provided parametric values (i.e. interpolation) or its derivatives. Splines
#'are clamped by default between the lowest and highest knot values. The current implementation selects 
#'the lowest and the highest knot values to build a new knot sequence of evenly distributed values that 
#'fulfils the proper number of knots, so deleting the original knot sequence and throwing a warning.
#'@param x The parameter values, as a numeric vector.
#'@param v Control points of the spline, as a row-vector matrix of coordinates.
#'@param degree The degree of the polynomials that build the spline.
#'@param derivate The derivative order to be evaluated. With 0 it just evaluates the
#'spline at the parameter values (i.e. it does not calculate the derivative)
#'@param knots The sequence of knots used for the support, as a numeric vector.
#'This vector must have a proper length (i.e. number of control points minus degree plus one).
#'Knots are always ordered in an ascending fashion. (see note below about periodic splines)
#'@param is.periodic It specifies if the spline is periodic (\code{TRUE}); otherwise the spline
#'would be clamped (default).
#'@param unitary If \code{TRUE} unitary vectors are returned.
#'@param native If \code{TRUE} native code is used (default),
#'@param ... Parameters passed to the function \code{\link{bspline}}. NOTE that
#'that function is inner, hence the parameters that can be passed should be described here.
#'@return A matrix containing a set of row-vectors, representing curve positions at the
#'interpolated parametric values or its derivatives
#'@author Máximo Sánchez-Aragón
#'@example examples/interpolate_bspline_example1.R
#'@example examples/interpolate_bspline_example2.R
#'@example examples/interpolate_bspline_example3.R
#'@export
interpolate_bspline<-function(x, v, degree = 3,
                              derivate = 0,
                              knots = seq(from = 0, to = 1,
                                          length.out=nrow(v) - degree + 1),
                              is.periodic=FALSE,
                              unitary = FALSE, native = TRUE, ...){

  if (!native){
    
    knots <- knots[order(knots, decreasing=F)] # order the knots
    
    # ensures that the equivalence between knots limits
    # and the control points
    if(is.periodic){
      x[x == max(knots)] <- min(knots)
      vc <- rbind(v,v[1:degree,])
      knots <- seq(from=knots[1], to=knots[length(knots)], length.out=nrow(vc)-degree+1)
      warning("[Max] knots have been evenly distributed over the interval to fit the periodic spline operation. The original knot sequence is lost.")
    }else
      vc <- v
    
    
    #print(vc)
    # controls that derivation is feasible
    if(derivate < 0 || derivate > degree)
        stop("[Max] the derivate cannot be negative or higher than the degree")
    
    if(length(knots) <= 1)
        stop("[Max] the minimum length for the knot sequence is 2. Either increase the number of control points or decrease the current degree.");
    
    if(nrow(vc) <= 1)
        stop("[Max] the minimum number of control points is 2, don't be silly!");
  
  
    #copies the control points
    vt <- vc
  
    # original support
    if(derivate > 0)
      for(i in 1:derivate){
        sup_orig <- bspline_support(degree-i+1, knots = knots, native=TRUE)
  
        if(is.periodic)
          sup_orig <- transform_support_periodic(degree=degree-i+1,
                                                 U=sup_orig)
  
        vt <- recompute_control_points(vt, sup_orig, degree=degree-i+1)
      }
  
    # calcula el soporte de acuerdo con la derivada (new knot sequence)
    sup_orig <- bspline_support(degree - derivate, knots = knots, native=TRUE)
  
    if(is.periodic)
      sup_orig <- transform_support_periodic(degree=degree - derivate,
                                             U = sup_orig)
  
    B <- bspline_eval(x, degree = degree - derivate, U = sup_orig, ...
                      )
  
    R <- B %*% vt
  
    # normalizes the derivative vectors ...
    if(derivate > 0 && unitary)
      for(j in 1:nrow(R)   )
        R[j,] <- R[j,] / sqrt( sum(R[j,]^2) )
  
    return(R)
  
  }else{

    #if(is.periodic)
    #  stop("[Max] periodic spline implementation is not yet available for native code.");
    vc <- v;
    if(unitary)
      warning("[Max] unitary vectors will not be returned in native mode")
    
    errcode <- 0
    R <- rep(0, times=length(x)*ncol(vc))
    result  <- .C("interpolate_bspline",as.numeric(x),
                                        as.integer(length(x)),
                  
                                        as.numeric(vc),
                                        as.integer(nrow(vc)),
                                        as.integer(ncol(vc)),

                                        as.integer(degree),
                  
                                        as.numeric(knots),
                                        as.integer(length(knots)),
                  
                                        as.integer(derivate),
                                        as.numeric(R),
                  
                                        as.integer(is.periodic),
                                        errcode
                                        )
    

    if(result[[12]] == 1)#  & !is.periodic) # periodic 
      stop ("[Max] knots must be unique and increasing.")
    
    R <- matrix( result[[10]], nrow=length(x), ncol=ncol(vc) )
    #print(R)
    return(R)
    
  }

}
