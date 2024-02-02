#'@title Fitting of a spline curve with parameter correction
#'@description \code{bspline_fit} performs a spline fit to a set of data points by least squares. The only method provided for now is
#''point distance' (see method), i.e. it minimizes euclidean distances of the data points to the curve.
#'@param x The data points, as a matrix of row-vector coordinates
#'@param v Initial set of control points used, as a matrix of row-vector coordinates.
#'@param sup The support of knots. By default the support is set for 'clamped splines'
#'@param method (experimental) The method used. Currently the method is 'point distance', i.e. the routine ensures
#'that the squared euclidean distances of the data points to its footpoints is minimal. NOTE: changing this parameter does
#'not have any effect. It is there only as a reminder for further development.
#'@param tolerance The mean of squared differences (MSE) between control points of contiguous rounds below which
#'the routine exits. Lower values provide higher precision by performing more rounds with the risk of further numerical instability.
#'@param fp.tolerance The tolerance in the parameter correction step (see \code{bspline_footpoint})
#'@param reg.factor The Tikhonov regularization factor used during the resolution step (see \code{bspline_solve})
#'@param verb.progress Shows the mean squared deviation of the control points in each iteration.
#'@param ... other parameters passed to \code{bspline_footpoint}
#'@return A matrix of row-vector representing the control points of the fitted spline curve.
#'@example examples/bspline_fit_example1.R
#'@example examples/bspline_fit_example2.R
#'@author Máximo Sánchez-Aragón
#'@seealso \code{bspline_footpoint}
bspline_fit <- function(x, 
                        degree=3, 
                        v=bounding_diagonal_spline(x,degree,f=1.2),
                        sup=bspline_support(degree, knots=seq(from=0, to=1, length.out=nrow(v)-degree + 1)),
                        method="PD", 
                        tolerance=1e-9, 
                        fp.tolerance=1e-2, reg.factor=0, verb.progress=FALSE,
                        ...){

  if(ncol(x) != ncol(v))
    stop("[Max] the data point matrix and the control point matrix must have the same
          number of columns (i.e. equal dimensions)")
  
  nv <- v
  ftp <- bspline_footpoint(x, nv, degree = degree, 
                           tolerance=fp.tolerance, ...)
  
  tol <- Inf
  n <- nrow(nv)
  current_it <- 0
  while(tol > tolerance) {
        
    current_it <- current_it + 1
    nv_prev <- nv #makes a copy of the control points
    
    #makes a first approximation
    nv <- bspline_solve(x, ftp, nv, degree=degree, 
                        reg.factor=reg.factor, sup=sup)
    
    ftp <- bspline_footpoint(x, nv, u=ftp, degree = degree, 
                             tolerance=fp.tolerance, ...)    
    
    #tolerance is computed as the shift between computed control points
    #and the previous control points
    sm <- (nv_prev - nv)^2          
    tol <- sum(sm) / n    
    
    if(verb.progress)
    print( paste("i:",current_it,"mse:",tol) )  
    
  }
  
  nv
  
}
 
