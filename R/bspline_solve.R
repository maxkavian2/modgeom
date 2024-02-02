#'@title Finding the control points for a fitted spline curve.
#'@description \code{bspline_solve} finds the control points whose spline curve minimizes the squared distance to specified data points by 
#'QR decomposition.
#'@param x The data points in the form of a matrix of row-vectors.
#'@param u The sequence of parameter values corresponding to the best footpoints guesses.
#'@param v The control points of the spline curve.
#'@param degree The spline polynomial degree.
#'@param sup The support of knots. By default it is adjusted to a clamped spline.
#'@param reg.factor  It introduces a regularization factor (Tikhonov's) in order to limit the 
#'space of solutions. The Tikhonov matrix is made out of the product of this regularization factor and the identity matrix.
#'@param native if \code{TRUE} it executes the native library code.
#'@return A matrix of row-vectors representing the control points.
#'@seealso \code{bspline_fit}
#'@example examples/bspline_solve_example1.R
#'@author Máximo Sánchez-Aragón.
#'@export
bspline_solve <- function(x, u, v, degree=3,                        
                          sup=bspline_support(degree, knots=seq(from=0, to=1, length.out=nrow(v)-degree + 1)),
                          reg.factor=1, native=T){
  
  
  # A <- matrix(c(1,2,3,4,5,6), nrow=3, ncol=2)
  # B <- matrix(c(1,2,3,4,5,6), nrow=3, ncol=2)
  # 
  # print(A);
  # print(B);
  # 
  # print(crossprod(A,B));
  
  if(!native){
    
    N <- bspline_eval(u, U=sup, degree=degree)  # corrected line
    Nprod <- crossprod(N,N)
    Nprod <- Nprod + diag(nrow(Nprod))*(reg.factor^2)
    R <- qr.solve(Nprod, crossprod(N,x))
    return(R)
    
  }else{
    
    un <- length(sup)-degree-1
    B <- rep(0, times=un*length(u))
    R <- rep(0, times= (length(sup) - degree - 1)*ncol(x))
    
    result <- .C("bspline_solve",
                as.numeric(u),
                as.integer(length(u)),
                as.integer(degree),
                as.numeric(sup),
                as.integer(length(sup)),
                as.numeric(x),
                as.integer(ncol(x)),
                as.numeric(R), 
                as.numeric(reg.factor)
                )
    
    return(matrix(result[[8]], nrow=(length(sup) - degree - 1), ncol=ncol(x)) )
    
  }
}
