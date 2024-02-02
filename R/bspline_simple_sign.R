#' @section Auxiliary functions:
#' @title Assigning sign to data points by simple differential geometry operation
#' @description \code{bspline_simple_sign} computes the sign of a point projection as the sign of the z component 
#' of the vector product between the spline tanget vector and the projection vector at the footpoint. The sign is normalized by
#' an anterior-posterior axis vector (\code{ap_ref}), provided by the user.
#' 
#' This is only an auxiliary function for specific problems of modelling of the differentation wave in
#'  \emph{Drosophila} eye imaginal discs.
#'  
#' The current implementation only works with 3-dimensional vectors.
#' @param points Row-vector matrix of positions to be projected.
#' @param c.points Row-vector matrix of control points.
#' @param ap_ref Anterior-posterior axis vector.
#' @param proj.points The row-vector matrix of footpoints. If NULL they are calculated.
#' @param deg The spline polynomial degree.
#' @param fp.tolerance The tolerance applied during footpoint computation.
#' @author Máximo Sánchez-Aragón
#' @export
bspline_simple_sign <- function(points, c.points, ap_ref ,proj.points = NULL, deg=3, 
                                fp.tolerance=1e-3){
  
  # establish the orientation criterium 
  dv_m <- interpolate_bspline(c(0,1), c.points, degree=deg)
  dv_ref <- as.numeric(dv_m[2,] - dv_m[1,])
  mt <- rbind( 
                c(1,1,1), 
                ap_ref,
                dv_ref
             )
  mtsgn <- c(
        #   m[2,2]*m[3,3] - m[3,2]*m[2,3],
        # - m[2,1]*m[3,3] + m[2,3]*m[3,1],
        mt[2,1]*mt[3,2] - mt[3,1]*mt[2,2]
      ) >= 0
  
  # generate projection points upon request
  if( is.null(proj.points) ){
    ud <- bspline_footpoint(points, c.points, degree = deg, verbose=TRUE, tolerance = fp.tolerance)
    proj.points <- interpolate_bspline(ud, c.points, degree = deg, derivate=0)
  }
  
  # computes the first derivative 
  proj.points.d <- interpolate_bspline(ud, c.points, degree = deg, derivate=1)
  
  #ud   : projection parameter values
  #PM0  : matrix of cell positions (as row vector matrix)
  #PM1  : matrix of cell projection positions (as row vector matrix)
  #PM1d : matrix of cell projection derivatives (as row vector matrix)
  
  PM0 <- points
  PM1 <- proj.points
  PM1d <- proj.points.d
  
  
  ors <- rep(FALSE, times=nrow(PM0))
  for(i in 1:nrow(PM0)){
    p1 <- as.numeric( PM0[i,] )
    p2 <- as.numeric( PM1[i,] )
    p <- p2-p1;
    m <- rbind(c(1,1,1),PM1d[i,],p)
    p.cross <- c(
      #   m[2,2]*m[3,3] - m[3,2]*m[2,3],
      # - m[2,1]*m[3,3] + m[2,3]*m[3,1],
      m[2,1]*m[3,2] - m[3,1]*m[2,2]
    )
    if(mtsgn)
      ors[i] = p.cross < 0
    else
      ors[i] = p.cross >= 0
    
    #print(paste("ors test >>> ",ors[i]) )
  }
  rm(p1,p2,m,p, p.cross,i)
  
  result <- list()
  result[[1]] <- NULL #rep(0.5, times=length(ors))  # foo 
  result[[2]] <- ors
  
  #attr(result[[1]], "description") <- "foo, at the position of the u parameters (allegedly) - check parity_sign function for more info"
  #attr(result[[2]], "description") <- "the result of the xor operation (parity sign)"
  
  result
}

