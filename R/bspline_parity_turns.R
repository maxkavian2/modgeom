#' @title Finding turns for a set of projections onto a spline curve
#' @description  This function is intended to provide a point specific turn evaluation in the
#' spline using tensors, where we follow the definition of turn given in the PhD thesis, i.e. it is based on
#' the osculating plane method (OP). Projections are generated between
#' data points external to the curve and their footpoints, the latter defined as the closer points in the
#' curve that belong to 
#' each of the data points (see \code{\link{bspline_footpoint}}) in case the user does not specify otherwise.
#' 
#' This function specifically depends on tensorA package. The current evaluation only works
#' for three dimensional space, i.e. all row-vector matrices defining positions must have three columns.
#' 
#' @param c_points Row-vector matrix containing the positions of the three dimensional control points
#' @param points Row-vector matrix containing the positions of data points
#' @param footpoints  The precalculated footpoints (as a vector of parameter values) for each element in \code{points}, 
#' If NULL the algorithm tries to calculate them with parameters fp.*. (see \code{\link{bspline_footpoint}}).
#' @param degree The spline polynomial degree.
#' @param tolerance The threshold, as absolute dot product value, below which 
#' the routine ends for each point independently. It operates for each element in \code{points} independently
#' @param max.iterations Maximum number of allowed iterations before the algorithm exits if \code{tolerance} has not been
#' reached. It operates for each element in \code{points} independently.
#' @param initial.scan.points Initial number of evenly distributed starting positions from which the turns
#' are approximated afterwards.
#' @param knots The ordered knot sequence of the spline, not the support. By default the calculated spline is clamped
#' @param scan.interval (experimental) It specifies the bounds of the subdomain of the spline parameter where the turns 
#' will be searched. The argument must contain two elements.
#' @param fp.tolerance The tolerance passed to \code{\link{bspline_footpoint}}
#' @param fp.initial.scan.points The initial number of evenly distributed footpoints passed to \code{\link{bspline_footpoint}}
#' @param fp.max.iterations The maximum number of iterations passed to \code{\link{bspline_footpoint}}
#' @param fp.algorithm The algorithm mode passed to \code{\link{bspline_footpoint}}
#' @param fp.verbose If \code{TRUE} displays progress information.
#' @param ... Argument passed to \code{\link{interpolate_bspline}}
#' @return A \code{list} of \code{data.frame}s with the parameter values corresponding to the turns and additional information. 
#' Each \code{data.frame} also contains information about the total number of cycles and the final tolerance value, the former 
#' being an indication of the turn precision.
#' @example examples/bspline_parity_turns_example1.R
#' @author Máximo Sánchez-Aragón
#' @export
bspline_parity_turns <- function(  
                                c_points,
                                points,
                                footpoints = NULL,
                                degree = 3,
                                tolerance=1e-9,
                                max.iterations = 100,
                                initial.scan.points=nrow(c_points)*8,
                                knots = seq(from = 0, to = 1, length.out=nrow(c_points) - degree + 1),
                                scan.interval = c(min(knots),max(knots)),
                                fp.tolerance = 1e-12,
                                fp.max.iterations = 800,
                                fp.initial.scan.points=nrow(c_points)*8,
                                fp.algorithm = "directional",
                                fp.verbose = FALSE,
                                ...
){
  # here some SANITY check: dimensionality of the control points
  if(ncol(c_points) != 3)
    stop("[Max says] the control points must have three dimensions")
  
  # here some SANITY check: dimensionality of the test points
  if(ncol(points) != 3)
    stop("[Max says] the test points must have three dimensions")
  
  # here some SANITY check: the number of elements of the scanning interval is strictly 2
  if(length(scan.interval) != 2)
    stop("[Max says] the number of elements in the scanning interval must be 2")
  
  # here some SANITY check: scan.interval must preserve the order
  if(scan.interval[1] >= scan.interval[2])
    stop("[Max says] the elements in the scanning interval must preserve the order (lower first, higher second)")
  
  
  # here some SANITY check: scan.interval must be contained in the knot internal
  knot.int <- c(min(knots), max(knots))
  scan.int <- c(min(scan.interval), max(scan.interval))
  if( scan.int[1] < knot.int[1] | scan.int[2] > knot.int[2] )
    stop("[Max says] the parameter scanning interval should be contained in the knot interval!")
  
  rm(knot.int, scan.int)
  
  
  # finds the projection points if null
  if(is.null(footpoints)){
    print("solving footpoints ... ")
    ftu <- bspline_footpoint(points, 
                      c_points, 
                      degree = deg, 
                      algorithm=fp.algorithm, 
                      knots = knots,
                      tolerance=fp.tolerance, 
                      max.iterations = fp.max.iterations, 
                      initial.scan.points=fp.initial.scan.points, 
                      verbose=fp.verbose)
    ftpoints <- interpolate_bspline(ftu, c_points, degree = deg, knots=knots, unitary=F, ...)
  }else
    ftpoints <- footpoints
  
  
  print("finding approximations to turns ...")
  
  # this function transform the first term of a vectorial product into a matrix suitable for dot product
  # leaving the result unchanged
  transform.T <- function(U){
    "
    Transforms matrix of row-vectors U in a equivalent array-tensor for cross-product 
    multiplication.
    "
    if( ncol(U) != 3 ) stop("[Max] number of columns must be three")
    
    M <- array(rep(NA, times=prod(dim(U))), dim=c(dim(U)[2],dim(U)[2],dim(U)[1]))
    z <- rep(0, times=dim(U)[1])
    for(i in 1:dim(U)[2]) M[i,i,] <- z # sets the zero diagonal
    
    M[2,3,] <- U[,1]
    M[3,2,] <- -U[,1]
    
    M[3,1,] <- U[,2]
    M[1,3,] <- -U[,2]
    
    M[1,2,] <- U[,3]
    M[2,1,] <- -U[,3]
    
    M
  }
  # this function computes the matrix of cross products as a tensor product
  get.crossp.M <- function(u,v){
    U <- transform.T(u)
    uT <- to.tensor(U)
    dim(uT) <- c(x=3,y=3,z=dim(U)[3])
    
    vT <- to.tensor(v)
    dim(vT) <- c(z=dim(v)[1],x=3)
    
    R4 <- einstein.tensor(uT, vT, by="z")
    to.matrix.tensor(R4, j="y", i="z")
  }
  # this function gets the matrix dot product as a tensor product
  get.dotp.M <- function(m, p){
    A <- to.tensor( m )
    dim(A) <- c(z=nrow(p),x=3 )
    
    B <- to.tensor( p )
    dim(B) <- c(z=nrow(p), x=3)
    
    R <- mul.tensor(A, i=c("x"),B,j=c("x"), by="z")
    #print(R)                                  # INSERTED LINES
    
    as.vector(to.matrix.tensor(R, i="z"))
  }
  # module computation function on a 3 dimensional row-vector matrix
  get.modulus.M <- function(m)  as.vector(sqrt( m[,1]*m[,1]+m[,2]*m[,2]+m[,3]*m[,3] ))
  
  # multiplies each column of matrix m by the vector v (one by one)
  # get.simple.mult <- function(v, m){
  #   R <- m
  #   for(k in 1:ncol(m)) R[,k] <- R[,k] * v
  #   R
  # }
  #get.simple.mult <- function(v, m) m * v
  
  # get.scalar.mult <- function(m, p, exponent=1){
  #   v <- ( get.modulus.M(m) )^exponent
  # 
  #   get.simple.mult(v, p)
  # }

  #get.scalar.mult <- function(m, p, exponent=1) (get.modulus.M(m)^exponent) * p
  
  # 1.) GENERATES the ab initio scan points 
  # parameters of the initial scan points
  # evenly distributed, and the corresponding derivatives.
  vu  <- seq(from=scan.interval[1], to=scan.interval[2], length.out=initial.scan.points)  
  vpd <- interpolate_bspline(vu, c_points, degree=degree, derivate=1, unitary=F, ...) 
  vpdd <- interpolate_bspline(vu, c_points, degree=degree, derivate=2, unitary=F, ...) 
  
  
  # creates the matrix of tested points
  CPM <- matrix(NA, nrow=nrow(points), ncol=initial.scan.points)
  
  #projv <- ftpoints - points     # the projection vectors
  projv_prec <- ftpoints-points
  projv <- (get.modulus.M(projv_prec)^-1) * projv_prec

  
  # NOTE: this loop could be optimized using tensor multiplication ***
  CPpre <- get.crossp.M(vpd, vpdd)                                
  CP    <- (get.modulus.M(CPpre)^-1) * CPpre                   # EDITING LINE *****
  for(j in 1:nrow(points))
    CPM[j,] <- t(CP %*% matrix( projv[j,] ))                            # MODIFIED LINE ****
  
  # gets the number of columns of the starting minima matrix
  # and obtains the proper number of rows
  a <- rep(0, times = nrow(CPM))
  for( i in 2:ncol(CPM) ){
    jinds  <-  which( (CPM[,i-1] > 0 & CPM[,i] <= 0) | (CPM[,i-1] < 0 & CPM[,i] >= 0)  )
    a[jinds] <- a[jinds] + 1
  }
  max.size <- max(a)
  # builds the starting minima matrix
  a.m <- matrix(NA, nrow=nrow(points), ncol=max.size)
  a <- rep(0, times=nrow(CPM))
  for( i in 2:ncol(CPM) ){
    jinds <- which( (CPM[,i-1] > 0 & CPM[,i] <= 0) | (CPM[,i-1] < 0 & CPM[,i] >= 0) )
    a[jinds] <- a[jinds] +1
    for(j in jinds)  a.m[j, a[j]] <- vu[i] # THIS LINE requires optimization *****
    #a.m[j+a[jinds]*(ncol(a.m)-1) ] <- vu[i] # REQUIRES MODIFICATION *****
  }

  rm(a, max.size, CPM, jinds)

  
  print ("finding turns (for real) ...")
  
  # 2.) INITIALIZATION ...
  # each item in this list is a registry of the turns found with each point in the set
  rL <- list()
  for(i in 1:nrow(points)) rL[[i]] <- data.frame()
  
  
  # NOTE that these intervals may require a modification *************
  # this block might require a modification
  #mink <- min(min(knots),min(scan.interval))
  #maxk <- max(max(knots),max(scan.interval))
  # ************************************
  mink <- min(scan.interval)
  maxk <- max(scan.interval)
  
  
  # 3.-) performs Newton-Gauss minimization on each suggested parameter found near the true value
  #for(j in 1:length(u.min)){
  for(j in 1:ncol(a.m)) {
    
    uL <- as.vector(a.m[,j])
    
    tol   <- rep(Inf, times=nrow(points))
    count <- rep(0 , times=nrow(points))
    b     <- rep(TRUE, times=nrow(points))
    
    # excludes the evaluation of a point if the value in u is missing
    # it is missing if and only if no more minima are supposed to be found at that position
    b[is.na(uL)] <- FALSE
    
    bb    <- TRUE       # the global flag
    indexes <- which(b) # initializes the operative indexes, i.e. the indexes of point that will be processed
    
    while(bb)    {
      # we only select the non fitted points to improve the performance
      uLT <- uL[b]
      projvT <- projv[b,]
      tolT <- tol[b]
      countT <- count[b]
      
      # check bug line
      if(length(projvT) <= 3){ projvT <- t(matrix(projvT))
        print("correcting class on projections ...")
      }
      
      # computes the important derivative matrixes 
      vpd  <- interpolate_bspline(uLT, c_points, degree=degree, derivate=1, unitary=F, ...)
      vpdd <- interpolate_bspline(uLT, c_points, degree=degree, derivate=2, unitary=F, ...)
      vpddd <- interpolate_bspline(uLT, c_points, degree=degree, derivate=3, unitary=F, ...)
      
      # CODE TEMPORALY SHUT-DOWN ***********************
      #fd  <- get.dotp.M(get.crossp.M(vpd, vpdd), projvT)
      #fdd <- get.dotp.M(get.crossp.M(vpd, vpddd), projvT)
      # ************************************************
      
      # REPLACEMENT CODE *******************************
      
      # calculates denominator for NR routine
      alpha <- get.crossp.M(vpd, vpdd)
      beta <- get.crossp.M(vpd, vpddd)

      r1 <-  (get.modulus.M(alpha)^-1)*beta         #get.scalar.mult(alpha, beta, exponent=-1)
      r2.1 <- (get.modulus.M(alpha)^-3)*alpha       #get.scalar.mult(alpha, alpha, exponent=-3)

      r2.2 <- get.dotp.M(beta, alpha)
      #return(r2.2)     # DEBUGGING LINE *************************+
      #r2 <- get.simple.mult(r2.2, r2.1)
      r2 <- r2.2 * r2.1

      fdd <- get.dotp.M((r1 - r2), projvT)   # it works if r1 + r2
      
      # calculates the numerator
      r3 <- (get.modulus.M(alpha)^-1)*alpha  #get.scalar.mult(alpha, alpha, exponent=-1)
      fd <- get.dotp.M(r3 , projvT)
      
      # ************************************************
      
      # computing the paramerter increments
      inc_uLT <- -fd / fdd
      countT <- countT+1
      tolT <- fd * fd #abs(fd)
      
      #print("*")      # TEST LINE
      #print(tolT) 
      
      uLT <- uLT + inc_uLT
      
      # here we have to update the original lists: uL, tol, count, and b (bb will be supressed
      # later on). The update should be performed when a b[i] turns to FALSE, due to 
      # a tolerance fail or other of the listed criteria (see below)
      
      exc.lims <- uLT < mink | uLT > maxk
      b[indexes] <- tolT > tolerance & countT < max.iterations & !exc.lims
      uLT[exc.lims] <- NA         # points exceding the limits are excluded
      countT[exc.lims] <- NA      # idem.
      uL[indexes] <- uLT
      count[indexes] <- countT
      tol[indexes] <- tolT
      
      # updates the operative indexes
      indexes <- which(b)
      bb <- length(indexes) > 0 # full termination condition for the first suggested minimum
    }
    
    # updates the list of data frames for each point, indicating the points
    for(i in 1:nrow(points)) {
      if( !is.na(uL[i]) ) rL[[i]] <- rbind(rL[[i]], c(uL[i], count[i], tol[i]) )
      if( j >= ncol(a.m) && nrow(rL[[i]]) != 0 ) names(rL[[i]]) <- c("u.found","iterations", "tolerance")
    }
    
  }
  
  
  
  rL
  
  
}
