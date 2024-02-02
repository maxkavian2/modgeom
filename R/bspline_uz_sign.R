#'@title Signed positions of points projections to spline curves (TG method)
#'@description This function defines a sign according to the parity of turns and the plane defined by the footpoint projections
#' and the spline tangent vectors (at the footpoint), described in the PhD thesis as the tangent plane method (TG).
#' 
#' This function specifically depends on tensorA package. The current evaluation only works
#' for three dimensional space, i.e. all row-vector matrices defining positions must have three columns.
#' @param c_point The current control points, expressed as a matrix of row vectors.
#' @param points The projected points, as a matrix of row vectors.
#' @param fp.u Parameter vector corresponding to the footpoints. If \code{NULL} the algorithm would try to generate
#' them with parameters fp.* (below). Check the function \code{\link{bspline_footpoint}} for more information about
#' these parameters.
#' @param degree The maximum degree of the spline curve.
#' @param tolerance The threshold, as the absolute value of the dot product, below which
#' the routine ends for each point. The tolerance is applied to each tested point independently.
#' @param max.iterations The maximum number of allowed iterations prior to termination. Each test point is evaluated
#' independently.
#' @param initial.scan.points Number of initial points that are scanned when \code{fp.u} is not \code{NULL}.
#' @param knots The knot sequence of the spline.
#' @param ap_ref The direction (vector of length 3) provided for sign correction, i.e. this vector is used to normalize the 
#' sign criterium to one direction. This way results from differnet splines may be compared afterwards.
#' @param uz The unitary z vector used as criterium of projection for the sign,
#' @param scan.interval  \emph{[EXPERIMENTAL]} The interval within which the turns are to be found. This interval
#' would be divided as many times as indicated by \code{initial.scan.points}. 
#' @param sign.op If 'normal' the full sign operation is retrieved. 'parity' shows only the boolean for
#' the turn parity only. 'other' represents the other component in the operation, generally
#' a differential geometry operation. In this case, it is the sign of the vector product of the unitary z vector and 
#' the tangent.
#' @param ... Argument passed to \code{\link{interpolate_bspline}}.
#' @return A \code{list} of two elements. The first element contains a list of \code{data.frame}s with the parameter values 
#' corresponding to the turns and additional information. Each \code{data.frame} also contains information about the total 
#' number of cycles and the final tolerance value, the latter 
#' being an indication of the turn precision. The second element is a vector with the final sign of all projected points.
#' @example examples/bspline_uz_sign_example1.R
#' @author Máximo Sánchez-Aragón
#' @export

bspline_uz_sign <- function(  
                  c_points,
                  points,
                  fp.u = NULL,
                  degree = 3,
                  tolerance=1e-9,
                  max.iterations = 100,
                  initial.scan.points=nrow(c_points)*8,
                  knots = seq(from = 0, to = 1, length.out=nrow(c_points) - degree + 1),
                  ap_ref = NULL,
                  uz = c(0,0,1),
                  scan.interval = c(min(knots),max(knots)),
                  fp.tolerance = 1e-12,
                  fp.max.iterations = 800,
                  fp.initial.scan.points=nrow(c_points)*8,
                  fp.algorithm = "directional",
                  fp.verbose = FALSE,
                  sign.op = "normal",
                  ...
){

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
  if(is.null(fp.u)){
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
  }else{
    ftu <- fp.u
    ftpoints <- interpolate_bspline(ftu, c_points, degree = deg, knots=knots, unitary=F, ...)
  }
  
  if(length(ftu) != nrow(points))
    stop("[Max] the number of footpoint parameters must be the same as the number of points")
  
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
    dim(A) <- c(z=nrow(p), x=3 )
    
    B <- to.tensor( p )
    dim(B) <- c(z=nrow(p), x=3)
    
    R <- mul.tensor(A, i=c("x"),B,j=c("x"), by="z")
    
    
    MR <- to.matrix.tensor(R, i="z")
    
    as.vector(MR)
  }
  # module computation function on a 3 dimensional row-vector matrix
  get.modulus.M <- function(m)  as.vector(sqrt( m[,1]*m[,1]+m[,2]*m[,2]+m[,3]*m[,3] ))
  
  # [ explain later ]
  make.test.space <- function(u, v, uz_vector=c(0,0,1) ){
    
    U <- transform.T(u)
    uT <- to.tensor(U)
    dim(uT) <- c(x=3,y=3,z=dim(U)[3])
    
    vT <- to.tensor(v)
    dim(vT) <- c(l=dim(v)[1],x=3)
    
    UvT <- einstein.tensor(uT, vT, by=c("z","l") )
    
    #print(UvT)
    
    zT <- to.tensor(matrix(uz_vector))
    dim(zT) <- c(y=length(uz_vector))
    
    zUvT <- einstein.tensor(UvT, zT, by=c("z","l"))
    R <-  to.matrix.tensor(zUvT, i="z")
    
    
    
    R
    
  }
  
  
  # 1.) GENERATES the ab initio scan points 
  # parameters of the initial scan points
  # evenly distributed, and the corresponding derivatives.
  vu  <- seq(from=scan.interval[1], to=scan.interval[2], length.out=initial.scan.points)  
  vpd <- interpolate_bspline(vu, c_points, degree=degree, derivate=1, unitary=F, ...) 
  # vpdd <- interpolate_bspline(vu, c_points, degree=degree, derivate=2, unitary=F, ...) 
  
  
  
  # creates the matrix of tested points
  CPM <- matrix(NA, nrow=nrow(points), ncol=initial.scan.points)
  
  #projv <- ftpoints - points     # the projection vectors
  #projv_prec <- ftpoints-points
  #projv <- (get.modulus.M(projv_prec)^-1) * projv_prec
  projv <- ftpoints-points
  uzM <- t(matrix(rep(uz, times=nrow(projv)), nrow=ncol(projv), ncol=nrow(projv)))
  CPM <- make.test.space(projv, vpd, uz_vector = uz)

  # NOTE: this loop could be optimized using tensor multiplication ***
  # THIS BLOCK WAS COMMENTED, does not seem to have required variables in what follows
  #CPpre <- get.crossp.M(vpd, vpdd)                                
  #CP    <- (get.modulus.M(CPpre)^-1) * CPpre                   
  #for(j in 1:nrow(points))
  #  CPM[j,] <- t(CP %*% matrix( projv[j,] ))                   
  
  # NOTE: this loop could be optimized using tensor multiplication ***
  #projM <- t(matrix(NA, ncol=ncol(projv), nrow=length(vu)))

  #for(j in 1:nrow(points))
    #projM[] <- as.vector(projv[j,])
    #CPM[j,] <- get.dotp.M( get.crossp.M(vpd, vpdd), t(projM) )
    #CPM[j,] <- t(get.crossp.M(vpd, vpdd) %*% matrix( projv[j,] ))
    #CPM[j,] <- t( get.crossp.M(projv, vpd) %*% matrix(uzM[j,]) )
  # ***********
  
  
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
    
    uL0 <- uL
    
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
      uzMT <- uzM[b,]  # dull selection
      
      # check bug line
      if(length(projvT) <= 3){ projvT <- t(matrix(projvT))
        print("correcting class on projections ...")
      }
      
      if(length(uzMT) <= 3){ uzMT <- t(matrix(uzMT))  
        print("correcting class on projections ...")
      }
      
      # computes the important derivative matrixes 
      vpd  <- interpolate_bspline(uLT, c_points, degree=degree, derivate=1, unitary=F, ...)
      vpdd <- interpolate_bspline(uLT, c_points, degree=degree, derivate=2, unitary=F, ...)
      #vpddd <- interpolate_bspline(uLT, c_points, degree=degree, derivate=3, unitary=F, ...)
      
      #fd  <- get.dotp.M(get.crossp.M(vpd, vpdd), projvT)
      #fdd <- get.dotp.M(get.crossp.M(vpd, vpddd), projvT) 
      
      # REPLACEMENT CODE *******************************
      
      # calculates denominator for NR routine
      alpha <- get.crossp.M(projvT, vpd)
      beta <- get.crossp.M(projvT, vpdd)
      
      r1 <-  (get.modulus.M(alpha)^-1)*beta         #get.scalar.mult(alpha, beta, exponent=-1)
      r2.1 <- (get.modulus.M(alpha)^-3)*alpha       #get.scalar.mult(alpha, alpha, exponent=-3)
      
      r2.2 <- get.dotp.M(beta, alpha)
      #return(r2.2)     # DEBUGGING LINE *************************+
      #r2 <- get.simple.mult(r2.2, r2.1)
      r2 <- r2.2 * r2.1
      
      fdd <- get.dotp.M((r1 - r2), uzMT)   # it works if r1 + r2
      
      # calculates the numerator
      r3 <- (get.modulus.M(alpha)^-1)*alpha  #get.scalar.mult(alpha, alpha, exponent=-1)
      fd <- get.dotp.M(r3 , uzMT)
      
      # ************************************************
      
      # computing the parameter increments
      inc_uLT <- -fd / fdd
      countT <- countT+1
      tolT <- fd * fd #abs(fd) 
      
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
      if( !is.na(uL[i]) ) rL[[i]] <- rbind(rL[[i]], c(uL[i], count[i], tol[i], uL0[i]) )
      if( j >= ncol(a.m) && nrow(rL[[i]]) != 0 ) names(rL[[i]]) <- c("u.found","iterations", "tolerance", "u0")
    }
    
  }
  
  
  # gets the dot product sign
  #vp   <- ftpoints interpolate_bspline(ftpt_u, c_points, degree=deg, derivate=0, unitary=F)
  vpd  <- interpolate_bspline(ftu, c_points, degree=deg, knots=knots, derivate=1, unitary=F, ...)
  #vpdd <- interpolate_bspline(ftu, c_points, degree=deg, knots=knots, derivate=2, unitary=F, ...)
  
  p_x_pd <- get.crossp.M(projv, vpd)
  p_x_pd <- (get.modulus.M(p_x_pd)^-1) * p_x_pd 
  fd  <- get.dotp.M(p_x_pd, uzM)
  
  # it operates to find the sign, through the concept of curve parity.
  signv_1 <- c()
  signv_2 <- c()
  signv_1 <- fd >= 0
  for(i in 1:length(rL) ){
    ct <- rL[[i]]$u.found 
    signv_2[i] <- length(ct[ct <= ftu[i]]) %% 2 == 1
  }
  if(sign.op == "normal")
    signv_xor <- xor(signv_1, signv_2)
  else if(sign.op == "parity")
    signv_xor <- signv_2
  else if (sign.op == "other")
    signv_xor <- signv_1
  else
    stop("[Max says] the specified sign operation does not exist!")
    
  
  # the parity sign has to be modified to account for the actual relative orientation of the curve versus
  # the ap_axis. These extra lines modifies the result accordingly
  
  if(!is.null(ap_ref)){
    
    # TEST 1.)  works except for 1 sample
    
    # dv_v1 <- interpolate_bspline(c(0,.5,1), c_points, degree=deg, knots=knots, derivate=0, unitary=F, ...)
    # dv_v  <- t(matrix(dv_v1[3,] - dv_v1[1,]))
    # dv_vd <- dv_v1[2,] - dv_v1[1,]
    # ap_v  <- t(matrix(ap_ref))
    # 
    # # WARNING: we are deciding that the sign is mainly aligned to the z dimension; this is no necessarily so
    # warning("AP reference correction will be taken using z coordinate cross product. Note that this has not necessarily be so.")
    # picked.dimension <- 3
    # if( xor(get.crossp.M( dv_v, ap_v)[picked.dimension] <= 0, dv_vd[picked.dimension] <= 0)   ) 
    # signv_xor <- !signv_xor
    
    # TEST 2.)
    dv_v1 <- interpolate_bspline(c(0,0.5,1), c_points, degree=deg, knots=knots, derivate=0, unitary=F, ...)
    dv_v  <- t(matrix(dv_v1[3,] - dv_v1[1,]))
    dv_vd <- t(matrix(dv_v1[2,] - dv_v1[1,]))
    ap_v  <- t(matrix(ap_ref))
    
    dv_x_ap <- get.crossp.M(dv_v, ap_v)
    if( as.vector( get.dotp.M(dv_vd, dv_x_ap) ) <=0 )
      signv_xor <- !signv_xor
    
    
  }
  
  
  
  result <- list()
  result[[1]] <- rL
  result[[2]] <- signv_xor
  
  attr(result[[1]], "description") <- "turns found the curve, as a parameter"
  attr(result[[2]], "description") <- "the result of the xor operation (parity sign)"
  
  result
  
}

