#' @title Footpoint computation
#' @description \code{bspline_footpoint} computes the parameter values for the
#' footpoints of a point set onto the specified spline curve by numerical approximation.
#' @details The function initializes the parameter vector \code{u} by brute-force search of closer footpoint candidates whose
#' number is specified by the user. The candidate footpoints are evenly distributed throughout the curve by default and a Newton-Raphson minimization 
#' routine ensues in which \code{u} is updated each round through the computation of an increment using different types of functionals.
#' For a data point and its footpoint these functionals consider the orthogonality of the projection, their distance or 
#' a combination of both. When both criteria are mixed, their relative contribution is controled by the parameter \code{alpha}, 
#' e.g. a zero \code{alpha} indicates that orthogonality is not taken into account. 
#
#' Initialization can be controled by the user through the \code{prompt.mode}. 
#' @param x The matrix of row-vector positions that correspond to the data points.
#' @param c_point The matrix of row-vector positions that correspond to the control points.
#' @param u The initial parameter values that will be improved, equal to the number of rows in \code{x}. If 
#' NULL one is provided by default.
#' @param degree The spline polynomial degree.
#' @param tolerance The value below which the numerical approximation will exit. 
#' @param step This factor controls how much of the parameter increase resulting from each round in the Newton-Raphson
#' minimization is actually being applied for the next round.
#' @param initial.scan.points The number of evenly distributed points onto the spline curve that will be checked to 
#' generate to set up the parameter vector \code{u}. It only affects whenever \code{u} is NULL.
#' @param knots sequence of knots. (NOTE that this is not the support as in other functions found throughout 
#' the package \code{modgeom})
#' @param verbose if \code{TRUE} it displays progress information.
#' @param prompt.mode the mode of projection prompt (i.e. intial parameter values for the algorithm). Possible values are:
#' \describe{
#' \item{\strong{distance}}{Candidates footpoints are assigned by their least distances to the data points}
#' \item{\strong{orthogonal}}{Candidates footpoints are assigned by their degree of orthogonality to the data points}
#' \item{\strong{directional}}{Candidate footpoints are assigned by the minimum value of the directional functional (see \code{algorithm}) to
#' the data points}
#' }
#' @param prompt.begin The initial parameter value for the prompting points.
#' @param prompt.end The ending parameter value for the prompting points.
#' @param algorithm The functional used for the Newton-Raphson approximation. Possible values are:
#' \describe{
#'  \item{\strong{Rogers}}{It makes an approximation which employs an expasion series around the current projection (Rogers 1989). }
#'  \item{\strong{orthogonal}}{It finds the orthogonal points only. The functional tests the inner product of the 
#'  tangent and the projection vector.}
#'  \item{\strong{directional}}{The directional algorithm minimizes a function of the distance and orthogonality which is
#'  controled by a factor
#'  \code{alpha} (e.g. \code{alpha} = 0 only distance is considered). NOTE: values of \code{alpha} different from zero
#'  are experimental.}
#'  \item{\strong{none}}{Only prompt values are kept.}
#' }
#' @param max.iterations Maximum number of allowed iterations.
#' @param alpha A value between 0 and 1 that controls the degree of orthogonality that will be included in the directional modes (see \code{algorthim}) and
#' \code{prompt.mode}). For other modes this parameter has no effect.
#' @param vz (experimental) when orthogonality is considered, this vector marks the z direction.
#' @param ... More arguments passed to \code{\link{interpolate_bspline}}
#' @return A vector of parameter values corresponding to the computed footpoints.
#' @author Máximo Sánchez-Aragón
#' @example examples/bspline_footpoint_example1.R
#' @example examples/bspline_footpoint_example2.R
#' @export

bspline_footpoint <- function(x, c_point, 
                               u=NULL, 
                               degree = 3,
                               tolerance=1e-9,
                               step=0.2,
                               initial.scan.points=nrow(c_point)*8,
                               knots = seq(from = 0, to = 1, length.out=nrow(c_point) - degree + 1),
                               verbose=FALSE,
                               prompt.mode="distance",
                               prompt.begin=0,
                               prompt.end=1,
                               algorithm="Rogers",
                               max.iterations = 100,
                               alpha = 0,
                               vz = c(0,0,1),
                               ...){
  
  
  # takes native and is.periodic properties from the arguments passed to interoplate_bspline
  L <- as.list( sys.call() )
  
  fp.native <- as.logical(as.character(L$native))
  fp.is.periodic <- as.logical(as.character(L$is.periodic))
  
  if(is.null(L$native)){ fp.native <- TRUE; }
  if(is.null(L$is.periodic)) { fp.is.periodic <- FALSE; }

  if(!fp.native){
    

  # INITIALIZATION ...
  # if u is null, a default knot sequence is provided
  # with a number of points equal to 'initial.scan.points'    

  if(is.null(u)){ 
    
    if(verbose){
      print("No parameter vector found (u = NULL). Prompting points ...")
      print(paste("prompting mode:",prompt.mode))
    }
    
    
    fu.distance <- function(){
      param <- seq(from=prompt.begin, to=prompt.end, length.out=initial.scan.points)
      
      # computes first guess for u ...  
      v <- interpolate_bspline(param, c_point, degree = degree, 
                               knots = knots,...)

      for(j in 1:nrow(x)){
        #print(paste(j,"/",nrow(x))) # test line
        # inner loop -->> finds the minimum distance on the current point
        m_dis <- Inf
        m_ind <- -1
        for(i in 1:nrow(v)){
          dis <- sum((v[i,]-x[j,])^2)
          if(dis < m_dis){
            m_dis <- dis 
            m_ind <- i
          }                
        }    
        u[j] <- param[m_ind]   
      }
      u
    }
    fu.orthogonal <- function(){
      param <- seq(from=prompt.begin, to=prompt.end, length.out=initial.scan.points)
      
      # computes first guess for u ... 
      v <- interpolate_bspline(param, c_point, degree = degree, 
                               knots = knots, derivate=0, unitary=F,...)
      vd <- interpolate_bspline(param, c_point, degree = degree, 
                                knots = knots, derivate=1, unitary=T, ...)
      for(j in 1:nrow(x)){    
        # inner loop -->> finds the minimum distance on the current point
        m_dis <- Inf
        m_ind <- -1
        for(i in 1:nrow(v)){
          v0 <- v[i,]-x[j,]        
          dis <- sum( (vd[i,] * v0)^2 )           
          if(dis < m_dis){
            m_dis <- dis 
            m_ind <- i
          }                
        }    
        u[j] <- param[m_ind]   
      }
      u
    }
    fu.directional <- function(){
      param <- seq(from=prompt.begin, to=prompt.end, length.out=initial.scan.points)
      
      # computes first guess for u ...  
      v <- interpolate_bspline(param, c_point, degree = degree, 
                               knots = knots,  ...)
      for(j in 1:nrow(x)){    
        # inner loop -->> finds the minimum distance on the current point
        m_dis <- Inf
        m_ind <- -1
        for(i in 1:nrow(v)){
          dis <- (1-alpha) * sum((v[i,]-x[j,])^2) + 
                    alpha  * sum((v[i,]-x[j,]) * vz)
          if(dis < m_dis){
            m_dis <- dis 
            m_ind <- i
          }                
        }    
        u[j] <- param[m_ind]   
      }
      u
    }
    
    if(prompt.mode == "distance")
      u <- fu.distance()
    else if(prompt.mode == "orthogonal") 
      u <- fu.orthogonal()
    else if (prompt.mode == "directional")
      u <- fu.directional()
    else
      stop("[Max] unrecognized prompt.mode selection")
  }
    
  #print(u)
  
  if(length(u) != nrow(x))
    stop("[Max] the length of the data points and the parameter sequence is not equal")
    
  if(verbose){
    print("projecting points onto curve")
    print(paste("prompt mode:",prompt.mode))
    print(paste("algorithm:",algorithm))
    print("point: <index of the datum>, it.: <number of iterations>")
  }
  
  # qs stores the initial interpolations  
  mink <- min(knots)
  maxk <- max(knots)
  
  # ALGORITHM 1 : by default is the distance
  if(algorithm == "Rogers")
  for(i in 1:nrow(x)){ # test each point individually    
            
    t <- u[i]
    tol <- Inf 
    count <- 0
    while(tol > tolerance){ # this loop must test the iteration quality
                  
      # evaluates the point t and its first derivative
      c <- interpolate_bspline(t, c_point, degree = degree, knots = knots, ...)
      cd <- interpolate_bspline(t, c_point, degree = degree, derivate = 1,  ...)
                  
      # computes the intersection point with the tangent
      # to test closeness          
      qs <- c + ( sum( cd * ( x[i,] - c ) ) / sum( cd * cd ) ) * cd
            
      p = step * sum( cd * (x[i,] - c) ) / sum(cd * cd) 
      dt = p
      
      #print(dt);
      
      if(t + dt < mink || t + dt > maxk){
        if(verbose)
          print(paste("projected to tip, i:",i))
        break   
      }
        
      t <- t + dt       
      tol <- sqrt( sum( (qs - c)^2 ) )
      count <- count + 1  
      
      if(count > max.iterations)
        break
      
    }
    
    u[i] <- t  
    
    if(verbose)
    print(paste("point:",i,", it.:",count, "(",algorithm,")"))
  }        
  
  # ALGORITHM 2 : none
  else if(algorithm == "none")
  print("dull algorithm selected ...")
  
  # ALGORITHM 3 : orthogonal projection finding
  else if(algorithm == "orthogonal")
  for(i in 1:nrow(x)){
    
    t <- u[i]
    
    # RECURRENCE ...  
    # 2.-) performs Newton-Gauss minimization on each point
    tol <- Inf 
    count <- 0 
    b <- TRUE
    while(tol > tolerance){ # this loop must test the iteration quality
      
      # evaluates the point     
      x0   <- interpolate_bspline(t, c_point, degree = degree, derivate = 0, 
                                  knots = knots,  ...)
      xd   <- interpolate_bspline(t, c_point, degree = degree, derivate = 1, 
                                  knots = knots, ...)
      xdd  <- interpolate_bspline(t, c_point, degree = degree, derivate = 2, 
                                  knots = knots, ...)
      p <- x0 - x[i,]
                  
      # newtonian expression
      fd  <- as.numeric( p %*% t(xd) )
      fdd <- as.numeric( p %*% t(xdd) + xd %*% t(xd) ) 
      
      inc_t   <- -fd / fdd
      
      count <- count+1
      tol <- (p %*% t(xd))^2
      
      t <- t + inc_t  
      
      if(count > max.iterations || t < mink || t > maxk){
        warning("[Max] non-convergent minimum found. Discarding ... ")
        b <- FALSE
        break
      }
      
    }
    
    # adds the found minimum
    if(b) { 
        u[i] <- t    
        if(verbose)
          print(paste("point:",i,", it.:",count, "(",algorithm,")"))
      }
    # ENDS of RECURRENCE ...        
  }
  
  # ALGORITHM 4 : direction preferred
  else if(algorithm == "directional")
  for(i in 1:nrow(x)){
      
      t <- u[i]
      
      # RECURRENCE ...  
      # 2.-) performs Newton-Gauss minimization on each point
      tol <- Inf 
      count <- 0 
      b <- TRUE
      while(tol > tolerance){ # this loop must test the iteration quality
        
        # evaluates the point     
        x0   <- interpolate_bspline(t, c_point, degree = degree, derivate = 0, 
                                    knots = knots, ...)
        xd   <- interpolate_bspline(t, c_point, degree = degree, derivate = 1, 
                                    knots = knots, ...)
        xdd  <- interpolate_bspline(t, c_point, degree = degree, derivate = 2, 
                                    knots = knots, ...)

        x0 <- as.numeric(x0)
        xd <- as.numeric(xd)
        xdd <- as.numeric(xdd)
        
        p <- x0 - x[i,]
                
        # newtonian expression
        fd  <- (1-alpha) * sum( p * xd ) + 
                  alpha * sum(p * vz) * sum(xd * vz)
        fdd <- (1-alpha) * ( sum(xd * xd) + sum(p * xdd) ) + 
                  alpha  * ( sum( xd * vz )^2 + sum(p * vz) * sum(xdd * vz) )
        
        inc_t   <- -fd / fdd
        
        count <- count+1
        tol <- sum(p * xd)^2
        
        t <- t + inc_t  
        
        #print(paste(xd," ** "), sep="");
        
        if(count > max.iterations || t < mink || t > maxk){
          warning("[Max] non-orthogonal minimum found. discarding ... ")
          b <- FALSE
          break
        }
        #break;
      }
      
      # adds the found minimum
      if(b) { 
        u[i] <- t    
        if(verbose)
          print(paste("point:",i,", it.:",count, "(",algorithm,")"))
      }
      # ENDS of RECURRENCE ...        
    }
  
  else
  stop("[Max] unrecognized algorithm selected") 
  
  return(u)
  
  }else{
    
    # algorithm
    # Rogers R, orthogonal O, directional D, none N
    alg.native = 0 # none
    if(algorithm == "Rogers") alg.native = 1
    else if(algorithm == "orthogonal") alg.native = 2
    else if(algorithm == "directional") alg.native = 3
    
    # prompt.mode
    pm = 0 # distance
    if(prompt.mode == "orthogonal") pm = 1 
    else if(prompt.mode == "directional") pm = 2
    
    uu <- rep(0,times=nrow(x))
    uun <- length(u)
    nrowcp <- nrow(c_point)
    ncolcp <- ncol(c_point)
    knots_n <- length(knots)
    
    #warning("Periodic splines do not perfom well under native mode. Change mode for better performance.")
    #warning("alpha is set to 0 in native mode. Wait for new releases.")
    #warning("Previous values for the parameter are discarded in native mode. Wait for new releases.")
    
    #isp <- 0 # is.periodic
    alpha_alt <- 0
    
    result<- .C("bspline_footpoint", 
        as.numeric(uu),
        as.integer(uun),
        as.numeric(x),
        as.integer(ncol(x)),
        as.integer(nrow(x)),
        as.numeric(c_point),
        as.integer(nrowcp),
        as.integer(ncolcp),
        as.integer(degree),
        as.numeric(knots),
        as.integer(knots_n),
        as.integer(pm),
        as.numeric(prompt.begin),
        as.numeric(prompt.end),
        as.integer(initial.scan.points),
        as.integer(alg.native),
        as.integer(fp.is.periodic),
        as.numeric(tolerance),
        as.numeric(alpha_alt),
        as.numeric(vz), 
        as.integer(max.iterations),
        as.numeric(step)
       )
    
    return(as.numeric(result[[1]]) )
    
  }
  
}
