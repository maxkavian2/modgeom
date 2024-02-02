

# USER PARAMETERS -------------------------
require(modgeom)
#load("inst/extdata/bspline_parity_turns_data")
load(paste(.libPaths()[1],"/modgeom/extdata/bspline_parity_turns_data",sep=""))


deg <- 3


c_points <- rbind(c_points,c(2.5,-1,1.5))

u <- seq(from = 0, to = 1, length.out=100)
points <- interpolate_bspline(u, c_points, degree = deg)
ref_point <- c(4,1,1)    # the reference observation point
ref_point_alt <- ref_point[length(ref_point):1]

# we need to define also the points that will be projected
n <- 50 # the number of points that will be projected
dev <- 1  # the deviation in the test points
cell.points <- cbind( runif(2*n, min=-dev, max=dev), 
                      runif(2*n, min=-dev*2, max=dev*2), 
                      runif(2*n, min=-dev, max=dev) )

#rm(n)
rM <- c()
for(i in 1:n)
  rM <- rbind(rM, ref_point)
for(i in (n+1):(2*n))
  rM <- rbind(rM, ref_point_alt)

#cell.points <- rbind(cell.points, cell.points)
cell.points <- cell.points + rM
rm(rM,i)
rm(n)

# checks the turns function detection ----------------
ftpt_u <- bspline_footpoint(cell.points, c_points, degree = deg, algorithm="directional",
                            tolerance=1e-12, max.iterations = 800, initial.scan.points=nrow(c_points)*100)

ftpoints <- interpolate_bspline(ftpt_u, c_points, degree = deg)

multiplier <- 300
start.time1 <- Sys.time()

u_turn2_data <- bspline_parity_turns( # u_turn2_data
  c_points,
  cell.points,
  #ref_vec,   # TEST argument
  footpoints = ftpoints,
  degree = deg,
  tolerance=1e-10,
  max.iterations = 500,
  initial.scan.points=nrow(c_points)*multiplier,
  fp.verbose=TRUE
)
end.time1 <- Sys.time()
end.time1 <- start.time1 - end.time1
print(paste("procedure lasted: ",-end.time1))


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
  
  as.vector(to.matrix.tensor(einstein.tensor(A, B, by="z"), i="z"))
  #einstein.tensor(A, B, by="z.i")
}

# gets the dot product sign
vp   <- interpolate_bspline(ftpt_u, c_points, degree=deg, derivate=0, unitary=F)
vpd  <- interpolate_bspline(ftpt_u, c_points, degree=deg, derivate=1, unitary=F)
vpdd <- interpolate_bspline(ftpt_u, c_points, degree=deg, derivate=2, unitary=F)

p   <- ftpoints - cell.points#cell.points 

pd_x_pdd <- get.crossp.M(vpd, vpdd)
fd  <- get.dotp.M(pd_x_pdd, p)    

# tries to find the sign -----------------------------
signv_1 <- c()
signv_2 <- c()
# signv_xor <- c()
# pd <- interpolate_bspline(ftpt_u, c_points, degree=deg, derivate = 1)
# pa <- ftpoints - cell.points
signv_1 <- fd >= 0
for(i in 1:length(u_turn2_data) ){
  ct <- u_turn2_data[[i]]$u.found #u_turn2_data[[i]]$u.found
  signv_2[i] <- length(ct[ct <= ftpt_u[i]]) %% 2 == 1
  #signv_xor[i] <- xor(vp < 0, length(ct[ct < ftpt_u[i]]) %% 2 == 1) # finds the number of turns before the projection point
}
signv_xor <- xor(signv_1, signv_2)


# GRAPHICS --------------------------------------
d3graphics <- function (signv, tag, vpar = c(0,0)){
  #attach(loaded.data)
  require(rgl)
  
  cls <- rep("red", times=nrow(cell.points))
  cls[signv] <- "blue"
  

  
  plot3d(loaded.data$R_fuzzy, type="n", xlab="x", ylab="y", zlab="z")
  #points3d(R_fuzzy, col="yellowgreen", size=3, pch=21, bg="red" )
  par3d(windowRect=c(150, 150, 700, 700))
  #rgl.viewpoint(vpar[1], vpar[2])
  view3d(vpar[1],vpar[2])
  #par3d(zoom = 1)
  
  
  #points3d(t(as.matrix(ref_point)), col="red", size=8)
  points3d(c_points, col="lightblue4", size=4)
  lines3d(c_points, col="lightblue4", lty="dotted", lwd=.5)
  lines3d(points, col="lightblue4", lwd=2) # points of the spline curve
  sgs <- c();
  
  for(j in 1:length(u_turn2_data)){
    if(nrow(u_turn2_data[[j]]) > 0){
      
      tp <- interpolate_bspline (u_turn2_data[[j]]$u.found, c_points, degree=deg)
      
      clst <- cls[j]
      
      lines3d( rbind(cell.points[j,],ftpoints[j,]), col=clst, lwd=.01)
      
      points3d( tp,  col="green", size = 6  )
    }
  }
  
  points3d( cell.points, col=cls , size = 6 )
  v0 <- interpolate_bspline (0, c_points, degree=deg, derivate=0)
  v1 <- interpolate_bspline (0, c_points, degree=deg, derivate=1)
  v1 <- v1 / (2*sqrt(sum(v1^2)) )
  v1 <- v0 + v1
  arrow3d(v0, v1, 1/20, width=1/4, thickness=1/10)
  rm(tp)
  
  #detach(loaded.data)
}

cc <- 0
sq <- 50 #c(0:90,89:1)
for(i in sq[sq %% 1 == 0] ){
  cc <- cc+1
  vpar1 <- c(i,i)
  cst <- as.character(cc)
  if(nchar(cst) == 1) cst <- paste("00",cst, sep="")
  if(nchar(cst) == 2) cst <- paste("0",cst, sep="")
  
  #d3graphics(signv_1, paste("only_scalar_",cst,sep=""), vpar=vpar1)
  #d3graphics(signv_2, paste("only_parity_",cst,sep=""),vpar=vpar1)
  d3graphics(signv_xor, paste("xor_scalar_parity_",cst,sep=""),vpar=vpar1)
}
rm(i,cc)


