

# USER PARAMETERS -------------------------
require(modgeom)
# load("data/bspline_parity_turns_data") # comes from the processing in bspline2_test8 example.
# 
# 
# c_points[1,3] <- c_points[1,3] +3
# c_points[3,2] <- c_points[3,2] -15
# c_points[3,1] <- c_points[3,1] +3
# c_points[5,1] <- c_points[5,1] +10
# 
# deg <- 3
# c_points <- rbind(c_points,c(2.5,-1,1.5))

#load("inst/extdata/bspline_uz_turns_data") # comes from the processing in bspline2_test8 example.
load(paste(.libPaths()[1],"/modgeom/extdata/bspline_uz_turns_data",sep=""))

deg <- 3
#c_points <- rbind(c_points,c(2.5,-1,2))
vec_shift_value=2
c_points <- rbind( c(0,0,-.5),c(-1,0,-.5),c(-1,1,0), c(2,1,1) ,c(4,1,0),c(4,0,-.5),c(3,0,-.5) )

for(j in 1:nrow(c_points))
  c_points[j,3] <- c_points[j,3] + vec_shift_value

# load("data/bspline_parity_turns_data")
# deg <- 3
# c_points <- rbind(c_points,c(2.5,-1,1.5))


u <- seq(from = 0, to = 1, length.out=100)
points <- interpolate_bspline(u, c_points, degree = deg)
ref_point <- c(4,1,1)    # the reference observation point
ref_point_alt <- ref_point[length(ref_point):1]

# we need to define also the points that will be projected
n <- 200 # the number of points that will be projected
dev <- 1  # the deviation in the test points
cell.points <- cbind( runif(2*n, min=-dev, max=dev), 
                 runif(2*n, min=-dev, max=dev), 
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


multiplier <- 100
start.time1 <- Sys.time()
signv <- bspline_parity_sign (
          c_points,
          cell.points,
          fp.u = NULL,
          degree = deg,
          tolerance=1e-15,
          fp.tolerance=1e-12,
          max.iterations = 1500,
          initial.scan.points=nrow(c_points)*multiplier,
          fp.verbose=TRUE
)
signv_xor <- signv[[2]]
u_turn2_data <- signv[[1]]
end.time1 <- Sys.time()
end.time1 <- start.time1 - end.time1
print(paste("procedure lasted: ",-end.time1))


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
    clso <- cls[j]
    if(nrow(u_turn2_data[[j]]) > 0){
      
      tp <- interpolate_bspline (u_turn2_data[[j]]$u.found, c_points, degree=deg)
      
      clst <- cls[j]
      
      
      
      points3d( tp,  col="green", size = 6  )
    }
    lines3d( rbind(cell.points[j,],ftpoints[j,]), col=clso, lwd=.01)
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
sq <- c(0:90,89:1)
i = sq[sq %% 10 == 0][1]
cc <- cc+1
vpar1 <- c(i,i)
cst <- as.character(cc)
if(nchar(cst) == 1) cst <- paste("00",cst, sep="")
if(nchar(cst) == 2) cst <- paste("0",cst, sep="")
  
d3graphics(signv_xor, paste("xor_scalar_parity_",cst,sep=""),vpar=vpar1)

rm(i,cc)


