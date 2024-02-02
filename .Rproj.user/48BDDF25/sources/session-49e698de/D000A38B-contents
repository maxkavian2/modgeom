

# tests the projection quality of data points
# onto a B-spline curve, user defined

# USER PARAMETERS ---------------------------------
require(modgeom)

# the positions which are interpolated
position <- seq(from = 0, to = 1, length.out=800)

# the degree
deg <- 3

c_point <- t( rbind(c(0,1,2,3,3), c(0,3,3,3,3), c(2,1,1,3,4)) )
#c_point <- rbind(c_point, c_point[1:deg,])

# the number of data points
n_data <- 40

# the data points
d_point <- t( rbind( runif(n_data)*(max(c_point[,1])-min(c_point[,1]))+min(c_point[,1]), 
                     runif(n_data)*(max(c_point[,2])-min(c_point[,2]))+min(c_point[,2]), 
                     runif(n_data)*(max(c_point[,3])-min(c_point[,3]))+min(c_point[,3])) )


# COMPUTATION -------------------------
# generates the interpolation
R1 <- interpolate_bspline(position, c_point, degree = deg, is.periodic = TRUE)

# approximated parameters
ap_u <- bspline_footpoint(d_point, c_point, 
                          degree=deg, step=1, tolerance=1e-9, verbose=TRUE, algorithm="directional",
                          initial.scan.points=50,
                          max.iterations = 5000, is.periodic = TRUE, native = FALSE) #ap_u
Rd <- interpolate_bspline(ap_u, c_point, degree = deg, is.periodic=TRUE)

# TEST PROJECTION QUALITY -----------------------
# perpendicularity at the footpoint
proj_angles <- projection_angles(ap_u, d_point, c_point, degree = deg, 
                                 is.periodic = TRUE)
print("projection angles ***")
print(proj_angles)


# DRAWING -----
require(rgl)
draw_projs <- function(x, v){
  
  points3d(x, pch=2, col=rgb(0,1,0,0), size = 6)
  
  points3d(v, pch=2, col=rgb(0,1,1,0), size = 6)
  
  for(i in 1:nrow(x))
  lines3d(rbind(x[i,],v[i,]), col=rgb(0,1,1,0), lty=3, lwd=.5)
  
}



plot3d(R1, type="l", 
       aspect=FALSE, xlab ="x", ylab="y", zlab ="z", 
       col=rgb(0,0,.3), lwd=2, main = "cubic spline [periodic]" )

par3d(windowRect=c(184, 130, 813, 666))
#par3d(zoom = 1)

#draws the control points 
points3d(c_point, pch=2, col=rgb(1,0.5,0.5), size = 3)
lines3d(c_point, col=rgb(1,0,0), lty=2, lwd=.5)

draw_projs(d_point, Rd)


